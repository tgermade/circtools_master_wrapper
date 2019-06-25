#!/bin/bash

# @Author: Tomas Germade <tgermade>
# @Date:   Monday, June 24, 2019 12:08
# @Email:  tgermade@student-net.ethz.ch
# @Project: ETH Zuerich semester project Laboratory of Systems Neuroscience
# @Last modified by:   tgermade
# @Last modified time: Monday, June 24, 2019 12:08


# $1 -> work directory
# $2 -> define which parts to run (all, trim, star, dcc, fuchs)
# $3 -> read directory (fastq.gz files)
# $4 -> organism (mouse, rat)
# $5 -> define if read data is paired or not (paired, unpaired)
# $6 -> maximal number of threads that should be used (default, [num])
# $7 -> batch name (name for the collection of samples to run)

# One wrapper to rule them all and in a bash script bind them.
############################################################################################
# Usage
#######

usage() {
    cat <<EOF

Usage:
	$0 [work dir] [parts to run] [read files dir] [organism] [form of read data] [thread number] [batch name]

Options:
	- [parts to run]    	can assume values: 'all', 'star', 'dcc', 'fuchs'
	- [organism]        	can assume values: 'mouse', 'rat', 'rat_assembled'
	- [form of read data]	can assume values: 'paired', 'unpaired'
	- [thread number]		if value is 'default', will use max. number of cores - 8,
					otherwise input any number of cores you want to use (integer)
	- [batch name]		input unique name for the overall read collection;
					if no name is given, it will be automatically generated
					based on the supplied [read files dir]
					(only affects DCC module output)

Examples:
	~/script.sh ~/testing/ fuchs /mnt/schratt/mouse_RNAseR/raw mouse paired default "mouse_GRCm38_RNAseR"
	~/script.sh ~/testing dcc /mnt/schratt/mouse_control_hippocampus/ mouse unpaired 12 ""

EOF
}

if [ ! $# = 6 ] && [ ! $# = 7 ]; then
    usage
    exit 1
fi

############################################################################################
# Allocation
############

if [ $4 = "mouse" ]; then
	star_index="/reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes_STARIndex"
	genome_reference="/reference/Mus_musculus/GENCODE/GRCm38.p5/Sequence/WholeGenomeFasta/genome.fa"
	gene_annotation="/reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes.gtf"
	exon_annotation="/mnt/schratt/tgermade_test/mouse_GRCm38_p5.GENCODE.exons.bed"
elif [ $4 = "rat" ] || [ $4 = "rat_assembled"]; then
	star_index="/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/shortRNA/starINDEX"
	genome_reference="/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/WholeGenomeFasta/genome.fa"
	if [ $4 = "rat" ]; then
		gene_annotation="/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.gtf"
		exon_annotation="/mnt/schratt/tgermade_test/rat_Rnor_6.0.Ensembl.exons.bed"
	elif [ $4 = "rat_assembled" ]; then
		gene_annotation="/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes_assembled.gtf"
		exon_annotation="/mnt/schratt/tgermade_test/rat_Rnor_6.0.Ensembl.exons_assembled.bed"
	fi
fi

# get full script path name
script_path=$(dirname `realpath $0`)

# create work directory if not yet existing
mkdir -pv $1 &&
# get full work directory path name
work_dir=`realpath $1` &&

# get full read file path name
cd $3 &&
read_dir=`pwd` &&
# define read markers for naming convention
if [ $5 = "paired" ]; then
	if [[ `echo ./*` = *_R1.fastq.gz* ]]; then
		r1_marker=_R1 &&
		r2_marker=_R2
	else
		echo "INPUT DEMAND: Input read 1 marker of raw read data files. Examples: '_R1', '_read1', etc." &&
		read r1_marker &&
		r2_marker=`echo $r1_marker | tr 1 2`
	fi
fi

# setup temporary directory
mkdir -pv $work_dir/.tmp &&

# specify set of read files to analyze (prevent computation of sets already computed)
if [ $5 = "paired" ]; then
	for file in ./*.fastq.gz; do
		if [[ $file = *$r1_marker* ]]; then
	            	echo $(basename $file) | sed "s%\.fastq\.gz$%%;s%$r1_marker%%"
		fi
	done > $work_dir/.tmp/names.tmp

elif [ $5 = "unpaired" ]; then
	for file in ./*.fastq.gz; do
		echo $(basename $file) | sed "s%\.fastq\.gz$%%"
	done > $work_dir/.tmp/names.tmp
fi

# setup maximum threads to use
if [ $6 = "default" ]; then
	procs=$(expr `nproc` - 8)
else
	procs=$6
fi

# create batch name (name of directory which will contain information for all read samples)
if [ -z $7 ]; then
	if [ `echo $(basename $read_dir)` = "raw" ]; then
      	batch_name=`echo $(dirname $read_dir) | rev | cut -d "/" -f1 | rev`
	else
      	batch_name=`echo $(basename $read_dir)`
	fi
else
	batch_name=$7
fi

############################################################################################
# Tools
#######

# ALIGNMENT

if [ $2 = "star" ] || [ $2 = "all" ]; then

	# setup alignment output directory
	mkdir -pv $work_dir/star &&
	cd $read_dir &&

	# different wrapper options depending on paired / unpaired read data
	if [ $5 = "paired" ]; then
		align_option=$(echo "{1} {2} $r1_marker ::: *$r1_marker* ::: *$r2_marker*")
      elif [ $5 = "unpaired" ]; then
		align_option=$(echo "{} ::: *")
	fi

	# run STAR
	# PAIRED: parallel slurm_circtools_detect_mapping [STAR index] [target dir] [gene annotation GTF file] [thread number] [Read 1 file] [Read 2 file] [Read 1 marker, e.g. R1]
      # UNPAIRED: parallel slurm_circtools_detect_mapping [STAR index] [target dir] [gene annotation GTF file] [thread number] [Read file]
	echo "Running STAR."
	parallel -j1 --xapply $script_path/slurm_circtools_detect_mapping.sh $star_index $work_dir/star/ $gene_annotation $procs $align_option
fi


# DETECTION

if [ $2 = "dcc" ] || [ $2 = "all" ]; then
	cd $work_dir &&
	# create directories
	mkdir -pv circtools &&
	mkdir -pv circtools/01_detect &&

	# setup links to /star folder
	echo "Setting up links."
	cd $work_dir/star &&
	## generate batch directory
	mkdir -pv ../circtools/01_detect/$batch_name
	## link aligned reads (bam files) and indexing files for the aligned reads
	parallel --plus ln -s `pwd`/{}/Aligned.noS.bam ../circtools/01_detect/$batch_name/{}.bam :::: ../.tmp/names.tmp
	parallel --plus ln -s `pwd`/{}/Aligned.noS.bam.bai ../circtools/01_detect/$batch_name/{}.bam.bai :::: ../.tmp/names.tmp
	## link chimeric junction files of the main mapping
	parallel --plus ln -s `pwd`/{}/Chimeric.out.junction ../circtools/01_detect/$batch_name/{}.Chimeric.out.junction :::: ../.tmp/names.tmp
	## mate links for paired reads
	if [ $5 = "paired" ]; then
		parallel --plus ln -s `pwd`/{}/mate1/Chimeric.out.junction ../circtools/01_detect/$batch_name/{}.mate1.Chimeric.out.junction :::: ../.tmp/names.tmp
		parallel --plus ln -s `pwd`/{}/mate2/Chimeric.out.junction ../circtools/01_detect/$batch_name/{}.mate2.Chimeric.out.junction :::: ../.tmp/names.tmp
		## (important for circtools reconstruct (FUCHS))
		parallel --plus ln -s `pwd`/{}/mate1/Aligned.noS.bam ../circtools/01_detect/$batch_name/{}.mate1.bam :::: ../.tmp/names.tmp
		parallel --plus ln -s `pwd`/{}/mate2/Aligned.noS.bam ../circtools/01_detect/$batch_name/{}.mate2.bam :::: ../.tmp/names.tmp
	fi

	cd ../circtools/01_detect/$batch_name

	# create txt files
	if [ $5 = "paired" ]; then
			ls *bam | grep -v mate > bam_files.txt &&
			ls *junction | grep  mate1 > mate1 &&
			ls *junction | grep  mate2 > mate2 &&
			ls *junction | grep -v mate > samplesheet

	elif [ $5 = "unpaired" ]; then
			ls *bam | grep -v mate > bam_files.txt &&
			ls *junction | grep -v mate > samplesheet
	fi


	# check if data is unstranded, first-stranded or second-stranded
	echo "Checking sample strandedness (.85 cutoff):"

	decideStrandedness () {

		if (( $(echo "$1 >= 0.85" | bc -l) )); then
			echo "firststrand";
		elif (( $(echo "$1 <= 0.15" | bc -l) )); then
			echo "secondstrand";
		else
			echo "unstranded";
		fi
	}

	## for each sample, calculate ratio btw. first strand alignments and total alignments
	sarr=()
	parr=()

	for dir in `cat $work_dir/.tmp/names.tmp`; do
            cd $work_dir/star/$dir
		p=`awk '{if(NR >= 5) { sumF+=$3; sumS+=$4 }} END{print sumF/(sumF+sumS)}' ReadsPerGene.out.tab;`
		strandcall=`decideStrandedness $p`
		parr+=($p)
		sarr+=($strandcall)
		echo "firststrand ratio: " $f $p "-->" $strandcall
		cd ..
	done

	# check if all samples of batch share the same strandedness
	if [ `printf '%s\n' "${sarr[@]}" | sort | uniq | wc -l` -gt 1 ]; then
		echo "WARNING: not all libraries appear to have the same strand specificity!"
		# check if (rough) median ratio passes a threshold
		n=${#parr[@]}
		strandcall=decideStrandedness `printf '%s\n' "${parr[@]}" | sort -n | tail -n $(( $n / 2 + 1 )) | head -1`
		echo "Falling back on $strandcall mode, but consider splitting the libraries."
	else
		# all calls are identical, take the first
		strandcall=${sarr[0]}
	fi

	echo "Running DCC on $strandcall mode."


	# allocate DCC options
	cd $work_dir/circtools/01_detect/$batch_name

	if [ $strandcall = "unstranded" ]; then
		strand_option="-N"
	elif [ $strandcall = "secondstrand" ]; then
		strand_option="-ss"
	else strand_otpion=""
	fi

	if [ $5 = "paired" ]; then
		paired_option=$(echo "-mt1 @mate1 -mt2 @mate2 -Pi \\")
	else paired_option=""
	fi

	# run DCC
	source activate python27

	circtools detect 	@samplesheet \
			$paired_option
	            -D \
	            -an $gene_annotation \
	            -F \
	            -Nr 2 1 \
	            -fg \
	            -G \
	            -A $genome_reference \
	            -O ./output \
	            -T $procs \
	            -B @bam_files.txt \
			$strand_option

	conda deactivate
fi


# RECONSTRUCTION

if [ $2 = "fuchs" ] || [ $2 = "all" ]; then

	mkdir -pv $work_dir/circtools
	mkdir -pv $work_dir/circtools/03_reconstruct

	cd $work_dir/circtools/01_detect


############################################################################################

	circtools_reconstruct () {

		source activate python27

		# allocation
		sample_name=$1
	      main_out=$2/
		bed_file=$3
		dcc_dir=$4
		dcc_out_dir=$5
		format=$6
		procs=$7
		tmp_folder=/tmp/global_tmp/

		main_bam=$dcc_dir/${sample_name}.bam &&
		main_junction=$dcc_dir/${sample_name}.Chimeric.out.junction &&

		mkdir -p $main_out/${sample_name}

		if [ $format = "paired" ]; then
		      mate1_bam=$dcc_dir/${sample_name}.mate1.bam &&
		      mate1_junction=$dcc_dir/${sample_name}.mate1.Chimeric.out.junction &&

		      mate2_bam=$dcc_dir/${sample_name}.mate2.bam &&
		      mate2_junction=$dcc_dir/${sample_name}.mate2.Chimeric.out.junction.fixed &&

		      merged_bam=$main_out/${sample_name}/${sample_name}_merged.bam
		fi

		# merge both mate BAM files into one new BAM file
		if [ $format = "paired" ]; then
		      samtools merge -l 9 -@ 8 $merged_bam $main_bam $mate1_bam $mate2_bam &&
		      # re-index the newly aggregated BAM file
		      samtools index $merged_bam &&
		      # run FUCHS
		      FUCHS -N $sample_name -D $dcc_out_dir/CircRNACount -B $merged_bam -A $bed_file -O $main_out/${sample_name} -F $mate1_junction -R $mate2_junction -J $main_junction -T $tmp_folder -p ensembl -r 2 -e 1 -q 2 -P $procs

		elif [ $format = "unpaired" ]; then
		      # run FUCHS
		      FUCHS -N $sample_name -D $dcc_out_dir/CircRNACount -B $main_bam -A $bed_file -O $main_out/${sample_name} -J $main_junction -T $tmp_folder -p ensembl -r 2 -e 1 -q 2 -P $procs
		fi

		conda deactivate
	}

############################################################################################

	export -f circtools_reconstruct

	# Parallel reconstruction:
	# parallel slurm_circtools_reconstruct.sh [Sample name] [target dir] [exon annotation BED file] [path to DCC dir] [CircRNACount dir] [format: paired reads, unpaired reads] [thread number]
	echo "Running FUCHS."
	parallel -j1 --xapply circtools_reconstruct {} ../03_reconstruct $exon_annotation ./$batch_name ./$batch_name/output $5 $procs :::: ../../.tmp/names.tmp

fi
