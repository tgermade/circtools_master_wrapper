#!/bin/bash

# @Author: Tomas Germade <tgermade>
# @Date:   Monday, June 3, 2019 11:57
# @Email:  tgermade@student-net.ethz.ch
# @Project: ETH Zuerich semester project Laboratory of Systems Neuroscience 
# @Last modified by:   tgermade
# @Last modified time: Monday, June 3, 2019 11:57


# $1 -> work directory
# $2 -> define which parts to run (all, trim, star, dcc, fuchs)
# $3 -> read directory (fastq.gz files)
# $4 -> organism (mouse, rat)
# $5 -> define if read data is paired or not (paired, unpaired)
# $6 -> maximal number of threads that should be used (default, [num])

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
					(only affects DCC module!)

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
		# parallel slurm_circtools_detect_mapping [STAR index] [target dir] [gene annotation GTF file] [thread number] [Read 1 file] [Read 2 file] [Read 1 marker, e.g. R1]
		parallel -j1 --xapply ~/slurm_circtools_detect_mapping.sh $star_index $work_dir/star/ $gene_annotation $procs {1} {2} $r1_marker ::: *$r1_marker* ::: *$r2_marker*

	elif [ $5 = "unpaired" ]; then
		# parallel slurm_circtools_detect_mapping [STAR index] [target dir] [gene annotation GTF file] [thread number] [Read file]
		parallel -j1 --xapply ~/slurm_circtools_detect_mapping.sh $star_index $work_dir/star/ $gene_annotation $procs {} ::: *
	fi
fi


# DETECTION

if [ $2 = "dcc" ] || [ $2 = "all" ]; then
	cd $work_dir &&
	# create directories
	mkdir -pv circtools &&
	mkdir -pv circtools/01_detect &&

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

	# setup links to /star folder
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
	for dir in `cat $work_dir/.tmp/names.tmp`; do
		cd $work_dir/star/$dir
		cat ReadsPerGene.out.tab | awk 'NR >= 5 { print }' | awk '{sum+=$2} END{print sum}'
		cat ReadsPerGene.out.tab | awk 'NR >= 5 { print }' | awk '{sum+=$3} END{print sum}'
		cat ReadsPerGene.out.tab | awk 'NR >= 5 { print }' | awk '{sum+=$4} END{print sum}'

	done
	#cd ..
	# Parallel detection:
	# parallel slurm_circtools_detect.sh [sample names] [paired/unpaired data] [# of cores]
	#parallel -j1 --xapply ~/slurm_circtools_detect.sh {} $5 $procs $gene_annotation $genome_reference :::: $work_dir/.tmp/names.tmp
	~/slurm_circtools_detect.sh "" $5 $procs $gene_annotation $genome_reference
fi


# RECONSTRUCTION

if [ $2 = "fuchs" ] || [ $2 = "all" ]; then
	mkdir -pv circtools
	mkdir -pv circtools/02_reconstruct
	# download the wrapper scrips for the reconstruct module (choose the reconstruct wrapper in the home dir! it's adapted to be unpaired compatible)
	#wget -nc https://raw.githubusercontent.com/dieterich-lab/bioinfo-scripts/master/slurm_circtools_reconstruct.sh
	# add execute permission
	#chmod 755 slurm_circtools_reconstruct.sh

	cd $work_dir/circtools/01_detect

	# Parallel reconstruction:
	# parallel slurm_circtools_reconstruct.sh [Sample name] [target dir] [exon annotation BED file] [path to DCC dir] [CircRNACount dir]
	parallel -j5 --xapply  ~/slurm_circtools_reconstruct.sh {} ../03_reconstruct $exon_annotation ./{} ./{}/output $5 3 :::: ../../.tmp/names.tmp
fi
