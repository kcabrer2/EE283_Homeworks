###########################
#### ATAC-seq Analysis ####
###########################
##########################################
## Make a list of all you file prefixes ##
##########################################

### Code to change ATAC-seq filenames. It creates a symlink where the name = ATAC_barcode_ForR.fq.gz ###
for i in $(ls *.gz); do echo ATAC_$(echo $i | cut -f5,6 -d"_") | xargs ln -s $i; done
ls ATAC_*_R1.fq.gz | sed 's/_R1.fq.gz//' > ATAC.prefixes.txt

########################################
## Align your samples to a ref genome ##
########################################
#!/bin/bash
#$ -N ATAC_BT2align
#$ -q bio,abio
#$ -pe openmp 8
#$ -ckpt restart
#$ -R y
#$ -t 1-24

module load bowtie2/2.2.7  
module load samtools/1.3

ref="/pub/kcabrer2/EE283/Bioinformatics_Course/ref/dmel-all-chromosome-r6.13.fasta"
output="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs"
files="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC.prefixes.txt"
path="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq"
bt2index="/pub/kcabrer2/EE283/Bioinformatics_Course/ref/dmel-all-chromosome-r6.13.index"

prefix=`head -n $SGE_TASK_ID ${files} | tail -n 1`
mkdir ${output}/BAM_files
mkdir ${output}/QC

#Align with Bowtie2 & convert to BAM#
bowtie2 --maxins 2000 -p 8 -x ${bt2index} -1 ${path}/${prefix}_R1.fq.gz -2 ${path}/${prefix}_R2.fq.gz | samtools view -bS - > ${output}/BAM_files/${prefix}.bam
#Sort your BAM file#
samtools sort ${output}/BAM_files/${prefix}.bam -o ${output}/BAM_files/${prefix}.sort.bam
##Remove duplicate reads 
java -Xmx20g -jar /data/apps/picard-tools/1.87/MarkDuplicates.jar M=${output}/QC/${prefix}.nodups.matrix.txt INPUT=${output}/BAM_files/${prefix}.sort.bam OUTPUT=${output}/BAM_files/${prefix}.nodups.bam REMOVE_DUPLICATES=true  
#Adds read groups to your data so you can do downstream stuff with your data in GATK#
java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=${output}/BAM_files/${prefix}.nodups.bam O=${output}/BAM_files/${prefix}.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=${prefix} RGSM=${prefix} VALIDATION_STRINGENCY=LENIENT
#Indexes your read group-modded BAM files#
samtools index ${output}/BAM_files/${prefix}.RG.bam

############################
##Get a TSS file from UCSC##
############################
- Click on the Tools tab & choose Table Browser
- Choose the following Parameters:
	- Clade: 'Insect'
	- Genome: 'D. melanogaster'
	- Assembly: 'Aug. 2014 (BDGP Release 6 + ISO1 MT/dm6)'
	- Group: 'Genes and Gene Predictions'
	- Track: 'NCBI RefSeq'
	- Table: 'UCSC RefSeq (refGene)'
	- Region: 'Genome'
	- Output Format: 'BED'
- Click 'get output'
- Choose the following parameter for Create one BED record per
	- Upstream by '1 bases'

#Removing duplicates and sorting the TSS BED file
#Awk command says input and output are tab-delimited and returns columns 1,2,3,6 (chrom, start, end, strand)
#The -u in the sort command returns the first unique iteration of a genomic position
cat dm6_tss.bed|awk '{OFS="\t";FS="\t"}{print $1, $2, $3, $6}'|sort -u > dm6_nodups_sorted_tss.bed

#For the QC stuff you're going to have to remove the 'chr' from the chromosome name in the TSS BED file. I modded a script to do just that
python Remove_chr.py --f1 dm6_nodups_sorted_tss.bed --f2 dm6_nodups_sorted_nochr_tss.bed

######################################################
##Run QC on your reads & Get ready for Visualization##
######################################################

#Download the pyMakeVplot.py script from https://github.com/jinxu9/ATACseq/blob/master/libs/pyMakeVplot.py 
#Variables
	# -a is your input BAM file
	# -b is your sorted, no dups, no chr TSS BED file
	# -o is your output, 
	# -c is the number of threads


#!/bin/bash
#$ -N ATAC_QC
#$ -q bio,abio
#$ -pe openmp 1
#$ -ckpt restart
#$ -R y
#$ -t 1-24

module load enthought_python/7.3.2
module load samtools/1.3
module load bedtools/2.25.0

code="/pub/kcabrer2/EE283/Bioinformatics_Course/Code"
bam_loc="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs/BAM_files"
files="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC.prefixes.txt"
tss_bed="/pub/kcabrer2/EE283/Bioinformatics_Course/ref/dm6_nodups_sorted_nochr_tss.bed"
output="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs"
output_qc="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs/QC"
bed_loc="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs/BED_files"

prefix=`head -n $SGE_TASK_ID ${files} | tail -n 1`
mkdir ${output_qc}/tssplots
mkdir ${output}/BED_files

#Collect TSS-proximal reads
python ${code}/pyMakeVplot_mod.py -a ${bam_loc}/${prefix}.RG.bam -b ${tss_bed} -e 2000 -p ends -v -u -o ${output_qc}/tssplots/${prefix}_tssplot -c 1

##getting q30 reads  
samtools view -f2 -q30 -b ${bam_loc}/${prefix}.RG.bam > ${bam_loc}/${prefix}.q30.nodups.bam  

##indexing q30 bam  
samtools index ${bam_loc}/${prefix}.q30.nodups.bam

##converting bam to bed files  
bedtools bamtobed -i ${bam_loc}/${prefix}.q30.nodups.bam |gzip -c > ${bed_loc}/${prefix}.bed.gz  

##converting to just the cut sites  
gzip -dc  ${bed_loc}/${prefix}.bed.gz | awk 'BEGIN{OFS="\t";FS="\t"};{if($6=="+"){print $1, $2+1, $2+1, $4, $5, $6} else {print $1, $3, $3, $4, $5, $6}}' > ${bed_loc}/${prefix}.cuts.bed  
gzip ${bed_loc}/${prefix}.cuts.bed  

##removing mitochondrial reads  
gzip -dc ${bed_loc}/${prefix}.cuts.bed.gz | awk 'BEGIN{OFS="\t";FS="\t"};{if($1!="chrM"){print $1, $2, $3, $4, $5, $6}}' > ${bed_loc}/${prefix}.cuts.noM.temp.bed   

##Sorting  
sort -k 1,1 ${bed_loc}/${prefix}.cuts.noM.temp.bed | gzip -c > ${bed_loc}/${prefix}.sorted.bed.gz  

###############################################
## Calculate Coverage Per Sample & Call Peaks##
###############################################

#!/bin/bash
#$ -N ATAC_CovNPeak
#$ -q bio,abio
#$ -pe openmp 1
#$ -ckpt restart
#$ -R y
#$ -t 1-24

module load samtools/1.3
module load bedtools/2.25.0
module load enthought_python/7.3.2
module load homer/4.10

files="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC.prefixes.txt"
code="/pub/kcabrer2/EE283/Bioinformatics_Course/Code"
ref="/pub/kcabrer2/EE283/Bioinformatics_Course/ref/dmel-all-chromosome-r6.13.fasta"
sizes="/pub/kcabrer2/EE283/Bioinformatics_Course/ref/dm6.chrom.sizes"
bed_loc="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs/BED_files"
output="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs"
bw_loc="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs/BigWigs"


prefix=`head -n $SGE_TASK_ID ${files} | tail -n 1`
mkdir ${output}/BigWigs
mkdir ${output}/peakcalls

gzip -dc ${bed_loc}/${prefix}.sorted.bed.gz > ${bed_loc}/${prefix}.sorted.bed

##This command is powerful since it not only adds the string 'chr' to the beginning of your chromosome numbers, but also removes mitochondrial reads and reads mapping to non-canonical chromosomes
python ${code}/Chromosome_Rename.py --f1 ${bed_loc}/${prefix}.sorted.bed --f2 ${bed_loc}/${prefix}.chr.sorted.bed

#Calculate bedgraph coverage of 5’ positions of your BED files genome-wide
bedtools genomecov -bg -5 -i ${bed_loc}/${prefix}.chr.sorted.bed -g ${sizes} > ${bed_loc}/${prefix}.bedgraph  

bedGraphToBigWig ${bed_loc}/${prefix}.bedgraph ${sizes} ${bw_loc}/${prefix}.q30.bw  


###Peak Calling using HOMER###

##converting bam to sorted sam  
samtools view ${output}/BAM_files/${prefix}.RG.bam -o ${output}/BAM_files/${prefix}.RG.sam  

##making tag directory  
makeTagDirectory ${output}/BAM_files/${prefix}.homer ${output}/BAM_files/${prefix}.RG.sam -format sam  
touch ${output}/BAM_files/${prefix}.tagsdone 

##calling peaks with homer  
findPeaks ${output}/BAM_files/${prefix}.homer -o ${output}/peakcalls/${prefix}.homer.peaks -style dnase  

grep -v "^#" ${output}/peakcalls/${prefix}.homer.peaks | awk 'BEGIN{OFS="\t";FS="\t"};{print $2, $3, $4, $1, $8}' > ${output}/peakcalls/${prefix}.homer.bed 



##########################################################################################
##########################################################################################
##########################################################################################

##The following was used to make BigWig files of the entire coverage, not just the cut sites

#!/bin/bash
#$ -N ATAC_test
#$ -q bio,abio
#$ -pe openmp 1
#$ -ckpt restart
#$ -R y
#$ -t 1-24

module load samtools/1.3
module load bedtools/2.25.0
module load enthought_python/7.3.2

files="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC.prefixes.txt"
code="/pub/kcabrer2/EE283/Bioinformatics_Course/Code"
sizes="/pub/kcabrer2/EE283/Bioinformatics_Course/ref/dm6.chrom.sizes"
bed_loc="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs/BED_files"
output="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs"
test_loc="/pub/kcabrer2/EE283/Bioinformatics_Course/ATACseq/ATAC_Outputs/test"


prefix=`head -n $SGE_TASK_ID ${files} | tail -n 1`
mkdir ${output}/test


gzip -dc  ${bed_loc}/${prefix}.bed.gz | sort -k 1,1 - >  ${test_loc}/${prefix}.test.sort.bed

python ${code}/Chromosome_Rename.py --f1 ${test_loc}/${prefix}.test.sort.bed --f2 ${test_loc}/${prefix}.test.chr.sort.bed

bedtools genomecov -bg -i ${test_loc}/${prefix}.test.chr.sort.bed -g ${sizes} > ${test_loc}/${prefix}.test.bedgraph 

bedGraphToBigWig ${test_loc}/${prefix}.test.bedgraph ${sizes} ${test_loc}/${prefix}.test.bw