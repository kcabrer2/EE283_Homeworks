###Organizing anf converting DNA-seq files to Illumina format###
rename _1.fq.gz _F.fq.gz ADL*
rename _2.fq.gz _R.fq.gz ADL*
rename ADL06 A4 ADL06*
rename ADL09 A5 ADL09*
rename ADL10 A6 ADL10*
rename ADL14 A7 ADL14*

#!/bin/bash
#$ -N IlluminaConvert
#$ -q bio,abio
#$ -pe openmp 1
#$ -ckpt restart
#$ -R y
#$ -t 1-24

names="/pub/kcabrer2/EE283/Bioinformatics_Course/DNAseq/DNAseq_Samples.txt"

files=`head -n $SGE_TASK_ID ${names} | tail -n 1`

gunzip -c ${files} | /data/apps/user_contributed_software/kcabrer2/seqtk/seqtk seq -Q64 -V - | gzip -c - > IlluminaDNAseq_${files}

ls IlluminaDNAseq*_F.fq.gz | sed 's/_F.fq.gz//' > DNAseqIllumina.prefixes.txt

########################################
## Align your samples to a ref genome ##
########################################
#!/bin/bash
#$ -N DNAseq_BT2align
#$ -q bio,abio
#$ -pe openmp 8
#$ -ckpt restart
#$ -R y
#$ -t 1-12

module load bowtie2/2.2.7  
module load samtools/1.3

path="/pub/kcabrer2/EE283/Bioinformatics_Course/DNAseq"
ref="/pub/kcabrer2/EE283/Bioinformatics_Course/ref/dmel-all-chromosome-r6.13.fasta"
files="/pub/kcabrer2/EE283/Bioinformatics_Course/DNAseq/DNAseqIllumina.prefixes.txt"
bt2index="/pub/kcabrer2/EE283/Bioinformatics_Course/ref/dmel-all-chromosome-r6.13.index"
output="/pub/kcabrer2/EE283/Bioinformatics_Course/DNAseq/DNAseq_Outputs"
bam="/pub/kcabrer2/EE283/Bioinformatics_Course/DNAseq/DNAseq_Outputs/BAM_Files"
qc="/pub/kcabrer2/EE283/Bioinformatics_Course/DNAseq/DNAseq_Outputs/QC"

prefix=`head -n $SGE_TASK_ID ${files} | tail -n 1`
mkdir ${path}/DNAseq_Outputs
mkdir ${output}/BAM_Files
mkdir ${output}/QC

#Align with Bowtie2 & convert to BAM#
bowtie2 -p 8 -x ${bt2index} -1 ${path}/${prefix}_F.fq.gz -2 ${path}/${prefix}_R.fq.gz | samtools view -bS - > ${bam}/${prefix}.bam
#Sort your BAM file#
samtools sort ${bam}/${prefix}.bam -o ${bam}/${prefix}.sort.bam
##Remove duplicate reads 
java -Xmx20g -jar /data/apps/picard-tools/1.87/MarkDuplicates.jar M=${qc}/${prefix}.nodups.matrix.txt INPUT=${bam}/${prefix}.sort.bam OUTPUT=${bam}/${prefix}.nodups.bam REMOVE_DUPLICATES=true  
#Adds read groups to your data so you can do downstream stuff with your data in GATK#
java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=${bam}/${prefix}.nodups.bam O=${bam}/${prefix}.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=${prefix} RGSM=${prefix} VALIDATION_STRINGENCY=LENIENT
#Indexes your read group-modded BAM files#
samtools index ${bam}/${prefix}.RG.bam




bam="/pub/kcabrer2/EE283/Bioinformatics_Course/DNAseq/DNAseq_Outputs/BAM_Files"

java -Xmx20g -jar /data/apps/picard-tools/1.87/MergeSamFiles.jar $(printf 'I=%s ' ${bam}/*.RG.bam) SO=coordinate AS=true VALIDATION_STRINGENCY=SILENT O=${bam}/mergedbams.bam
