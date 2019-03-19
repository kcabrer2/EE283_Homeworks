Project Description

The plan for the project was to call peaks on the ATAC-seq founder line data and overlap them with non-synonymous SNPs identified in those lines from the DNA-seq data.

Due to my project in lab being really labor intensive lately I got to complete up to the peak-calling part of the ATAC-seq pipeline and only part of the DNA-seq pipeline.

For ATAC-seq I was able to align the reads to the genome, do some QC by visualizing how many reads fall close to TSSs (results of which are included in the TSS plots), make tracks of just the cut sites and of whole coverage, and ultimately call peaks with Homer.

For DNA-seq I was able to align to the genome and do part of the SNP calling pipeline. I did ultimately get a VCF file of SNPs from Tony but didn't get to do the whole gamut of finding non-synonymous SNPs

My plan once I was able to call SNPs using GATK in DNA-seq was to use the interectBed function from Bedtools to overlap my ATAC-seq peaks BED file with the VCF file (I didn't have time to finish this)