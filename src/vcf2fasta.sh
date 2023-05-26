#!/bin/bash

reference_genome="/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"
compressed_vcf_file="/data/antwerpen/208/vsc20886/finaloutput.vcf.gz"
gtf_file="/scratch/antwerpen/208/vsc20886/CompBioProject_2324/Data/GCF_900246225.1_fAstCal1.2_genomic_FIXED_chromnames.gtf"
#output_file="/scratch/antwerpen/208/vsc20886/CompBioProject_2324/GCF_900246225.1_fAstCal1.2_genomic_FIXED_chromnames.gtf"

# Index the reference
samtools faidx ${reference_genome}

# And the VCF file should be tabix indexed and compressed:
#bgzip my_vcf_file.vcf
tabix "${compressed_vcf_file}"

# Convert the VCF file to FASTA
python3 /scratch/antwerpen/208/vsc20886/CompBioProject_2324/libs/vcf2fasta.py -f ${reference_genome} -v ${compressed_vcf_file}.tbi -g ${gtf_file} -e gene


