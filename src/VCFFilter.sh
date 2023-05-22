#!/bin/bash

selected_individuals="Data/individuals.tsv"
selected_individuals_txt="Data/individuals.txt"

vcf_file="/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/malawi_cichlids_v3_phase.biallelic_snps.pass.ancestral_as_sample.benthic.chr1.vcf.gz"
bed_file="/scratch/antwerpen/208/vsc20886/CompBioProject_2324/Data/chr1_windows.bed"
output_file="${VSC_DATA}/filtered.vcf.gz"
finaloutput_file="${VSC_DATA}/finaloutput.vcf.gz"

echo "filter VCF file using BED file"
# Filter VCF using BED file and redirect output to bgzip-compressed VCF file
bedtools intersect -a "$vcf_file" -b "$bed_file" -header | bgzip > "$output_file"
echo "VCF file filtered using BED file and saved as compressed VCF"

echo "convert TSV file with selected individuals to TXT format"
awk -F'\t' 'NR>1 {print $2; print $3}' "$selected_individuals" > "$selected_individuals_txt"
echo "TSV file with selected individuals converted to TXT format"

echo "filter VCF file further based on selected individuals"
# Filter VCF further based on selected individualsbcftools view -S "$selected_individuals_txt" "$output_file" | gzip "$finaloutput_file"

echo "VCF file filtered further based on selected individuals"
echo "Final output VCF file compressed using bgzip"
