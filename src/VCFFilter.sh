#!/bin/bash

selected_individuals="selected_individuals.txt"

vcf_file="malawi_cichlids_v3_phase.biallelic_snps.pass.ancestral_as_sample.benthic.chr1.vcf.gz"
bed_file="goodwindows.bed"
output_file="filtered.vcf"

# Filter VCF using BED file
bedtools intersect -a "$vcf_file" -b "$bed_file" -header > "$output_file"

# Filter VCF further based on selected individuals
bcftools view -S "$selected_individuals" "$output_file" -o finaloutput.vcf