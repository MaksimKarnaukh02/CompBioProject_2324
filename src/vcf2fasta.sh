vcf_file=".vcf.gz"
fasta_file="/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"

#cat "$fasta_file" | vcf-consensus "$vcf_file" > output.fa

bcftools consensus -f "$fasta_file" "$vcf_file" > out.fa