#!/bin/bash
ref="/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"
vcf="/data/antwerpen/208/vsc20886/finaloutput.vcf.gz"
cut_ref="./Data/genome_chr1.fa"
head -c 100M  ${ref} | sed -n '/>chr1/,/chr2/p' | head -n-1 > ${cut_ref}

temp_dir=`mktemp -d`
if [ -z "$temp_dir" ]
then
  echo "Cannot create temp dir!"
  exit 1;
fi

#if [ -z "$1" ]
#then
#  echo "Need an argument <window_number>"
#  exit 1;
#fi


# !!! doesn't affect rehydration part (TODO: fix)
window_size=100000

window_n=${1}
window_start=$((window_size*window_n))
window_end=$((window_start + window_size))
window_start_1_based=$((window_start + 1))

window_VCF="${temp_dir}/window.chr1.${window_start}.vcf.gz"
vcf2phylip_prefix="${temp_dir}/window.chr1.${window_start}.aln"

gatk SelectVariants -R ${ref} -V ${vcf} -L chr1:${window_start_1_based}-${window_end} -O ${window_VCF}
python3 ../scripts/vcf2phylip.py -i ${window_VCF} --output-prefix ${vcf2phylip_prefix} -f -w
python3 ../scripts/rehydrate_vcf_to_fas.py ${vcf2phylip_prefix}.min4.fasta ${vcf2phylip_prefix}.min4.used_sites.tsv ${cut_ref} ${window_start}

cp "${vcf2phylip_prefix}.min4.rehydrated.fasta" .

echo "About to rm ${temp_dir}"
#rm -rf "${temp_dir}"