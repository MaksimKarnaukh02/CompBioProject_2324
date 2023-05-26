import subprocess
from pprint import pprint

# import sys, os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
from pyfaidx import FastaVariant
import vcf


def fasta_alignment_from_vcf(vcf_file, ref):
    """
    Get a fasta alignment for all snp sites in a multi
    sample vcf file, including the reference sequence.
    """
    print('reading vcf', flush=True)
    # index vcf
    cmd = 'tabix -p vcf -f {i}'.format(i=vcf_file)
    tmp = subprocess.check_output(cmd, shell=True)
    # get samples from vcf
    vcf_reader = vcf.Reader(filename=vcf_file)
    samples = vcf_reader.samples
    print('%s samples' % len(samples), flush=True)
    result = []

    # reference sequence
    reference = Fasta(ref)
    chrom = list(reference.keys())[0]

    # get the set of all sites first
    sites = []
    print('getting sites', flush=True)
    for sample in samples:
        print(sample, flush=True)
        variant = FastaVariant(ref, vcf_file,
                               sample=sample, het=True, hom=True)
        print('got variant', flush=True)
        print(variant, flush=True)
        print(dir(variant[chrom]), flush=True)
        print(variant[chrom].variant_sites, flush=True)
        pos = list(variant[chrom].variant_sites)
        print('got %s pos' % len(pos), flush=True)
        sites.extend(pos)
        print('got %s sites' % len(sites), flush=True)

    sites = sorted(set(sites))
    print('using %s sites' % len(sites), flush=True)
    # get reference sequence for site positions
    refseq = []
    print('getting reference sequence', flush=True)
    for p in sites:
        refseq.append(reference[chrom][p - 1].seq)
    refseq = ''.join(refseq)
    # seqrecord for reference
    print('making reference sequence', flush=True)
    refrec = SeqRecord(Seq(refseq), id='ref')
    result.append(refrec)

    sites_matrix = {}
    print('getting sample sequences', flush=True)
    # iterate over variants in each sample
    for sample in samples:
        seq = []
        variant = FastaVariant(ref, vcf_file,
                               sample=sample, het=True, hom=True)
        for p in sites:
            rec = variant[chrom][p - 1:p]
            seq.append(rec.seq)
        seq = ''.join(seq)
        # make seqrecord for the samples sites
        seqrec = SeqRecord(Seq(seq), id=sample)
        result.append(seqrec)
        sites_matrix[sample] = list(seqrec)
    print('done', flush=True)
    df = pd.DataFrame(sites_matrix)
    df.index = sites
    output_file = vcf_file.replace('.vcf.gz', 'testout.fasta')
    SeqIO.write(result, output_file, "fasta")
    return result, df

def unitTest():
    test_vcf_file = "./testData/test.vcf.gz"
    test_ref_fasta = "./testData/test_ref.fa"

    # vcf_file = "/data/antwerpen/208/vsc20886/finaloutput.vcf.gz"
    # ref_fasta = "/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"

    test_output_fasta = './testOutput/output.fasta'

    result, df = fasta_alignment_from_vcf(test_vcf_file, test_ref_fasta)
    # vcf_to_fasta(test_vcf_file, test_ref_fasta, test_output_fasta)
    pprint(result)
if __name__ == "__main__":
    unitTest()
    # vcf_file = "/data/antwerpen/208/vsc20886/finaloutput.vcf.gz"
    # ref = "/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"
    # print("Running fasta_alignment_from_vcf", flush=True)
    # result, df = fasta_alignment_from_vcf(vcf_file, ref)
    # print(result)
    # print(df)