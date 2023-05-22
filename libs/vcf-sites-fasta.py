import sys,os
import pandas as pd
from Bio import SeqIO
from pyfaidx import Fasta
from pyfaidx import FastaVariant
import vcf

def fasta_alignment_from_vcf(vcf_file, ref):
    """
    Get a fasta alignment for all snp sites in a multi
    sample vcf file, including the reference sequence.
    """

    #index vcf
    cmd = 'tabix -p vcf -f {i}'.format(i=vcf_file)
    tmp = subprocess.check_output(cmd,shell=True)
    #get samples from vcf
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
    samples = vcf_reader.samples
    print ('%s samples' %len(samples))
    result = []

    #reference sequence
    reference = Fasta(ref)
    chrom = list(reference.keys())[0]

    #get the set of all sites first
    sites=[]
    for sample in samples:
        variant = FastaVariant(ref, vcf_file,
                                 sample=sample, het=True, hom=True)
        pos = list(variant[chrom].variant_sites)
        sites.extend(pos)

    sites = sorted(set(sites))
    print ('using %s sites' %len(sites))
    #get reference sequence for site positions
    refseq=[]
    for p in sites:
        refseq.append(reference[chrom][p-1].seq)
    refseq = ''.join(refseq)
    #seqrecord for reference
    refrec = SeqRecord(Seq(refseq),id='ref')
    result.append(refrec)

    sites_matrix = {}
    #iterate over variants in each sample
    for sample in samples:
        seq=[]
        variant = FastaVariant(ref, vcf_file,
                                 sample=sample, het=True, hom=True)
        for p in sites:
            rec = variant[chrom][p-1:p]
            seq.append(rec.seq)
        seq = ''.join(seq)
        #make seqrecord for the samples sites
        seqrec = SeqRecord(Seq(seq),id=sample)
        result.append(seqrec)
        sites_matrix[sample] = list(seqrec)
    df = pd.DataFrame(sites_matrix)
    df.index=sites
    return result, df
