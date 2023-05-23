import vcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pprint import pprint
import inspect

def vcf_to_fasta(vcf_file, ref_fasta, output_fasta): #TODO test if actually works
    # Read the reference genome FASTA file
    print('Reading reference genome FASTA file')
    ref_records = SeqIO.to_dict(SeqIO.parse(ref_fasta, 'fasta'))
    print('Finished reading reference genome FASTA file')
    vcf_reader = vcf.Reader(open(vcf_file, 'rb')) # Open the VCF file

    seq_records = [] # Create a list to store the resulting SeqRecord objects

    print('Processing VCF file')
    for record in vcf_reader: # Process each record in the VCF file
        # print(f"record: {record}")

        chrom = record.CHROM
        # print(f"chrom: {chrom}")
        pos = record.POS
        # print(f"pos: {pos}")

        ref_allele = record.REF
        # print(f"ref_allele: {ref_allele}")
        ref_sequence = str(ref_records[chrom].seq) # Get the reference sequence for the current chromosome
        # print(f"ref_sequence: {ref_sequence}")
        seq = Seq(ref_sequence)
        seq = seq.tomutable()

        seq[pos - 1] = ref_allele # set to reference allele by default

        for sample in record.samples:
            sample_name = sample.sample
            mutated_seq = seq.toseq()

            alleles = sample.allele_indices # e.g. 1/0 => [1, 0] for a heterozygous sample
            for allele_index in alleles:
                if allele_index == 0: # no variation (no change compared to reference genome), so no change back needed
                    pass
                elif allele_index == "." or allele_index == "None" or allele_index is None or allele_index == "":
                    mutated_seq[pos - 1] = ""
                else: # variation (change compared to reference genome), so change back to original genome allele
                    alt_allele = record.ALT[allele_index - 1]
                    mutated_seq[pos - 1] = str(alt_allele)

                seq_record = SeqRecord(mutated_seq, id=f'{chrom}_{pos}_{sample_name}', description='')
                seq_records.append(seq_record)
                break # break because we only look at the first number, so before the | or /

    print('Finished processing VCF file')
    SeqIO.write(seq_records, output_fasta, 'fasta') # Save the SeqRecord objects to a FASTA file
    print('Finished writing FASTA file')
if __name__ == '__main__':
    vcf_file = "/data/antwerpen/208/vsc20886/finaloutput.vcf.gz"
    ref_fasta = "/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"

    output_fasta = '/data/antwerpen/208/vsc20886/output.fasta'

    vcf_to_fasta(vcf_file, ref_fasta, output_fasta)