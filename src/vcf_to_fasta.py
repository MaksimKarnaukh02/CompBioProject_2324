import copy

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

    seq_records = {} # Create a list to store the resulting SeqRecord objects
    mutated_seqs = {}
    print('Processing VCF file')

    for record in vcf_reader: # Process each record in the VCF file
        # print(f"record: {record}")
        print(dir(record))
        chrom = record.CHROM
        # print(f"chrom: {chrom}")
        pos = record.POS
        # print(f"pos: {pos}")

        ref_allele = record.REF
        # print(f"ref_allele: {ref_allele}")
        ref_sequence = str(ref_records[chrom].seq) # Get the reference sequence for the current chromosome
        # print(f"ref_sequence: {ref_sequence}")
        if chrom not in seq_records.keys():
            seq_records[chrom] = []
            seq_records[chrom] = [Seq(ref_sequence).tomutable() for i in range(len(vcf_reader.samples))]

        for sample in record.samples:
            sample_name = sample.sample
            sample_index = record._sample_indexes[sample_name]
            # pprint(inspect.getmembers(sample.site))
            # print(f"sample_attributes: {sample}", flush=True)
            # print(f"sample alleles: {sample['GT']}", flush=True)
            # alleles = sample.site.allele_indexes # e.g. 1/0 => [1, 0] for a heterozygous sample
            alleles = sample['GT']
            allele_list = alleles.split("|")


            for allele_index in [ int(allele_index) for allele_index in allele_list]:
                """
                If you have different length of ALT and REF, 
                    you need to insert gaps ('-') in your alignment. 
                In your case it depends on what is in ALT. 
                    If REF is 'ATG' and ALT is 'A' you should insert 'A--' 
                    for 1|1, 'ATG' for 0|0 and 'NNN' or 'N--' for .|. , 
                    the last is arbitrary somehow and will be threated 
                    just the same downstream in analysis, 
                    but 'NNN' is more clear I think. 
                
                Bad things happen if REF is 'A' and ALT is 'ATG', then you need to insert 'A--' in all other sequences, 
                or you end up with unaligned sequences and the tree will be totally wrong. And that's why I recommended 
                to take readymade code for this.
                """

                if allele_index == 0: # no variation (no change compared to reference genome), so no change back needed
                    pass
                elif allele_index == "." or allele_index == "None" or allele_index is None or allele_index == "":
                    pass
                else: # variation (change compared to reference genome), so change back to original genome allele
                    if record.ALT[int(allele_index) - 1] == ".":
                        seq_records[chrom][sample_index][pos - 1] = ""
                    else:
                        alt_allele = record.ALT[int(allele_index) - 1]
                        seq_records[chrom][sample_index][pos - 1] = str(alt_allele)


                break  # break because we only look at the first number, so before the | or /
    final_seq_records = []
    for chrom, mutated_seqs in seq_records.items():
        for idx, mutated_seg in enumerate(mutated_seqs):
            id = f"{chrom}_{vcf_reader.samples[idx]}"
            final_seq_records.append(SeqRecord(mutated_seg, id=id))

    print('Finished processing VCF file')
    SeqIO.write(final_seq_records, output_fasta, 'fasta') # Save the SeqRecord objects to a FASTA file
    print('Finished writing FASTA file')


def unitTest():
    test_vcf_file = "./testData/test.vcf.gz"
    test_ref_fasta = "./testData/test_ref.fa"

    # vcf_file = "/data/antwerpen/208/vsc20886/finaloutput.vcf.gz"
    # ref_fasta = "/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"

    test_output_fasta = './testOutput/output.fasta'

    vcf_to_fasta(test_vcf_file, test_ref_fasta, test_output_fasta)

if __name__ == '__main__':
    unitTest()
    # vcf_file = "/data/antwerpen/208/vsc20886/finaloutput.vcf.gz"
    # ref_fasta = "/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"
    #
    # output_fasta = '/data/antwerpen/208/vsc20886/output.fasta'
    #
    # vcf_to_fasta(vcf_file, ref_fasta, output_fasta)