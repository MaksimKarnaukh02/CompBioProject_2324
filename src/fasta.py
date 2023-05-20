import subprocess
import os
# from pyfaidx import FastaVariant

def clip_reference_genome(ref_genome_file: str, vcf_file:str, output_file:str):
    """
    Clip the reference genome to the same extent as the VCF file.

    Args:
        ref_genome_file (str): Path to the reference genome FASTA file (Input).
        vcf_file (str): Path to the VCF file (Input).
        output_file (str): Path to the output file for the clipped reference genome (Output).
    """
    # Run bcftools to generate the consensus sequence
    print("Generating consensus sequence...")
    # Run bcftools to generate the consensus sequence
    cmd = ["bcftools", "consensus", "-f", ref_genome_file, vcf_file, "-o", output_file]
    print(" ".join(cmd))
    subprocess.run(cmd)
    print("\tDone.")
    # consensus = FastaVariant(ref_genome_file, vcf_file, het=True, hom=True)
    # out = open("variants.fasta", "w")
    # for chrom in consensus.keys:
    #     for var in consensus[chrom].variants_sites:
    #         record = consensus[chrom][var - 1:var]
    #         print(record, file=out)
    # out.close()

def split_fasta_file(fasta_file, window_size, output_directory):
    """
    Split a large FASTA file into smaller windows of specified size.

    Args:
        fasta_file (str): Path to the input FASTA file containing aligned sequences for all samples (Input).
        window_size (int): Size of the window in base pairs.
        output_directory (str): Path to the output directory to store the windowed FASTA files (Output).
    """

    # Split the large FASTA file into smaller windows
    with open(fasta_file, 'r') as input_file:
        sequence = ''
        for line in input_file:
            if line.startswith('>'):
                if sequence:
                    write_fasta_window(sequence, window_size, output_directory)
                    sequence = ''
                write_header(line, output_directory)
            else:
                sequence += line.strip()
        if sequence:
            write_fasta_window(sequence, window_size, output_directory)


def write_header(header, output_directory):
    """
    Write FASTA header to each window file.

    Args:
        header (str): Header line of the FASTA sequence.
        output_directory (str): Path to the output directory (Output).
    """
    with open(output_directory + '/header.txt', 'w') as header_file:
        header_file.write(header)


def write_fasta_window(sequence, window_size, output_directory):
    """
    Write a window of FASTA sequence to a separate file.

    Args:
        sequence (str): FASTA sequence.
        window_size (int): Size of the window in base pairs.
        output_directory (str): Path to the output directory (Output).
    """
    # Write a window of FASTA sequence to a separate file
    for i in range(0, len(sequence), window_size):
        window = sequence[i:i + window_size]
        window_file = output_directory + '/window_{}.fa'.format(i)
        with open(window_file, 'w') as output_file:
            output_file.write('>Window {}\n'.format(i))
            output_file.write(window)


if __name__ == '__main__':
    # Define the input and output file paths
    # Reference genome FASTA file (Input)
    ref_genome_file = "/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"
    vcf_file = os.path.expandvars('${VSC_DATA}/filtered.vcf.gz')   # VCF file (Input)
    output_file = os.path.expandvars('${VSC_DATA}/clipped_reference.fa')  # Output file for the clipped reference genome (Output)

    # Clip the reference genome
    clip_reference_genome(ref_genome_file, vcf_file, output_file)

    # Define the input and output file paths
    window_size = 100000  # Size of the window in base pairs
    output_directory = os.path.expandvars('${VSC_DATA}/windowed_fastas')  # Output directory to store the windowed FASTA files (Output)

    # Split the FASTA file into smaller windows
    split_fasta_file(output_file, window_size, output_directory)
