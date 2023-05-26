import shutil
import subprocess
import os

# from pyfaidx import FastaVariant

def clip_reference_genome(ref_genome_file: str, vcf_file: str, output_file: str):
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



def split_fasta_file(fasta_file, window_size, output_directory):
    """
    Split a large FASTA file into smaller windows of specified size.

    Args:
        fasta_file (str): Path to the input FASTA file containing aligned sequences for all samples (Input).
        window_size (int): Size of the window in base pairs.
        output_directory (str): Path to the output directory to store the windowed FASTA files (Output).
    """
    print(f"Splitting FASTA file into windows with windows size {window_size}...")
    # Split the large FASTA file into smaller windows

    sequences = {}  # Dictionary to store sequences per header
    with open(fasta_file, 'r') as input_file:

        # Process each line in the file
        for line in input_file:
            line = line.strip()

            if line.startswith('>'):
                # If it is a header line, store the header
                current_header = line
                sequences[current_header] = []
            else:
                # If it is a sequence line, append the sequence to the corresponding header
                sequences[current_header].append(line)

    # Create the output directory if it doesn't exist
    # os.makedirs(output_directory, exist_ok=True)

    length_first_sequence = 0

    for letter_seq in sequences[list(sequences.keys())[0]]:
        length_first_sequence += len(letter_seq)

    # Split the sequences into windows and write to separate FASTA files
    for i in range(0, length_first_sequence, window_size):
        windowed_sequences = {header: ''.join(sequence)[i:i + window_size] for header, sequence in sequences.items()}

        output_filename = os.path.join(output_directory, f'window_{i}-{i + window_size}.fasta')
        with open(output_filename, 'w') as output_file:
            for header, sequence in windowed_sequences.items():
                output_file.write(f'{header}\n{sequence}\n')

    print("\tDone.")


VSC_DATA = '/data/antwerpen/208/vsc20886'
VSC_SCRATCH = '/scratch/antwerpen/208/vsc20886'

if __name__ == '__main__':
    # Define the input and output file paths
    # Reference genome FASTA file (Input)
    ref_genome_file = f"/scratch/antwerpen/208/vsc20811/2024-05_compbio_project/genome/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"
    vcf_file = f'{VSC_DATA}/finaloutput.vcf.gz'  # VCF file (Input)
    output_file = f'{VSC_DATA}/clipped_reference.fa'  # Output file for the clipped reference genome (Output)

    # Clip the reference genome
    clip_reference_genome(ref_genome_file, vcf_file, output_file)

    # Define the input and output file paths
    window_size = 100000  # Size of the window in base pairs
    output_directory = f'{VSC_DATA}/windowed_fastas'  # Output directory to store the windowed FASTA files (Output)
    # Create the output directory if it does not exist
    if os.path.exists(output_directory):
        # Remove the directory if it already exists
        shutil.rmtree(output_directory)
    os.makedirs(output_directory)

    # Split the FASTA file into smaller windows
    split_fasta_file(output_file, window_size, output_directory)
