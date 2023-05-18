import sys
import os
import subprocess
import pandas as pd


def get_accessible_size(accessible_fn, chrom, start=None, end=None):
    """
    Returns accessible size in bp of a given genomic region from a tabix
    indexed bed.gz file.

    Requires tabix to be installed and in the PATH.
    """
    region_str = str(chrom)
    if start is not None:
        assert end is not None, "start and end most both be given or none"
        region_str += f":{start}-{end}"
    elif end is not None:
        raise Exception("start and end most both be given or none")

    p = subprocess.Popen(f"tabix {accessible_fn} {region_str}", shell=True, stdout=subprocess.PIPE)
    d = pd.read_csv(p.stdout, sep='\t', header=None,
                    names=['chrom', 'start', 'end'], usecols=['start', 'end'])
    # check that intervals are non overlapping
    assert (d.iloc[:-1]['end'].values <= d.iloc[1:]['start'].values).all()
    try:
        d.iloc[0]['start'] = start
        d.iloc[-1]['end'] = end
        access = (d['end'] - d['start']).sum()
    except IndexError:
        access = 0
    return access


def get_chr_length(reference_genome, chrom):
    """
    Returns the length of the given chromosome using a reference genome.
    """
    return 41162407
    with open(reference_genome, 'r') as file:
        for line in file:
            if line.startswith(f"@SQ\tSN:{chrom}"):
                length_field = line.split('\t')[1]
                chrom_len = int(length_field.split(':')[1])
                return chrom_len
    raise ValueError(f"Chromosome '{chrom}' not found in the reference genome file.")


def get_good_windows(accessible_fn, chrom, window_size=100000, reference_genome=None):
    """
    Returns a list of (start, end) tuples of windows of size window_size
    that have at least 90% of sites accessible in the given chromosome.
    """
    chrom_len: int = get_chr_length(reference_genome,
                                    chrom)  # 41162407 in file 'GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.dict' for chr1
    print(chrom_len)
    good_windows = []
    for start in range(0, chrom_len, window_size):
        end = min(start + window_size, chrom_len)
        if get_accessible_size(accessible_fn, chrom, start, end) / (
                end - start) > 0.9:  # # of accessible sites / window size
            good_windows.append((start, end))
    return good_windows


def toBED(data: list, output_file):
    """
    Converts list of pairs of positions to BED file
    """
    with open(output_file, 'w') as f:
        for entry in data:
            pos1, pos2 = entry
            chrom = "chr1"
            bed_line = f"{chrom}\t{pos1}\t{pos2}\n"
            f.write(bed_line)


if __name__ == "__main__":
    window_list = get_good_windows("./Data/malawi_cichlids_callset_v3_qc_subset_chr1_pass.bed.gz", "chr1")
    toBED(window_list, "./chr1_windows.bed")
