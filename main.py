"""
Maybe it's better to subset windows, dropping those, which have less positions marked as "good"
than some threshold (say 90%). Otherwise you'll get some windows with only a few sites
available to analysis and it can affect topologies badly. I have a function in python
which computes the count of accessible positions in python (get_accessible_size in the attached file),
but you can do it any other way as well.

you should split your chromosome in windows.
For this you can generate a list of non-overlapping windows i.e. pairs of coordinates (start,end),
then check with this function how many positions are good within this window.
It will be a number in range from 0 to end-start.
Then you can get percentage of good sites by dividing it by window length.
The function calls tabix from SAMtools, so load module BioTools , it contains SAMtools and other useful things as well.
Chromosome parameter is just the chromosome name (look at the VCF file name its "chr[something]")

The idea was to split chromosome in 100k windows, i.e. the first window: 1 to 100000 ,
the second 100001 to 200000, the third 200001 to 300000,
but the window size can be changed if it turns out that there is not enough information
in 100kbp to build trees or if it will be too computationally
costly to build 41M(length of the chromosome)/100K(window size) trees.

VCF apart from it's header (lines starting with "#" ) is just a table -
rows are variable sites (a bit more complex in fact, but good enough so far), columns are individuals.
First 9 columns are not individuals, but info about sites - chromosome name,  coordinate, which nucleotide
changed to which etc. So it's basically just "Filter rows which have in field 'coordinate'
a value which is between S and E, in at least one pair of (S,E) in a list. list of pairs (S,E)
is a list of coordinates of good windows Starts and Ends".

"""

import sys
import os
import subprocess
import pandas as pd

import src.python_tree_traversing as ptt


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
    chrom_len: int = get_chr_length(reference_genome, chrom) # 41162407 in file 'GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.dict' for chr1
    good_windows = []
    for start in range(0, chrom_len, window_size):
        end = min(start + window_size, chrom_len)
        if get_accessible_size(accessible_fn, chrom, start, end) / (end - start) > 0.9: # # of accessible sites / window size
            good_windows.append((start, end))
    return good_windows


def main():
    individuals: pd.DataFrame = ptt.get_individuals()

    # get_good_windows("malawi_cichlids_callset_v3_qc_subset_chr1_pass.bed.gz", "chr1")

    return 0


if __name__ == "__main__":
    main()

