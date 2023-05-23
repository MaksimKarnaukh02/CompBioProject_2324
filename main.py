"""
subset windows, dropping those, which have less positions marked as "good"
than some threshold (say 90%). Otherwise we get some windows with only a few sites
available to analysis and it can affect topologies badly. (get_accessible_size in the attached file)

split chromosome in windows.
For this we generate a list of non-overlapping windows i.e. pairs of coordinates (start,end),
then check with this function how many positions are good within this window.
It will be a number in range from 0 to end-start.
Then get percentage of good sites by dividing it by window length.
The function calls tabix from SAMtools, so load module BioTools, it contains SAMtools and other useful things as well.

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

import src.python_tree_traversing as ptt
import src.windows as w


def main():
    # individuals: pd.DataFrame = ptt.get_individuals()

    window_list = w.get_good_windows("./Data/malawi_cichlids_callset_v3_qc_subset_chr1_pass.bed.gz", "chr1")
    print(window_list)

    return 0


if __name__ == "__main__":
    main()

