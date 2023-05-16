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
The function calls tabix from SAMtools, so  load module BioTools , it contains SAMtools and other useful things as well.
Chromosome parameter is just the chromosome name (look at the VCF file name its "chr[something]")
"""

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
