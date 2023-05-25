import sys
from os import path

from Bio import AlignIO
from Bio import SeqIO

assert len(
    sys.argv) == 5, "Not enough args given. Should be: <converted_fasta> <coord_tsv_file> <ref_genome_fn> <window_start>"

# fasta with variable sites only, produced by 
# https://github.com/edgardomortiz/vcf2phylip
converted_fasta = sys.argv[1]

out_base = path.splitext(converted_fasta)[0]
out_file = out_base + ".rehydrated.fasta"

# tsv of var sites' coordinates
# produced by vcf2phylip
# vcf2phylip/rec.chr1.10001-20000.fasta.min4.used_sites.tsv 
coord_tsv_file = sys.argv[2]

# ref genome
ref_genome_fn = sys.argv[3]

assert sys.argv[4].isdigit(), f"Window size should be an integer number, but given: {sys.argv[4]}"

# zero-based!
window_start = int(sys.argv[4])
window_size = 100000
print_every = 10000
reference = SeqIO.read(ref_genome_fn, "fasta")

aln = AlignIO.read(converted_fasta, "fasta")
samples = [s.name for s in aln]
sample_fas = ["" for _ in range(0, len(samples))]

old_chr = None
old_pos = window_start
last_printed = window_start
subst_site = 0
coordinate_csv = open(coord_tsv_file)
print('Rehydrating fasta:', end='', flush=True)
for n, line in enumerate(coordinate_csv):
    if line.startswith("#CHROM"):
        continue
    # if n>10:
    #  break
    chr, pos, _ = line.rstrip().split("\t")
    assert pos.isdigit(), f"Wrong value in position field of coordinate csv: {pos}"
    pos = int(pos)
    if pos - last_printed > print_every:
        print('.' * int((pos - last_printed) / print_every), end='', flush=True)
        last_printed = pos
    if old_chr is not None:
        assert old_chr == chr, "Only one chromosome allowed"
    else:
        old_chr = chr
    if old_pos is not None and pos - old_pos > 1:
        insert = str(reference[old_pos:pos - 1].seq)
        for sn in range(0, len(samples)):
            sample_fas[sn] += insert
    else:
        insert = ""
    for sn, nuc in zip(range(0, len(samples)), str(aln[:, subst_site])):
        # if sn==0:
        #  print (f"{insert} {nuc}\t{n} {subst_site} {old_pos} {pos}")
        sample_fas[sn] += nuc
    old_pos = pos
    subst_site += 1

insert = str(reference[old_pos:window_start + window_size].seq)
for sn in range(0, len(samples)):
    sample_fas[sn] += insert

with open(out_file, "w") as out:
    for sample_name, sequence in zip(samples, sample_fas):
        _ = out.write(f">{sample_name}\n")
        _ = out.write(sequence + "\n")

print('Done!', flush=True)
