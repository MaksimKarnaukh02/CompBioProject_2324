#!/bin/bash
# The trees will be made using the fasttree algorithm.
# The trees will be made using the fasta files located in the data/windowed_fastas folder.
# These trees will be saved in the data/trees folder.

# CD to the data folder
cd "/scratch/antwerpen/208/vsc20886/CompBioProject_2324" || exit

if [ -z "$1" ]
then
  echo "Need an argument <window_number>"
  exit 1;
fi
window_n=${1}
window_size=100000
window_start=$((window_size*window_n))
window_end=$((window_start + window_size))
window_start_1_based=$((window_start + 1))

echo "Making trees using the fasttree algorithm"

fasta_file="./Data/windowed_fastas/window.chr1.${window_start}.aln.min4.rehydrated.fasta"
echo "    Making tree for $fasta_file"
# Get the name of the fasta file
fasta_file_name=$(basename "$fasta_file")
# Get the name of the tree file
tree_file_name="${fasta_file_name%.*}.nwk"
# Make the tree using the fasttree algorithm
/scratch/antwerpen/208/vsc20886/CompBioProject_2324/libs/FastTree -gtr -nt  "$fasta_file" > "./Data/trees/$tree_file_name"
echo "    Done"

