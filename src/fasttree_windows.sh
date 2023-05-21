#!/bin/bash
# The trees will be made using the fasttree algorithm.
# The trees will be made using the fasta files located in the data/windowed_fastas folder.
# These trees will be saved in the data/trees folder.

# CD to the data folder
cd "/data/antwerpen/208/vsc20886" || exit
# Remove the trees folder if it exists
if [  -d "./trees" ]; then
    rm -r "./trees"
fi
# Make the trees folder
mkdir "./trees"

ls -lah
# Loop over all the fasta files in the windowed_fastas folder
for fasta_file in ./windowed_fastas/*.fa;
do
    # Get the name of the fasta file
    fasta_file_name=$(basename "$fasta_file")
    # Get the name of the tree file
    tree_file_name="${fasta_file_name%.*}.nwk"
    # Make the tree using the fasttree algorithm
    FastTree -gtr -nt  "$fasta_file" > "./trees/$tree_file_name"
done

