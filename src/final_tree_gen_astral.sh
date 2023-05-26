#!/bin/bash

tree_folder="./Data/trees/"
one_tree_file="./Data/one_tree_file.nwk"
output_astral_tree="./Data/output_astral_tree.nwk"

# Step 1: Collect all newick trees into a single file
cat ${tree_folder}*.nwk > "${one_tree_file}"

#touch input_trees.nwk
#for file in ${tree_folder}*.nwk; do
#    sed -e '/^$/d' -e '/^\s*$/d' "$file" >> input_trees.nwk
#    echo >> input_trees.nwk
#done

# Step 2: Run ASTRAL
java -jar "./libs/Astral/astral.5.7.1.jar" -i ${one_tree_file} -o ${output_astral_tree}