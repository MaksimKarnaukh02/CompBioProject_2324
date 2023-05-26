import os
import subprocess
from typing import List, Tuple

import Bio.Phylo.Newick
import pandas as pd
from Bio import Phylo
Clade = Bio.Phylo.Newick.Clade

"""
check all nodes in a tree:
Get a list of OTUs under a node. 
Check if they belong to the same taxon. 
If yes - write down node number and count of OTUS(leaves) under this node (= LC).
Sort this list by LC. Iterate through list in a high to low LC order. For each row check if nodes lower in the list are 
descendants of this node, if they are - remove them from list.

info is in metadata csv.
"""


def load_tree():
    tree = Phylo.read(r"data/inversion_samples_noninverted_pwd_nj_tree_ancestral_rooted.newick.nwk", "newick")
    return tree


def load_tsv_in_dataframe(filename, seperator='\t'):
    return pd.read_csv(filename, sep=seperator)


def load_metadata():
    filename = "./Data/inversion_man_mt_2022-11-16.tsv"
    # load the metadata from inversion_man_mt_2022-11-16.tsv into a pandas dataframe
    return load_tsv_in_dataframe(filename)


def load_individuals():
    """
    Loads the individuals from the metadata.
    """
    filepath = 'individuals.tsv'
    return load_tsv_in_dataframe(filepath)


def get_species_name(fish_id, dataframe):
    """
    Receives an id of a fish and looks up the species name in the metadata dataframe.
    """
    species_name = dataframe.loc[dataframe['id'] == fish_id, 'species'].iloc[0]
    return species_name


def get_genus_name(fish_id, dataframe):
    """
    Receives an id of a fish and looks up the genus in the metadata dataframe.
    """
    genus_name = dataframe.loc[dataframe['id'] == fish_id, 'genus'].iloc[0]
    return genus_name


def get_taxus_identifier(fish_id, dataframe):
    """
    Receives an id of a fish and returns the taxus.
    """
    return get_genus_name(fish_id, dataframe) + " / " + get_species_name(fish_id, dataframe)


def check_taxon(OTUS, dataframe):
    """
    Check if a list of OTUs belong to the same taxon.
    """
    fish_name = OTUS[0].name
    # Get the species name of the first fish in the list
    taxus_identifier1 = get_taxus_identifier(fish_name, dataframe)
    # Loop over all the OTUs
    for OTU in OTUS[1:]:
        # Get the species name of the fish
        taxus_identifier2 = get_taxus_identifier(OTU.name, dataframe)
        # If the species name of the fish is different from the first fish in the list
        if taxus_identifier1 != taxus_identifier2:
            # Break out of the loop
            return False
    return True


def filter_doubles_LC(sorted_LC_list):
    """
    Removes double entries from the list of LCs.
    """
    i = 0
    while i < len(sorted_LC_list):  # Iterate through list in a high to low LC order.
        children = sorted_LC_list[i][0].get_nonterminals()  # Get a list of children of the node
        j = i + 1
        while j < len(sorted_LC_list):
            # For each entry check if nodes lower in the list are descendants of this node
            node_to_check = sorted_LC_list[j][0]

            if node_to_check in children:
                sorted_LC_list.pop(j)  # Remove them from list
            else:
                j += 1
        i += 1
    return sorted_LC_list


def isAnsestor(node, dataframe):
    OTUS = node.get_terminals()  # Get a list of OTUs under a node.
    return check_taxon(OTUS, dataframe)  # Check if they belong to the same taxon.



def get_existing_individuals():
    """
    Returns a list of individuals that exist in the VCF.
    """
    command = """zgrep -m 1 '^#CHROM' "../malawi_cichlids_v3_phase.biallelic_snps.pass.ancestral_as_sample.benthic.chr1.vcf.gz" | awk -F "FORMAT\t" '{print $2}' | awk -F "\t
ancestral" '{print $1}'"""
    existing_individuals_string = subprocess.check_output(command, shell=True, encoding="utf-8")
    existing_individuals_list = set(existing_individuals_string.split("\t"))
    return existing_individuals_list


def filter_individuals_by_existence_in_VCF(individuals: List[Tuple[Clade, int]]):
    print(len(individuals))
    existing_individuals = get_existing_individuals()
    print(len(existing_individuals))
    indx = 0
    while indx < len(individuals):
        clade: Clade = individuals[indx][0]
        if clade.name not in existing_individuals:
            individuals.pop(indx)
        else:
            indx += 1
    print(len(individuals))
    return individuals

def filter_LC(LC_list):
    """
    Filters the list of LCs.
    """
    LC_list.sort(key=lambda x: x[1], reverse=True)  # Sort this list by LC, high to low.

    # Remove double entries from the list of LCs.
    LC_list = filter_doubles_LC(LC_list)

    return LC_list

def findTwoIndividuals(OTUS, existing_individuals):
    """
    Finds two individuals from the same taxon.
    """
    individuals = []
    for OTU in OTUS:
        if OTU.name in existing_individuals:
            individuals.append(OTU.name)
            if len(individuals) == 2:
                return individuals[0], individuals[1]
    return None, None
def select_individuals():
    """
    Selects individuals from the tree.
    """
    tree = load_tree()  # Load the tree
    dataframe = load_metadata()  # Load the metadata

    lc_list = []  # pairs of (node number, count of OTUS(leaves) under this node (= LC))
    species: str  # species name of the fish

    nonterminals = tree.get_nonterminals()
    for node in nonterminals:  # Iterate over all nodes in the tree
        if isAnsestor(node, dataframe):
            # Write down node number and count of OTUS(leaves) under this node (= LC).
            lc_list.append((node, len(node.get_terminals())))

    lc_list = filter_LC(lc_list)  # Filter function for the LC list.
    # Now we have a list of not nested nodes, each of them is a last common ancestor for a set of samples
    # which belong to the same taxon.
    # Then for each of these nodes you get a list of OTUs, pick two at random and it's done.
    # We pick the first two OTUs of each node.
    individuals_dataframe = pd.DataFrame(columns=['taxus', 'id_individual_1', 'id_individual_2'])
    existing_individuals_in_lake = get_existing_individuals()
    for node in lc_list:
        OTUS = node[0].get_terminals()
        # Get the species name of the fish
        taxus_id = get_taxus_identifier(OTUS[0].name, dataframe)
        id_individual_1, id_individual_2 = findTwoIndividuals(OTUS, existing_individuals_in_lake)
        # If 2 fishes exist in the lake
        if id_individual_1 and id_individual_2:
            # Add the species name and the first two OTUs of each node to the dataframe
            row = {'taxus': taxus_id, 'id_individual_1': id_individual_1, 'id_individual_2': id_individual_2}
            individuals_dataframe = pd.concat([individuals_dataframe, pd.DataFrame(row, index=[0])], ignore_index=True)
    toTSV(individuals_dataframe)
    return individuals_dataframe


def toTSV(dataframe):
    """
    Saves the dataframe to a tsv file.
    """
    dataframe.to_csv('individuals.tsv', sep='\t', index=False)

def get_individuals():
    """
    Returns a dataframe with the individuals.
    If the individuals.tsv file exists, it loads it. otherwise it creates it.
    """
    individuals_filepath = 'individuals.tsv'
    if os.path.isfile(individuals_filepath):
        return load_tsv_in_dataframe(individuals_filepath)
    else:
        return select_individuals()

def test():
    tree = Phylo.read(r"data/inversion_samples_noninverted_pwd_nj_tree_ancestral_rooted.newick.nwk", "newick")
    # now tree is a root node of the tree
    # let's get all internal nodes:
    nodes = tree.get_nonterminals()

    # we can take random node
    node = nodes[42]

    # and then we can get it's OTUs (leaves, terminals, tips here are synonyms)
    # as nodes
    OTUS = node.get_terminals()
    # or just labels/names of this terminals
    print([n.name for n in OTUS])

    # you can also get children of a node
    c1, c2 = node.clades
    print([n.name for n in c1.get_terminals()])
    print([n.name for n in c2.get_terminals()])

    # you can compare nodes to find if they're the same node
    print([n for n in nodes if n == c1])
    print([n for n in nodes if n == c2])


if __name__ == "__main__":
    individuals: pd.DataFrame = select_individuals()
    # individuals = get_existing_individuals()
    # print(individuals)

