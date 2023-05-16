import pandas as pd
from Bio import Phylo

"""
check all nodes in a tree:
Get a list of OTUs under a node. 
Check if they belong to the same taxon. 
If yes - write down node number and count of OTUS(leaves) under this node (= LC).
Sort this list by LC. Iterate through list in a high to low LC order. For each row check if nodes lower in the list are 
descendants of this node, if they are - remove them from list.
Once you finished - you have a list of not nested nodes, each of them is a last common ancestor for a set of samples  
which belong to the same taxon. 
Then for each of these nodes you get a list of OTUs, pick two at random and it's done.
Just be careful to not to miss boundary cases 
(e.g. node which has only two children and they're samples of different species).
"""


def load_tree():
    tree = Phylo.read(r"data/inversion_samples_noninverted_pwd_nj_tree_ancestral_rooted.newick.nwk", "newick")
    return tree


def load_metadata():
    filename = "./Data/inversion_man_mt_2022-11-16.tsv"
    # load the metadata from inversion_man_mt_2022-11-16.tsv into a pandas dataframe
    df = pd.read_csv(filename, sep='\t')
    return df


def get_species_name(fish_id, dataframe):
    """
    Receives an id of a fish and looks up the species name in the metadata dataframe.
    """
    species_name = dataframe.loc[dataframe['id'] == fish_id, 'species'].iloc[0]
    return species_name


def check_taxon(OTUS, dataframe):
    """
    Check if a list of OTUs belong to the same taxon.
    """
    fish_name = OTUS[0].name
    # Get the species name of the first fish in the list
    species1 = get_species_name(fish_name, dataframe)
    # Loop over all the OTUs
    for OTU in OTUS[1:]:
        # Get the species name of the fish
        species2 = get_species_name(OTU.name, dataframe)
        # If the species name of the fish is different from the first fish in the list
        if species1 != species2:
            # Break out of the loop
            return False
    return True


def remove_doubles_LC(sorted_LC_list):
    """
    Removes double entries from the list of LCs.
    """
    i=0
    while i < len(sorted_LC_list):  # Iterate through list in a high to low LC order.
        children = sorted_LC_list[i][0].get_nonterminals()  # Get a list of children of the node
        j = i + 1
        species1 = get_species_name(sorted_LC_list[i][0].get_terminals()[0].name, dataframe)
        while j < len(sorted_LC_list):
            # For each entry check if nodes lower in the list are descendants of this node
            node_to_check = sorted_LC_list[j][0]
            species2 = get_species_name(sorted_LC_list[j][0].get_terminals()[0].name, dataframe)

            if node_to_check in children:
                sorted_LC_list.pop(j)  # Remove them from list
            else:
                j += 1
        i += 1
    return sorted_LC_list

dataframe = load_metadata()  # Load the metadata

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
        OTUS = node.get_terminals()  # Get a list of OTUs under a node.
        if check_taxon(OTUS, dataframe):  # Check if they belong to the same taxon.
            # Write down node number and count of OTUS(leaves) under this node (= LC).
            lc_list.append((node, len(OTUS)))
    lc_list.sort(key=lambda x: x[1], reverse=True)  # Sort this list by LC, high to low.
    print(len(lc_list))
    remove_doubles_LC(lc_list)  # Remove double entries from the list of LCs.
    print(len(lc_list))
    remove_doubles_LC(lc_list)  # Remove double entries from the list of LCs.
    print(len(lc_list))

    # Now we have a list of not nested nodes, each of them is a last common ancestor for a set of samples
    # which belong to the same taxon.
    # Then for each of these nodes you get a list of OTUs, pick two at random and it's done.
    # We pick the first two OTUs of each node.
    individuals_dataframe = pd.DataFrame(columns=['species', 'species2', 'id_individual_1', 'id_individual_2'])
    for node in lc_list:
        OTUS = node[0].get_terminals()
        # Get the species name of the fish
        species_name = get_species_name(OTUS[0].name, dataframe)
        species_name2 = get_species_name(OTUS[1].name, dataframe)
        # Add the species name and the first two OTUs of each node to the dataframe
        row = {'species': species_name, 'species2': species_name2, 'id_individual_1': OTUS[0].name, 'id_individual_2': OTUS[1].name}
        individuals_dataframe = pd.concat([individuals_dataframe, pd.DataFrame(row, index=[0])], ignore_index=True)
    toTSV(individuals_dataframe)
    return


def toTSV(dataframe):
    """
    Saves the dataframe to a tsv file.
    """
    dataframe.to_csv('individuals.tsv', sep='\t', index=False)

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
