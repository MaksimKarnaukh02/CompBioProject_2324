from Bio import Phylo

"""
OK, easiest (not optimal, but performance here doesn't matter too much) way to do that is to check all nodes in a tree:
Get a list of OTUs under a node. Check if they belong to the same taxon. If yes - write down node number and count of OTUS(leaves) under this node (LC)
Sort this list by LC. Iterate through list in a high to low LC order. For each row check if nodes lower in the list are descendants of this node, if they are - remove them from list.
Once you finished - you have a list of not nested nodes, each of them is an last common ancestor for a set of samples  which belong to the same taxon. Then for each of these nodes you get a list of OTUs, pick two at random and it's done.
Just be careful to not to miss boundary cases (e.g. node which has only two children and they're samples of different species).
"""


def select_individuals():
    tree = Phylo.read(r"inversion_samples_noninverted_pwd_nj_tree_ancestral_rooted.newick.nwk", "newick")
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
