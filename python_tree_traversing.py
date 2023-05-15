from Bio import Phylo
tree=Phylo.read(r"inversion_samples_noninverted_pwd_nj_tree_ancestral_rooted.newick.nwk")
# now tree is a root node of the tree
# let's get all internal nodes:
nodes=tree.get_nonterminals()

# we can take random node
node=nodes[42]

# and then we can get it's OTUs (leaves, terminals, tips here are synonyms)
# as nodes
OTUS=node.get_terminals()
# or just labels/names of this terminals
print([n.name for n in OTUS])

# you can also get children of a node
c1,c2=node.clades
print([n.name for n in c1.get_terminals()])
print([n.name for n in c2.get_terminals()])

# you can compare nodes to find if they're the same node
print([n for n in nodes if n==c1])
print([n for n in nodes if n==c2])