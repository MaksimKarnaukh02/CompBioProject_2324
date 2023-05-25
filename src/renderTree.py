import os,sys
from ete3 import Tree
os.environ['QT_QPA_PLATFORM']='offscreen'

# if first argument does not exist
if len(sys.argv) < 2:
    exit(1)

# grab first argument
window_n = int(sys.argv[1])


# if directory Data/rendered_trees does not exist, create it
if not os.path.exists("Data/rendered_trees"):
    os.makedirs("Data/rendered_trees")


window_size=100000
window_start=window_size*window_n



newick_file=f"./Data/trees/window.chr1.{window_start}.aln.min4.rehydrated.nwk"
output_png = f"Data/rendered_trees/window.chr1.{window_start}.aln.min4.rehydrated.png"
# Hide leaf labels
tree = Tree(newick_file)
tree.render(output_png)
