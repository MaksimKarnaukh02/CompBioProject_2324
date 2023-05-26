import os
import numpy as np
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

from Bio import Phylo


def collect_tree_files(tree_directory):
    ls = os.listdir(tree_directory)

    treeFiles = [f"{tree_directory}{file}" for file in ls if file.endswith('.nwk')]
    return treeFiles


def compute_distance(tree_file1, tree_file2):

    tree1 = Phylo.read(tree_file1, 'newick')
    tree2 = Phylo.read(tree_file2, 'newick')

    rf_distance = tree1.compare(tree2, 'rf')
    return rf_distance


def compute_topo_distances(tree_files):
    distances = np.zeros((len(tree_files), len(tree_files)))
    for i, tree_file1 in enumerate(tree_files):
        for j, tree_file2 in enumerate(tree_files):
            # Compute the topological distance between tree_file1 and tree_file2
            distance = compute_distance(tree_file1, tree_file2)  # Implement your own distance calculation method
            distances[i, j] = distance
    return distances


def dim_reduction(distances):
    mds = MDS(n_components=2, dissimilarity='precomputed')
    reduced_distances = mds.fit_transform(distances)
    return reduced_distances


def plot_treespace(reduced_distances):
    plt.scatter(reduced_distances[:, 0], reduced_distances[:, 1])
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.title('Treespace Plot')
    plt.show()



if __name__ == '__main__':

    # Step 1: Collect Tree Files
    _tree_directory = f'./Data/trees/'
    _tree_files = collect_tree_files(_tree_directory)

    # Step 2 & 3: Compute Topological Distances and store them
    _distances = compute_topo_distances(_tree_files)

    # Step 4: Perform Dimensionality Reduction (MDS)
    _reduced_distances = dim_reduction(_distances)

    # Step 5: Plot Treespace
    plot_treespace(_reduced_distances)



