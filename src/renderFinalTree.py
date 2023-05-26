import random
from pprint import pprint

import pandas as pd
from ete3 import Tree, AttrFace, TreeStyle


def load_tree(tree_file):
    tree = Tree(tree_file)
    return tree


def load_metadata(metadata_file):
    metadata = pd.read_csv(metadata_file, sep='\t')
    return metadata


def customTreeStyle():
    ts: TreeStyle = TreeStyle()
    ts.show_leaf_name = False
    ts.tree_width = 1500
    ts.complete_branch_lines_when_necessary = True
    ts.draw_guiding_lines = True
    ts.guiding_lines_color = "black"
    ts.guiding_lines_type = 2
    return ts


def render_tree(tree, output_png):
    ts = customTreeStyle()
    tree.render(output_png, tree_style=ts)


def rename_labels(tree, metadata):
    """
    Receives a tree and a metadata dataframe and renames the labels of the tree with the genus and species names.
    """
    for leaf in tree.get_leaves():
        leaf.id = leaf.name
        # if the leaf name is in the metadata dataframe
        if leaf.name in metadata['id'].values:
            # Get the species name of the fish
            species_name = metadata.loc[metadata['id'] == leaf.name, 'species'].iloc[0]
            # Get the genus name of the fish
            genus_name = metadata.loc[metadata['id'] == leaf.name, 'genus'].iloc[0]
            # Rename the leaf
            leaf.name = genus_name + " " + species_name
        else:
            leaf.name = "unknown"

    return tree


def random_color():
    """
    Returns a random color in hex format, only light colors so that black text is readable on it.
    """
    color = '#%02X%02X%02X' % (random.randint(100, 255), random.randint(100, 255), random.randint(100, 255))
    return color


def getGenus(node, metadata):
    """
    Receives a node and a metadata dataframe and returns the species name of the node.
    Returns None if the node has different genusses in it's leaves
    """
    # take the first leaf of the node
    leaf = node.get_leaves()[0]
    # Get the genus name of the fish
    genus_name_1 = metadata.loc[metadata['id'] == leaf.id, 'genus'].iloc[0]
    # loop over all the leaves of the node
    for leaf in node.get_leaves():
        # Get the genus name of the fish
        genus_name_2 = metadata.loc[metadata['id'] == leaf.id, 'genus'].iloc[0]
        # If the genus name of the fish is different from the first fish in the list
        if genus_name_1 != genus_name_2:
            # Break out of the loop
            return None
    return genus_name_1


def recolorNodesBySupport(tree):
    for node in tree.traverse():
        node.img_style["fgcolor"] = red_black_gradient(node.support)
        node.img_style["size"] = 10
    return tree


def recolor_leaves_by_genus(tree, metadata):
    """
    Receives a tree and a metadata dataframe and colors the branches by genus.
    """
    colormap = getColorMap()

    for node in tree.get_leaves():
        genus_name = getGenus(node, metadata)

        attrFace: AttrFace = AttrFace("name")

        if genus_name in colormap.keys():
            attrFace.background.color = colormap[genus_name]
        node.add_face(attrFace, 1, position='aligned')
    pprint(colormap)
    return tree


def red_black_gradient(percentage):
    """
    Receives a percentage (0-1) and returns a color in hex format.
    """
    red = int((1 - percentage) * 255)
    green = 0
    blue = 0
    color = '#%02X%02X%02X' % (red, green, blue)
    return color


def getColorMap():
    return {
        'Copadichromis': '#ff5500',
        'Otopharynx': '#ff8b06',
        'Protomelas': '#ffba0b',
        'Lethrinops': '#ffeb11',
        'Aulonocara': '#e4ff17',
        'Alticorpus': '#bbff1c',
        'Taeniolethrinops': '#91ff22',
        'Tramitichromis': '#68ff27',
        'Mchenga': '#46ff2d',
        'Mylochromis': '#33ff44',
        'Placidochromis': '#38ff70',
        'Ctenopharynx': '#3eff98',
        'Fossorochromis': '#44ffc1',
        'Sciaenochromis': '#49ffe7',
        'Stigmatochromis': '#4ff6ff',
        'Buccochromis': '#54d4ff',
        'Champsochromis': '#5ab5ff',
        'Tyrannochromis': '#609aff',
        'Nimbochromis': '#657fff',
        'Dimidiochromis': '#706bff',
        'Hemitaeniochromis': '#9071ff',
        'Hemitilapia': '#af76ff',
        'Cheilochromis': '#cd7cff',
        'Chilotilapia': '#e682ff',
        'Trematocranus': '#ff87ff',

    }


def remove_leave_names(tree):
    for node in tree.get_leaves():
        node.name = ""
    return tree


def setVirginsTreeRoot(tree):
    """
        Sets the root of the tree to the common ancestors node for all the nodes
        with name "Copadichromis virginalis" or "Copadichromis sp-virginalis-kajose"
    """
    allowed_genus = ["Copadichromis virginalis", "Copadichromis sp-virginalis-kajose"]
    # loop over all the leaves
    found_virgins = []
    for leaf in tree.get_leaves():
        if leaf.name in allowed_genus:
            found_virgins.append(leaf)
    common_ancestor = found_virgins[0].get_common_ancestor(found_virgins[1:])
    tree.set_outgroup(common_ancestor)


def FinalTreeRender():
    newick_file = f"./Data/output_astral_tree.nwk"
    output_png = f"./Data/output_astral_tree.svg"
    metadata_file = "./Data/inversion_man_mt_2022-11-16.tsv"
    tree: Tree = load_tree(newick_file)
    metadata: pd.DataFrame = load_metadata(metadata_file)

    # Use metadata file to rename branches and color them by genus.
    tree = rename_labels(tree, metadata)

    setVirginsTreeRoot(tree)

    # Recolor branches by genus
    tree = recolor_leaves_by_genus(tree, metadata)
    tree = recolorNodesBySupport(tree)
    # tree = remove_leave_names(tree)

    render_tree(tree, output_png)


def InitialTreeRender():
    """
    Renders the initial tree with the leaf names colored by genus.
    """
    newick_file = f"./Data/inversion_samples_noninverted_pwd_nj_tree_ancestral_rooted.newick.nwk"
    output_png = f"./Data/initial_tree.svg"
    metadata_file = "./Data/inversion_man_mt_2022-11-16.tsv"

    tree: Tree = load_tree(newick_file)
    metadata: pd.DataFrame = load_metadata(metadata_file)
    tree = rename_labels(tree, metadata)

    tree = recolor_leaves_by_genus(tree, metadata)
    tree = recolorNodesBySupport(tree)
    render_tree(tree, output_png)


if __name__ == "__main__":
    FinalTreeRender()
    InitialTreeRender()

# Also you need to indicate supports (values indicating which percentage of quartets are concordant with your final topology).
