DESCRIPTION:
-----------
ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees.
ASTRAL is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS).
ASTRAL finds the species tree that has the maximum number of shared induced quartet trees with the set of gene trees, subject to the constraint that the set of bipartitions in the species tree comes from a predefined set of bipartitions. This predefined set is empirically decided by ASTRAL (but see tutorial on how to expand it). The current code corresponds to **ASTRAL-III** (see below for the publication).


The algorithm was designed by Tandy Warnow and Siavash Mirarab originally. ASTRAL-III incorporates many ideas by Chao Zhang and Maryam Rabiee.
[Code developers](https://github.com/smirarab/ASTRAL/graphs/contributors) are mainly Siavash Mirarab, Chao Zhang, Maryam Rabiee, and Erfan Sayyari.

Email: `astral-users@googlegroups.com` for questions.


### Other branches

**NOTE**: 
Several new  features of ASTRAL are not merged in this branch and are available in other branches. 
Please use those branches if you find these features useful. 

* **Multi-threaded ASTRAL**: A multi-threaded version of ASTRAL is now available on [this branch](https://github.com/smirarab/ASTRAL/tree/MP)
* **Astral with user constraints**: A version of ASTRAL that can satisfy user constraints is available [here](https://github.com/maryamrabiee/Constrained-search)
* **Tree updates**:  An ASTRAL-based algorithm called INSTRAL enables inserting  new species onto and existing ASTRAL tree is available [here](https://github.com/maryamrabiee/INSTRAL)

## Publications:

#### Papers on the current version:
- Since version 5.1.1, the code corresponds to **ASTRAL-III**, described in:
    * Zhang, Chao, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153. [doi:10.1186/s12859-018-2129-y](https://doi.org/10.1186/s12859-018-2129-y).
- For **multi-individual** datasets, the relevant paper to cite is:
	* Rabiee, Maryam, Erfan Sayyari, and Siavash Mirarab. 2019. “Multi-Allele Species Reconstruction Using ASTRAL.” Molecular Phylogenetics and Evolution 130 (January). 286–96. [doi:10.1016/j.ympev.2018.10.033](https://doi.org/10.1016/j.ympev.2018.10.033).
- Since version 4.10.0, ASTRAL can also compute branch length (in coalescent units) and a measure of support called **local posterior probability**, described here:
    * Sayyari, Erfan, and Siavash Mirarab. 2016. “Fast Coalescent-Based Computation of Local Branch Support from Quartet Frequencies.” Molecular Biology and Evolution 33 (7): 1654–68.  [doi:10.1093/molbev/msw079](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
- ASTRAL can also perform a **polytomy test** (`-t 10` option):
    * Sayyari, Erfan, and Siavash Mirarab. 2018. “Testing for Polytomies in Phylogenetic Species Trees Using Quartet Frequencies.” Genes 9 (3): 132. [doi:10.3390/genes9030132](http://www.mdpi.com/2073-4425/9/3/132).

- For practical tips on using ASTRAL see [this preprint book chapter](https://arxiv.org/pdf/1904.03826.pdf).

#### Papers on older versions:

- The original algorithm (ASTRAL-I) is described in:
    - Mirarab, Siavash, Rezwana Reaz, Md. Shamsuzzoha Bayzid, Théo Zimmermann, M. S. Swenson, and Tandy Warnow. 2014. “ASTRAL: Genome-Scale Coalescent-Based Species Tree Estimation.” Bioinformatics 30 (17): i541–48. [doi:10.1093/bioinformatics/btu462](doi.org/10.1093/bioinformatics/btu462).
- All the versions between 4.7.4  and 5.1.0 correspond to ASTRAL-II, described in:
    * Mirarab, Siavash, and Tandy Warnow. 2015. “ASTRAL-II: Coalescent-Based Species Tree Estimation with Many Hundreds of Taxa and Thousands of Genes.” Bioinformatics 31 (12): i44–52. [doi:10.1093/bioinformatics/btv234](http://bioinformatics.oxfordjournals.org/content/31/12/i44)


#### Papers with relevance to ASTRAL:
    
These papers do not describe features in ASTRAL, but are also relveant and we encourage you to read them:

1. **DiscoVista**: This paper shows how quartet scores (more broadly, genome discordance) can be visualized in interpretable ways. The visualization of quartet scores, in particular, is closely tied to the ASTRAL method. 
    - Sayyari, Erfan, J.B. James B. Whitfield, and Siavash Mirarab. 2018. “DiscoVista: Interpretable Visualizations of Gene Tree Discordance.” Molecular Phylogenetics and Evolution 122 (May): 110–15. [doi:10.1016/j.ympev.2018.01.019](https://doi.org/10.1016/j.ympev.2018.01.019).
- **Fragmentary data**: The following paper made the case that before inferring gene trees, removing fragmentary data (e.g., those that have uncharacteristically large numbers of gaps) should be removed. It also showed RAxML gene trees are preferable to FastTree trees. 
    - Sayyari, Erfan, James B Whitfield, and Siavash Mirarab. 2017. “Fragmentary Gene Sequences Negatively Impact Gene Tree and Species Tree Reconstruction.” Molecular Biology and Evolution 34 (12): 3279–91. [doi:10.1093/molbev/msx261](https://doi.org/10.1093/molbev/msx261).
- **Missing data**: The following paper showed that excluding genes because they have missing data is often detrimental to accuracy. 
	- Molloy, Erin K., and Tandy Warnow. 2018. “To Include or Not to Include: The Impact of Gene Filtering on Species Tree Estimation Methods.” Systematic Biology 67 (2): 285–303. [doi:10.1093/sysbio/syx077](https://doi.org/10.1093/sysbio/syx077).
- **TreeShrink**: This paper introduced a method for removing very long branches from gene trees in a statistically motivated way. These branches make gene trees less accurate. 
	-  Mai, Uyen, and Siavash Mirarab. 2018. “TreeShrink: Fast and Accurate Detection of Outlier Long Branches in Collections of Phylogenetic Trees.” BMC Genomics 19 (S5): 272. [doi:10.1186/s12864-018-4620-2](https://doi.org/10.1186/s12864-018-4620-2).
- **Sample Complexity**: This paper established the theoretical sample complexity (i.e., number of required genes) for ASTRAL.
    - Shekhar, Shubhanshu, Sebastien Roch, and Siavash Mirarab. 2018. “Species Tree Estimation Using ASTRAL: How Many Genes Are Enough?” IEEE/ACM Transactions on Computational Biology and Bioinformatics 15 (5): 1738–47. [doi:10.1109/TCBB.2017.2757930](https://doi.org/10.1109/TCBB.2017.2757930).
- **INSTRAL**: introduces an ASTRAL-based algorithm for adding new species unto an existing species tree; so, the phylogenetic placement problem but for species trees. 
	- Rabiee, Maryam, and Siavash Mirarab. 2018. “INSTRAL: Discordance-Aware Phylogenetic Placement Using Quartet Scores.” BioRxiv 432906. [doi:10.1101/432906](https://doi.org/10.1101/432906).
- **BestML:** This paper was published before ASTRAL but showed that using best ML gene trees is often preferable to using the consensus of running summary methods on bootstrapped gene trees. 
	- Mirarab, Siavash, Md Shamsuzzoha Bayzid, and Tandy Warnow. 2016. “Evaluating Summary Methods for Multilocus Species Tree Estimation in the Presence of Incomplete Lineage Sorting.” Systematic Biology 65 (3). Oxford University Press: 366–80. [doi:10.1093/sysbio/syu063](https://doi.org/10.1093/sysbio/syu063).


Documentations
-----------

- The rest of this README file
- **Our [tutorial](astral-tutorial.md)**.
- For practical tips on using ASTRAL, including on how to prepare input and interpret output, see [this paper](https://arxiv.org/pdf/1904.03826.pdf).
- The chapter of Siavash Mirarab's dissertation that describes ASTRAL in detail is provided [here](thesis-astral.pdf).
- Publications shown above have scientific details
- A [developer guide](developer-guide.md).

INSTALLATION:
-----------
* There is no installation required to run ASTRAL.
* Download using one of two approaches:
    * You simply need to download the [zip file](https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.1.zip) and extract the contents to a folder of your choice. 
    * Alternatively, you can clone the [github repository](https://github.com/smirarab/ASTRAL/). You then run `make.sh` to build the project or simply uncompress the zip file that is included with the repository.
* ASTRAL is a java-based application, and should run in any environment (Windows, Linux, Mac, etc.) as long as java is installed. 
  Java 1.5 or later is required. We have tested ASTRAL only on Linux and MAC.
* To test your installation, go to the place where you put the uncompressed ASTRAL, and run:

  ``` bash
   java -jar astral.5.7.1.jar -i test_data/song_primates.424.gene.tre
   ```

  This should quickly finish. There are also other sample input files under `test_data/` that can be used.

* ASTRAL can be run from any directory (e.g., `/path/to/astral/`). Then, you just need to run:

  ``` bash
  java -jar /path/to/astral/astral.5.7.1.jar
  ```

* Also, you can move `astral.5.7.1.jar` to any location you like and run it from there, but note that you need to move the `lib` directory with it as well.


EXECUTION:
-----------
ASTRAL currently has no GUI. You need to run it through the command-line. In a terminal, go the location where you have downloaded the software, and issue the following command:

```
  java -jar astral.5.7.1.jar
```

This will give you a list of options available in ASTRAL.

To find the species tree given a set of gene trees in a file called `in.tree`, use:

```
java -jar astral.5.7.1.jar -i in.tree
```

The results will be outputted to the standard output. To save the results in a file use the `-o` option (**Strongly recommended**):

```
java -jar astral.5.7.1.jar -i in.tree -o out.tre
```
To save the logs (**also recommended**), run:

```
java -jar astral.5.7.1.jar -i in.tree -o out.tre 2>out.log
```

###### Input: 
* The input gene trees are in the Newick format
* The input trees can have missing taxa, polytomies (unresolved branches), and also multiple individuals per species.
*  Taxon names cannot have quotation marks in their names (sorry!). This means you also cannot have weird characters like ? in the name (underscore is fine).
* When multiple individuals from the same species are available, you can ask ASTRAL to force them to be together in the species tree. To do this, a mapping file needs to be provided using the `-a` option. This mapping file should have one line per species, and each line needs to be in one of two formats:

```
species_name [number of individuals] individual_1 individual_2 ...

species_name:individual_1,individual_2,...
```
   Note that when multiple individuals exist for the same species, your species name should be different from the individual names.
   
###### Output: 
The output in is Newick format and gives: 

* the species tree topology, 
* branch lengths in coalescent units (only for internal branches or for terminal branches if that species has multiple individuals),
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1). 
* It can also annotate branches with other quantities, such as quartet support, as described in the [tutorial](astral-tutorial.md).

The ASTRAL tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if you tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py). 

### Other features (local posterior, bootstrapping):
Please refer to the [tutorial](astral-tutorial.md) for all other features, including bootstrapping, branch annotation, and local posterior probability.

### Memory:
For big datasets (say more than 5000 taxa), increasing the memory available to Java can result in speedups. Note that you should give Java only as much free memory as you have available on your machine. So, for example, if you have 3GB of free memory, you can invoke ASTRAL using the following command to make all the 3GB available to Java:

```
java -Xmx3000M -jar astral.5.7.1.jar -i in.tree
```

Acknowledgment
-----------
ASTRAL code uses bytecode and some reverse engineered code from PhyloNet package (with permission from the authors).


Bug Reports:
-----------
contact ``astral-users@googlegroups.com``
