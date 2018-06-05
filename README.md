# Topological and Evolutionary Analysis of a Metabolic Network

Perform a topological and evolutionary analysis of a metabolic network: 
creates a directed reaction graph from a metabolic model (MATLAB format).


### Prerequisites

For running in Ubuntu you need to install:

* Python 2.7.12 (Packages: cobrapy, networkx, pandas)
* R 3.4.1 (Packages: biomaRt, seqinr, ggplot2)
* awk
* gzip
* bedtools intersectBed http://bedtools.readthedocs.io/en/latest/index.html
* CodonW https://sourceforge.net/projects/codonw/


Paths to binary programs and input/output directories are set in `config.mk`

### Testing 

```
Give the example running with Erythrocyte model

```

## Description

### model2DRG	

Create a directed reaction graph from a metabolic model.
If the products of a node (REACTION) are the reactants of any other node,
it creates a directed link between them.
It removes currency metabolites and allows self-loops.

Returns a series of files with information about the graph:

* List of nodes (REACTIONS)
* List of edges (DIRECTED)
* List of genes (EntrezGene IDs)
* List of subsystems (METABOLIC PATHWAYS)
* Link between genes IDs and reactions
* Gene coordinates in GRCh37 (requires internet connection to connect to Ensembl).


### topology	

Calculate topology measures by connected components.
It calculates: in-degree, out-degree, degree, closeness, betweenness, ratio in/out-degree, 
source/sink, predecessors, successors.

### parseBoosting	

Format hierarchical boosting scores files already calculated in 1000GP data to use with bedtools intersect.

Added file of threshold of significance based on simulations: thresholds_of_significance.POP.txt

Downloaded from:
http://hsb.upf.edu/hsb_data/hierarchical_boosting/selection_scan/

### boostStats

Calculate maximum and average value for each boosting score per gene (+10kb up/downstream). If at least one window reaches the threshold of significance classifies the gene as being under positive selection. 

### getSequences

Retrieve CDS (coding sequence) in FASTA format for all genes. 

## Workflow

```

make parseBoosting	# already run 

make model2DRG	

make topology	

make boostStats

make getSequences


```

Alternatively `make all` run all the previous steps. 


## Authors

* **Begona Dobon** 
