# Topological and Evolutionary Analysis of a Metabolic Network

Perform a topological and evolutionary analysis of a metabolic network: 
creates a directed reaction graph from a metabolic model (MATLAB format).


### Prerequisites

For running in Ubuntu you need to install:

* Python 2.7.12
* R 3.4.1
* awk
* gzip
* bedtools intersectBed http://bedtools.readthedocs.io/en/latest/index.html

#### Python modules (install with pip install module):

* cobrapy https://github.com/opencobra/cobrapy
* networkx
* pandas


#### R libraries:

* biomaRt
* ggplot2

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



## Workflow






## Authors

* **Begona Dobon** 
