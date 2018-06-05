##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##                  FILE LOCATIONS                       ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Input & output destinations
## Only need to change the first two variables:
##	- MATFILE = Metabolic model in MATLAB format
##	- OUTPUT_DIR = General results folder
#MATFILE ?= data/Recon2.v04/Recon2.v04.mat
#OUTPUT_DIR ?= ./results/Recon2
#MATFILE ?= data/Recon3D/Recon3D_301/Recon3DModel_301.mat
#OUTPUT_DIR ?= ./results/Recon3DModel
MATFILE ?= data/Recon3D/Recon3D_301/Recon3D_301.mat
OUTPUT_DIR ?= ./results/Recon3D


DATABOOST_DIR ?= data/hierarchical_boosting

## Path to the scripts
PIPELINE_DIR ?= /home/bego/Documents/PROJECTS/METABOLOME/metabolic_evo-topo
SCRIPTS_DIR ?= $(PIPELINE_DIR)/src

## Output destinations
REACTIONGRAPH_DIR ?= $(OUTPUT_DIR)/reactionGraph
CCOMPONENTS_DIR ?= $(OUTPUT_DIR)/connectedComponents
SEQUENCES_DIR ?= $(OUTPUT_DIR)/sequences
STATSBOOST_DIR ?= $(OUTPUT_DIR)/boosting
PLOTS_DIR ?= $(OUTPUT_DIR)/plots


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##              BINARY PROGRAMS LOCATION                 ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

SHELL ?= /bin/bash
SED ?= /bin/sed
MKDIR_P ?= /bin/mkdir -p 
PYTHON ?= /usr/bin/python2
RSCRIPT ?= /usr/bin/Rscript
BEDTOOLS ?= /home/bego/Software/bedtools2/bin/intersectBed


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##                  COMMAND ALIASES                      ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## NOTE: add command-line options here

## Create a directed reaction graph (DRG)
MODEL2DRG_SRC=$(SCRIPTS_DIR)/create_reaction_graph.py 
MODEL2DRG_EXE=$(PYTHON) $(MODEL2DRG_SRC) 

## Extract gene coordinates & link to reactions
COORD_SRC=$(SCRIPTS_DIR)/get_genes_coordinates.R 
COORD_EXE=$(RSCRIPT) $(COORD_SRC)

## Extract gene sequences & calculate genomic features
SEQ_SRC=$(SCRIPTS_DIR)/get_genes_sequences.R 
SEQ_EXE=$(RSCRIPT) $(SEQ_SRC)

## Calculate topology measures of all connected components of a DRG
TOPOLOGY_SRC=$(SCRIPTS_DIR)/calculate_topology_RG.py 
TOPOLOGY_EXE=$(PYTHON) $(TOPOLOGY_SRC) 

## Convert boosting files into BED
PARSEBOOST_SRC=$(SCRIPTS_DIR)/HierBoosting2BED.sh
PARSEBOOST_EXE=$(SHELL) $(PARSEBOOST_SRC)
BOOST_FILES=$(wildcard $(DATABOOST_DIR)/*.scores)
BOOSTBED_FILES=$(patsubst $(DATABOOST_DIR)/%.scores, $(DATABOOST_DIR)/%.bed, $(BOOST_FILES))

## Calculate boosting stats per gene
STATSBOOST_SRC=$(SCRIPTS_DIR)/selection_score_stats.R
STATSBOOST_EXE=$(RSCRIPT) $(STATSBOOST_SRC)
STATSBOOST_FILES=$(patsubst $(DATABOOST_DIR)/%.bed, $(STATSBOOST_DIR)/%.intersect, $(BOOSTBED_FILES))

## Plot permutation postive genes boosting vs. centralities
PLOTBOOST_SRC=$(SCRIPTS_DIR)/plot_PS_permut_Boost.R
PLOTBOOST_EXE=$(RSCRIPT) $(PLOTBOOST_SRC)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



