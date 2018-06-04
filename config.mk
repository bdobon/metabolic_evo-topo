##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##                    PARAMETERS                         ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Metabolic model in MATLAB format
#MATFILE=data/erythrocyte/iAB_RBC_283.mat
MATFILE=data/Recon3D/Recon3D_301/Recon3DModel_301.mat
DATABOOST_DIR=data/hierarchical_boosting


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##                  FILE LOCATIONS                       ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Path to the scripts
PIPELINE_DIR ?= /home/bego/Documents/PROJECTS/METABOLOME/metabolic_evo-topo
SCRIPTS_DIR ?= $(PIPELINE_DIR)/src

## Output destinations
OUTPUT_DIR ?= ./results/Recon3DModel
REACTIONGRAPH_DIR = $(OUTPUT_DIR)/reactionGraph
CCOMPONENTS_DIR = $(OUTPUT_DIR)/connectedComponents
SEQUENCES_DIR = $(OUTPUT_DIR)/sequences
STATSBOOST_DIR = $(OUTPUT_DIR)/boosting
PLOTS_DIR = $(OUTPUT_DIR)/plots



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##              BINARY PROGRAMS LOCATION                 ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

SHELL=/bin/bash
SED = /bin/sed
MKDIR_P = /bin/mkdir -p 
PYTHON ?= /usr/bin/python2
RSCRIPT ?= /usr/bin/Rscript
BEDTOOLS ?= /home/bego/Software/bedtools2/bin/intersectBed


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##                  COMMAND ALIASES                      ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# NOTE: add command-line options here
# Create a directed reaction graph (DRG)
MODEL2DRG_SRC=$(SCRIPTS_DIR)/create_reaction_graph.py 
MODEL2DRG_EXE=$(PYTHON) $(MODEL2DRG_SRC) 

# Extract gene coordinates & link to reactions
COORD_SRC=$(SCRIPTS_DIR)/get_genes_coordinates.R 
COORD_EXE=$(RSCRIPT) $(COORD_SRC)

# Extract gene sequences & calculate genomic features
SEQ_SRC=$(SCRIPTS_DIR)/get_genes_sequences.R 
SEQ_EXE=$(RSCRIPT) $(SEQ_SRC)

# Calculate topology measures of all connected components of a DRG
TOPOLOGY_SRC=$(SCRIPTS_DIR)/calculate_topology_RG.py 
TOPOLOGY_EXE=$(PYTHON) $(TOPOLOGY_SRC) 

# Convert boosting files into BED
PARSEBOOST_SRC=$(SCRIPTS_DIR)/HierBoosting2BED.sh
PARSEBOOST_EXE=$(SHELL) $(PARSEBOOST_SRC)
BOOST_FILES=$(wildcard $(DATABOOST_DIR)/*.scores)
BOOSTBED_FILES=$(patsubst $(DATABOOST_DIR)/%.scores, $(DATABOOST_DIR)/%.bed, $(BOOST_FILES))

STATSBOOST_SRC=$(SCRIPTS_DIR)/selection_score_stats.R
STATSBOOST_EXE=$(RSCRIPT) $(STATSBOOST_SRC)
STATSBOOST_FILES=$(patsubst $(DATABOOST_DIR)/%.bed, $(STATSBOOST_DIR)/%.intersect, $(BOOSTBED_FILES))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




























