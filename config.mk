##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##                    PARAMETERS                         ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Metabolic model in MATLAB format
MATFILE=data/erythrocyte/iAB_RBC_283.mat

## Human Build GRCh37 ##**add this to src/get_genes_coordinates.R
HUMAN_BUILD = GRCh37
ADD_DUPS ?=FALSE ##**add this to src/get_genes_coordinates.R as an option

## Send output to STDOUT (leave empty) or to a file.
PROG_OUT = #&>/dev/null ##**check if this works or add


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##                                                       ##
##                  FILE LOCATIONS                       ##
##                                                       ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Path to the scripts
PIPELINE_DIR ?= /home/bego/Documents/PROJECTS/metabolic_evo-topo
SCRIPTS_DIR ?= $(PIPELINE_DIR)/src

## Output destinations
OUTPUT_DIR ?= ./results
REACTIONGRAPH_DIR = $(OUTPUT_DIR)/reaction_graph
CCOMPONENTS_DIR = $(OUTPUT_DIR)/cComponents
BOOSTING_DIR = $(OUTPUT_DIR)/boosting
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
COORD_EXE=$(RSCRIPT_PATH) $(COORD_SRC)

# Calculate topology measures of all connected components of a DRG
TOPOLOGY_SRC=$(SCRIPTS_DIR)/calculate_topology_RG.py 
TOPOLOGY_EXE=$(PYTHON) $(TOPOLOGY_SRC) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

