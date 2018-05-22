#
# Makefile for Evolutionary & Topologycal Analysis of a Metabolic Network
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include config.mk

.PHONY: all help clean-all variables model2DRG topology duplicates boost2BED

## all		: reaction graph,topology,selection,correlations,plots,archive.
all : model2DRG topology 


## model2DRG	: create a directed reaction graph from a model
model2DRG : $(REACTIONGRAPH_DIR)


$(REACTIONGRAPH_DIR) : $(MATFILE) $(MODEL2DRG_SRC) $(COORD_SRC) 
	$(MKDIR_P) $@
	$(MODEL2DRG_EXE) $(REACTIONGRAPH_DIR) $< 
	$(COORD_EXE) $@


## topology	: calculate topology measures by connected components
topology :	$(CCOMPONENTS_DIR)

$(CCOMPONENTS_DIR) : $(REACTIONGRAPH_DIR)
	$(MKDIR_P) $@
	$(TOPOLOGY_EXE) $< $(CCOMPONENTS_DIR)


## duplicates	: manually add duplicates/missing gene IDs
duplicates :	
	$(COORD_EXE) $(REACTIONGRAPH_DIR) 'addDUP'


## getBoosting	: convert boosting data into BED
getBoosting : $(BED_FILES)

$(DATABOOST_DIR)/%.bed : $(DATABOOST_DIR)/%.scores $(PARSEBOOST_SRC) $(REACTIONGRAPH_DIR)
	@$(PARSEBOOST_EXE) $< $@
	$(BEDTOOLS) -a $@.gz -b $(REACTIONGRAPH_DIR)/gene_coordinates.bed  -wb > $@.genes
	

## clean		: Remove results folder files.
clean-all: 
	-rm -vrf $(OUTPUT_DIR)
	@echo
	@echo 'Everything is gone.'

## clean-RG	: Remove reaction graph folder.
clean-RG:
	-rm -vrf $(REACTIONGRAPH_DIR)


## variables	: Print variables.
variables:
	@echo OUTPUT: $(OUTPUT_DIR)
	@echo MATFILE: $(MATFILE)
	@echo BOOSTING_DATA: $(DATABOOST_DIR)

## help		: Show this help 
help : Makefile
	@$(SED) -n 's/^##//p' $<









