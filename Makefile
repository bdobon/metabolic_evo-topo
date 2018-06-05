#
# Makefile for Evolutionary & Topologycal Analysis of a Metabolic Network
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include config.mk

.PHONY: all help clean-all variables model2DRG topology duplicates parseBoosting boostStats getSequences

## all		: reaction graph,topology,selection,correlations,plots,archive.
all : parseBoosting model2DRG topology boostStats


## model2DRG	: create a directed reaction graph from a model.
model2DRG : $(REACTIONGRAPH_DIR)

$(REACTIONGRAPH_DIR) : $(MATFILE) $(MODEL2DRG_SRC) $(COORD_SRC) 
	$(MKDIR_P) $@
	$(MODEL2DRG_EXE) $@ $< 
	$(COORD_EXE) $@


## topology	: calculate topology measures by connected components.
topology :	$(CCOMPONENTS_DIR)

$(CCOMPONENTS_DIR) : $(REACTIONGRAPH_DIR)
	$(MKDIR_P) $@
	$(TOPOLOGY_EXE) $< $(CCOMPONENTS_DIR)


## duplicates	: manually add duplicates/missing gene IDs.
duplicates :	
	$(COORD_EXE) $(REACTIONGRAPH_DIR) 'addDUP'


## getSequences	: retrieve fasta sequences & genomic features.
getSequences :	$(SEQUENCES_DIR)

$(SEQUENCES_DIR) : $(REACTIONGRAPH_DIR)/gene_coordinates.txt $(SEQ_SRC)
	$(MKDIR_P) $@
	$(SEQ_EXE) $< $(SEQUENCES_DIR)


## parseBoosting	: convert boosting data into BED.
parseBoosting : $(BOOSTBED_FILES)

$(DATABOOST_DIR)/%.bed : $(DATABOOST_DIR)/%.scores $(PARSEBOOST_SRC)
	@$(PARSEBOOST_EXE) $< $@


## boostStats	: calculate boosting stats per gene.
boostStats : $(STATSBOOST_FILES) 

$(STATSBOOST_DIR)/%.intersect : $(DATABOOST_DIR)/%.bed $(REACTIONGRAPH_DIR)/gene_coordinates.bed  $(STATSBOOST_SRC) 
	$(MKDIR_P) $(STATSBOOST_DIR)
	$(BEDTOOLS) -a $< -b $(REACTIONGRAPH_DIR)/gene_coordinates.bed -wb > $@
	$(STATSBOOST_EXE) $@

## plotBoostTopo	: plot permutations boosting vs. topology 
plotBoostTopo :	$(PLOTS_DIR)

$(PLOTS_DIR) : $(CCOMPONENTS_DIR) $(REACTIONGRAPH_DIR)/geneReactions.list $(STATSBOOST_FILES) $(PLOTBOOST_SRC)
	$(MKDIR_P) $@
	$(PLOTBOOST_EXE) $< 


## clean		: remove results folder files.
clean-all : 
	-rm -vrf $(OUTPUT_DIR)
	@echo
	@echo 'Everything is gone.'

## clean-RG	: remove reaction graph folder.
clean-RG :
	-rm -vrf $(REACTIONGRAPH_DIR)

## clean-topo	: remove topology measures folder.
clean-topo :
	-rm -vrf $(CCOMPONENTS_DIR)

## clean-plots	: remove plots folder.
clean-plots :
	-rm -vrf $(PLOTS_DIR)


## variables	: print variables.
variables :
	@echo MATFILE: $(MATFILE)
	@echo BOOSTING_DATA: $(DATABOOST_DIR)
	@echo OUTPUT: $(OUTPUT_DIR)


## help		: show this help. 
help : Makefile
	@$(SED) -n 's/^##//p' $<









