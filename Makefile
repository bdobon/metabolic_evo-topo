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


## boostStats	:	calculate boosting stats per gene
boostStats : $(STATSBOOST_FILES) 

$(STATSBOOST_DIR)/%.intersect : $(BOOSTBED_FILES) $(REACTIONGRAPH_DIR)/gene_coordinates.bed  $(STATSBOOST_SRC)
	$(MKDIR_P) $(STATSBOOST_DIR)
	$(BEDTOOLS) -a $< -b $(REACTIONGRAPH_DIR)/gene_coordinates.bed -wb > $@
	$(STATSBOOST_EXE) $@


## clean		: Remove results folder files.
clean-all : 
	-rm -vrf $(OUTPUT_DIR)
	@echo
	@echo 'Everything is gone.'

## clean-RG	: Remove reaction graph folder.
clean-RG :
	-rm -vrf $(REACTIONGRAPH_DIR)


## variables	: Print variables.
variables :
	@echo MATFILE: $(MATFILE)
	@echo BOOSTING_DATA: $(DATABOOST_DIR)
	@echo OUTPUT: $(OUTPUT_DIR)


## help		: Show this help. 
help : Makefile
	@$(SED) -n 's/^##//p' $<









