#!/bin/bash

# Format hierarchical boosting scores files already calculated in 1000GP data to use with bedtools intersect.
# Threshold of significance based on simulations: thresholds_of_significance.POP.txt
#
# Downloaded from:
# http://hsb.upf.edu/hsb_data/hierarchical_boosting/selection_scan/
#
# Input file format is:
# gridID chromosome start end score
# grid_3775000 1 3750000 3775000 -0.052914489027694
#
# Output file format is:
#   chromosome start end score
# chr1 3750000 3775000 -0.052914489027694
#
# Use only:
# whole_genome.POP.Complete.boosting.scores
# whole_genome.POP.Incomplete.boosting.scores  


ifile=$1
bedfile=$2

echo "Processing $ifile"

## skip header
awk 'FNR >1  {print "chr"$2"\t"$3"\t"$4"\t"$5}' $ifile > $bedfile
#gzip -f $bedfile


