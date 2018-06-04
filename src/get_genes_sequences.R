## Extract the CDS for a gene to calculate:
## - CDS length 
## - Effective number of codons with Codonw
## - GC content
## To use as a genomic features to correct for in the
## general linear model 


## functions ####
library(biomaRt)
library(foreach)
library(doParallel)
library(seqinr)
registerDoParallel(cores=detectCores()-1)

retrieveSequence <- function(gene, outpath) {
  getLongestSeq <- function(x) {
    ## if there are more than one transcript with GENCODE flag 
    lenSeq <- function(s) { length(unlist(strsplit(s,split = ""))) }
    tab <- sapply(x$coding, lenSeq)
    longest <- x[x$coding == names(tab[tab == max(tab)]),]
    return(longest[1,])
  }
  
  ## do these steps for every gene to handle exceptions
  # VALUES: set restrictions & data input (list of vectors == filters)
  val <- list('protein_coding',  c(1:22, 'X', 'Y', 'MT') , gene, TRUE)
  q <- try(getBM(attributes = att, filters= filt, values = val , mart = ensembl),silent = TRUE)
  if ( class(q) != "try-error") {
    q1 <- subset(q, q$transcript_gencode_basic == 1)
    if (dim(q1)[1] > 1){
      q2 <- getLongestSeq(q1)
      print(paste("Retrieved:",gene))
      write.fasta(sequences = q2['coding'],
                  names = paste(q2['hgnc_symbol'],q2['entrezgene'],q2['ensembl_gene_id'],sep = " | "), 
                  file.out = paste(outpath, '/',q2['hgnc_symbol'],'_human.fasta',sep="" ) ) 
    }
  }else{
    print(paste("Cannot retrieve:",gene))
    write.table(gene, paste(outpath,'/missing_sequences.txt', sep = ""),
                append = TRUE, row.names = F,col.names = F, quote = F)
  }
}


## MAIN ####
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]
gene_coord <- read.table(input, header = T, stringsAsFactors = F)

# Select BioMart database & dataset: GRCh37
ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl', GRCh=37)

# ATTRIBUTES: define the values the user is interested in
att <- c('entrezgene','ensembl_gene_id','hgnc_symbol', 'coding', 'transcript_gencode_basic')
filt <- c('biotype', 'chromosome_name', 'hgnc_symbol', 'with_ccds')

# remove file because it appends to it
if (file.exists(paste(output,'/missing_sequences.txt', sep = ""))){
  file.remove(paste(output,'/missing_sequences.txt', sep = ""))
}

# retrieve sequences
foreach(g = gene_coord$hgnc_symbol) %dopar% {
  retrieveSequence(gene = g, outpath=output)
}

