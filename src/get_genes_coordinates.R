## From a list of EntrezGene ids:
## Obtain gene coordinates in GRCh37 
## gene symbols and ensembl ids 
## Requires internet connection to connect to Ensmbl


## If addDUP set to TRUE 
## Create manually 2 files with the same format as genes_coordinates.txt:
## - dup_gene_coordinates2.txt
## - missing_ids2.txt

## functions ####
generateQuery <- function(list_genes, id_genes='EntrezGene', output) {
  require(biomaRt)
  # Select BioMart database & dataset: GRCh37
  ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl', GRCh=37)
  # ATTRIBUTES: define the values the user is interested in
  att <- c('entrezgene','ensembl_gene_id','hgnc_symbol', 'chromosome_name', 'start_position',
           'end_position', 'strand','transcript_gencode_basic')
  # FILTERS: restrictions (vector)
  # Extract only protein_coding genes, NO appris data, use GENCODE basic flag
  # by the geneID
  if (id_genes == 'Symbol'){
    gene_id <-'hgnc_symbol'
  }
  if (id_genes == 'Ensembl'){
    gene_id <- 'ensembl_gene_id'
  }
  if (id_genes == 'EntrezGene'){
    gene_id <- 'entrezgene'
  }else{
    print('geneID = EntrezGene | Symbol | define your own!')
  }

  filt <- c('biotype', 'chromosome_name', gene_id)
  # VALUES: set restrictions & data input (list of vectors == filters)
  val <- list('protein_coding',  c(1:22, 'X', 'Y', 'MT') , unique(list_genes))
  q <- getBM(attributes = att, filters= filt, values = val , mart = ensembl)
  q1 <- subset(q, q$transcript_gencode_basic == 1)
  q1$transcript_gencode_basic <- NULL

  ## Check the quality of the query result: number of retrieved genes differs
  print(paste('Genes:' ,length(list_genes)))
  ids_duplicated <- q1$entrezgene[which(duplicated(q1$entrezgene))]

  ## Number of unique ids
  ids_unique <-  q1[!q1$entrezgene %in% ids_duplicated,]

  ## Number of duplicated ids by entrez
  genes_duplicated <- q1[q1$entrezgene %in% ids_duplicated,]
  ## Remove duplicated ids without associated HGNC symbol
  nonempty_genes <- genes_duplicated[which(genes_duplicated$hgnc_symbol != ''),]
  ids_duplicated2 <- nonempty_genes$entrezgene[which(duplicated(nonempty_genes$entrezgene))]
  genes_duplicated2 <- nonempty_genes[nonempty_genes$entrezgene %in% ids_duplicated2,]
  ## Add the new unique ids
  ids_unique2 <-  nonempty_genes[!nonempty_genes$entrezgene %in% ids_duplicated2,]
  ids_unique <- rbind(ids_unique,ids_unique2)

  ## Number of duplicate ids by gene Symbol
  ids_duplicated3 <- ids_unique$hgnc_symbol[which(duplicated(ids_unique$hgnc_symbol))]
  genes_duplicated3 <- ids_unique[ids_unique$hgnc_symbol %in% ids_duplicated3,]

  ## Add to the other dups
  genes_duplicated2 <- rbind(genes_duplicated2, genes_duplicated3)
  print(paste('Duplicated ids to manually select:', length(unique(genes_duplicated2$entrezgene) )))
  ## remove these last dups from the unique ids
  ids_unique <-  ids_unique[!ids_unique$hgnc_symbol %in% ids_duplicated3,]
  print(paste('Total unique ids:', dim(ids_unique)[1]))
  write.table(ids_unique, paste(output,'/gene_coordinates.txt',sep=''),
              row.names = F, sep = '\t', quote = F)

  ## These are the ids to look for manually --> create dup_gene_coordinates2.txt
  write.table(genes_duplicated2, paste(output,'/dup_gene_coordinates.txt',sep=''),
              row.names = F, sep = '\t', quote = F)

  ## Missing ids to look for manually in the duplicated ids --> create missing_ids2.txt
  missing_ids <- setdiff(list_genes,
                         c(unique(genes_duplicated2$entrezgene),ids_unique$entrezgene))
  write.table(missing_ids, paste(output,'/missing_ids.list',sep=''),
              row.names = F, col.names = F, quote = F, sep = '\t')
  print(paste('Missing ids (if you manually select the dups):',
              length(missing_ids)))
  return(ids_unique)
}
coord2BED <- function(ifile, output) {
  # Create files for IntersectBed with genomic coordinates for each gene:
  # (start, end) and 10 kb up/downstream
  # Transform:  gene_coordinates.txt
  # entrezgene	ensembl_gene_id	hgnc_symbol	chromosome_name	start_position	end_position	strand
  # 100	ENSG00000196839	ADA	20	43248163	43280874	-1
  # Into this: gene_coordinates.bed
  # chr1	1157628	1180421	126792	B3GALT6
  # chr1	1672670	1721896	65220	NADK
  # NOTE: BED zero-based but Ensembl is one-based
  coord <- ifile[,c("chromosome_name", "start_position", "end_position", "entrezgene","hgnc_symbol")]
  coord <- coord[order(coord$chromosome_name, coord$start_position),]
  coord$chromosome_name <- paste("chr",coord$chromosome_name,sep = "")
  ## substract 1 to the start position: BED zero-based but Ensembl is one-based!! the way they're are grabbed is x____](____]
  ## add 10kb up/downstream
  coord$start_position <- coord$start_position - 10001
  coord$end_position <- coord$end_position + 10000
  ## make sure there are no negative coordinates

  ## NOTE: maybe check that they don't go outside chromosomes coordinates either**

  coord$start_position[which(coord$start_position < 0)] <- 0
  write.table(coord, paste(output,'/gene_coordinates.bed',sep=''),
              row.names = F, col.names = F,sep = '\t', quote = F)
}
gene2reaction <- function(ifile, output) {
  ## Make link between all gene IDs and REACTION
  geneIDs <- ifile[,c('entrezgene', 'ensembl_gene_id', 'hgnc_symbol' )]
  gene_reaction <- read.table(paste(output,'/geneReactions.list',sep=''), header = T)
  geneIDs_reaction <- merge(gene_reaction, geneIDs, by.x = 'GENE', by.y = 'entrezgene')
  write.table(geneIDs_reaction, paste(output,'/geneIDs_Reactions.txt',sep=''),
              row.names = F, sep = '\t', quote = F)
}
manually_addDUP <- function(list_genes,ifile, output) {
  ## Create manually 2 files with the same format as genes_coordinates.txt:
  ## - dup_gene_coordinates2.txt
  ## - missing_ids2.txt
  if (file.exists(paste(output,'/dup_gene_coordinates2.txt', sep=''))){
    ids_dups <- read.table(paste(output,'/dup_gene_coordinates2.txt', sep=''), header = T)
    print(paste('dups added:',dim(ids_dups)[1]))
    ifile <- rbind(ifile, ids_dups)
  }else{
    print('No duplicated ids recovered')
  }
  if (file.exists(paste(output,'/missing_ids2.txt', sep=''))){
    ids_miss <- read.table(paste(output,'/missing_ids2.txt', sep=''), header = T)
    print(paste('missing added:',dim(ids_miss)[1]))
    ifile <- rbind(ifile, ids_miss)
  }else{
    print('No missing ids recovered')
  }
  print(paste('Genes:' ,length(list_genes)))
  print(paste('Total ids recovered:', dim(ifile)[1] ))
  file.copy(from = paste(output,'/gene_coordinates.txt',sep=''),
            to = paste(output,'/gene_coordinates_ori.txt',sep=''))
  write.table(ifile, paste(output,'/gene_coordinates.txt',sep=''),
              row.names = F, sep = '\t', quote = F)
}


## MAIN ####

args <- commandArgs(trailingOnly = TRUE)
output <- args[1]
addDUP <- args[2] 

genelist <- unique(read.table(paste(output,'/gene.list',sep=''), header = F, colClasses='character')$V1)

if (isTRUE(addDUP)) {
  gene_coord_ori <- read.table(paste(output,'/gene_coordinates.txt', sep=''), header = T)
  manually_addDUP(genelist, gene_coord_ori,output)
  ## transform to BED file
  coord2BED(gene_coord_ori, output)
  ## link gene IDs with REACTION
  gene2reaction(gene_coord_ori, output)
}else{
  ## Query information from Biomart
  gene_info <- generateQuery(genelist, id_genes ='EntrezGene', output)
  ## transform to BED file
  coord2BED(gene_info, output)
  ## link gene IDs with REACTION
  gene2reaction(gene_info, output)
}
