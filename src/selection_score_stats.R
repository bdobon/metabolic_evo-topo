## Calculate mean and maximum boosting score by gene

# inputs ####
args <- commandArgs(trailingOnly = TRUE)

ifile <- args[1]
dat <- read.table(ifile)

filename <- gsub(replacement = '', pattern = '.*whole_genome.', ifile)
test <- paste(unlist(strsplit(x = filename, split = '\\.'))[1:2],collapse = '_')
output <- dirname(ifile)

table_cutoff <- read.table('data/hierarchical_boosting/thresholds_of_significance.txt', header = F)

# main ####

genes_link <- unique(dat[,c('V8', 'V9')])

mean_boost <- aggregate(dat$V4,by = list(dat$V9),FUN = mean)
max_boost <- aggregate(dat$V4,by = list(dat$V9),FUN = max)

merged <- merge(mean_boost, max_boost, by='Group.1')
merged2 <- merge(genes_link, merged, by.x = 'V9', by.y = 'Group.1')

names(merged2) <- c('hgnc_symbol','GENE', 'meanBoost', 'maxBoost')

write.table(merged2, paste(output,'/',test ,'.txt', sep=''), quote = F, sep = '\t', row.names =F)
