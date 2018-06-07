## 1.Plot permutation sampling distribution of the mean boosting score/gene 
## 2.Plot 95% Confidence intervals of the maximum boosting score/gene
## do plot for every centrality measure / test / ccomponent== 000

#NOTES: ####
# make plot with ggplot or basic plot
# revise permutations p-value
# do permuts with maximum value?

## Options & functions
options('scipen'=100, 'digits'=5)
set.seed(46268008) 
require(perm)
require(ggplot2)

permutPLOT <- function(centralfile,boostfile,genes_reactions,tab_cutoff,outplot) {
  boostscore <- read.table(boostfile, header = TRUE)
  central <- read.table(centralfile, header = TRUE)

  if ( names(central)[2] != "SOURCE_SINK"){
    
  
    ### need to merge with the centrality measure & with boosting stats
    g_react_central <- merge(genes_reactions, central, by = 'REACTION')
    
    g_r_central_score <- merge(g_react_central,boostscore[,(2:4)],by = 'GENE')
    
    if ( dim(g_r_central_score)[1] > 0 ) {
      ## classify genes as positive selected or not ####
      test <- unlist(strsplit(x = basename(boostfile), split = '\\.'))[1]
      cutoff <- tab_cutoff$V2[tab_cutoff$V1 == test]
      
      g_r_central_score$select1 <- 0
      g_r_central_score$select1[g_r_central_score$maxBoost >= cutoff]  <- 1
      
      selected <- as.numeric(subset(g_r_central_score[,3], g_r_central_score$select1 == 1))
      neutral <- as.numeric(subset(g_r_central_score[,3], g_r_central_score$select1 == 0))
      
      # only calculate & plot if there are genes under positive selection
      if ( length(selected) > 0 ) {
        n_size <- length(selected)
        mean_var <- mean(selected)
        
        ## 1. Lib Permutation function
        #     methodRuleTS1(x=c(selected,neutral),group=c(rep(1,12),rep(0,length(neutral))),exact=TRUE)    
        ## to estimate best method -> exact.mc
        stat <- permTS(selected, neutral, alternative = 'two.sided', method = 'exact.mc', 
                       control = permControl(nmc = 2000, tsmethod = 'abs',setSEED = TRUE))
        
        pv2 <- stat$p.value
        
        # 2. Permutation function  
        n_permut <- 10000 # 10000 sampling 
        plotmeans <- numeric(n_permut+1)
        plotmeans[1] <- mean_var
        for (z in 2:length(plotmeans)) { 
          plotmeans[z] <- mean(sample(as.numeric(g_r_central_score[,3]), 
                                      n_size,replace = TRUE),na.rm = TRUE)      
        } # for random set 
        
        # calculate p-value ####
        datap <- data.frame(score=plotmeans)
        datap$pvalue<- rank(datap$score,ties.method='average')
        rankpos <- max(subset(datap$pvalue, datap$score<=mean_var)) # low tail
        
        pv <- (n_permut - rankpos)/n_permut
        ## putting pv2 (permutation library) in legend instead of pv (manually)
        
        # plot ####
        png(filename = paste(outplot, '/',test,'.',names(g_r_central_score)[3], '.permutHist.png', sep=''),res = 200, bg = 'white',width = 1000,height = 1000)
        text_legend <- paste('n =',n_size,'\nmean =',round(mean_var,4),
                             '\n p-value=',round(pv2,4),'\n p-permut=',round(pv,4))
        
        ## plot with ggplot2
        pretty_bin <- function(x){
          bin_val <- 2 * IQR(x) / (length(x)^(1/3))
          if (sign(bin_val) != 1 ) {
            bin_val <- diff(range(x))/30
          }
          return(bin_val)
        }
        
        g1 <- ggplot() + aes(plotmeans)+ geom_histogram(binwidth = pretty_bin(plotmeans) ) +
        geom_point(data=data.frame(x=mean_var, y=0),aes(x=x, y=0) ,size=3,col='red') +
        geom_vline(xintercept = mean_var,col='red', lty=2)  +
        annotate(geom='text',hjust = 1, vjust = 1,x =Inf, y = Inf, label=text_legend) +
        theme_classic() + xlab(names(g_r_central_score)[3]) + 
        ggtitle(test)
        print(g1)
        ## basic plotting
        # hist(plotmeans, main = test, xlab = names(g_r_central_score)[3])
        # points(c(mean_var,mean_var), c(0,0), col = 'red', pch = 21,cex=3,lwd=2) # mean VARIABLE
        # points(c(mean_var,mean_var), c(0,0), col = 'red', pch = 18,cex=2,lwd=3) # mean VARIABLE
        # legend('topright', legend= text_legend, bty = 'n')
        dev.off()
      }# if there're genes under positive selection
    } # if there're genes with scores/centralities
  } # don't do anything with categorical variables
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
  ##   data: a data frame.
  ##   measurevar: the name of a column that contains the variable to be summariezed
  ##   groupvars: a vector containing names of columns that contain grouping variables
  ##   na.rm: a boolean that indicates whether to ignore NA's
  ##   conf.interval: the percent range of the confidence interval (default is 95%)
  
  
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

ciPLOT <- function(centralfile,boostfile,genes_reactions,outplot) {
  boostscore <- read.table(boostfile, header = TRUE)
  central <- read.table(centralfile, header = TRUE)
  test <- unlist(strsplit(x = basename(boostfile), split = '\\.'))[1]
  
  ### need to merge with the centrality measure & with boosting stats
  g_react_central <- merge(genes_reactions, central, by = 'REACTION')
  
  g_r_central_score <- merge(g_react_central,boostscore[,(2:4)],by = 'GENE')
  
  if ( names(g_r_central_score)[3] == "SOURCE_SINK"){
    g_r_central_score$central <- factor(g_r_central_score[,3], levels = c("source", "intermediate", "sink"))
  }else{
    # classify central measure by quantiles  
    qs <- quantile(g_r_central_score[,3],probs = seq(0.25,0.75,0.25), type=7, na.rm = TRUE,names = FALSE)
    # classify central measure by them
    g_r_central_score$central[g_r_central_score[,3] < qs[1]] <- "<25%"
    g_r_central_score$central[g_r_central_score[,3] >= qs[1] & g_r_central_score[,3] < qs[2] ] <- "25-50%"
    g_r_central_score$central[g_r_central_score[,3] >= qs[2] & g_r_central_score[,3] < qs[3] ] <- "50-75%"
    g_r_central_score$central[g_r_central_score[,3] >= qs[3]] <- ">75%"
  }
  # summarySE provides the standard deviation, standard error of the mean, and a 95% confidence interval
  tgc <- summarySE(g_r_central_score, measurevar="maxBoost", groupvars="central",na.rm = T)
  ## Plot #### 
  # LINE GRAPH: error bars represent the 95% CI
  png(filename = paste(outplot, '/',test,'.',names(g_r_central_score)[3], '.95CI.png', sep=''),res = 200, bg = 'white',width = 1000,height = 1000)
  
  p1 <- ggplot(tgc, aes(x=central, y=maxBoost, group=central)) +
    geom_errorbar(width=.1, aes(ymin=maxBoost-ci, ymax=maxBoost+ci)) +
    geom_point(shape=21, size=3, fill="white")+
    theme_classic()+ # theme(panel.background = element_rect(fill = '#dee8ed'))
    xlab(names(g_r_central_score)[3]) +
    ylab(test)
  print(p1)
  dev.off()
}


## inputs ####
args <- commandArgs(trailingOnly = TRUE)
output <- args[1]
table_cutoff <- read.table('data/hierarchical_boosting/thresholds_of_significance.txt', header = F)

## need to link GENE with REACTION
genes_reactions_link <- read.table(paste(dirname(output),'/reactionGraph/geneReactions.list',sep=''), header = T)

## boosting scores
boostingfiles <- list.files(paste(dirname(output),'/boosting' ,sep=''), pattern = '*.txt', full.names = TRUE)


## do this for every connected component
## most have no edges or genes under positive selection
components_folders <- list.dirs(paste(output,sep=''), full.names = TRUE, recursive = F)[1]

for (ccomp in components_folders) {
  
  if (dir.exists(paste(ccomp,'/topology', sep=""))){ 
    centralfiles <- list.files(paste(ccomp,'/topology' ,sep=''), pattern = '*.list',full.names = TRUE)

    cat(c('Plotting component:',basename(ccomp),'\n'))
    output_folder <- paste(dirname(output),'/plots/',basename(ccomp), sep='')
    dir.create(output_folder,showWarnings = F)
    
    for (bfile in boostingfiles) {
      for (cfile in centralfiles) {
        permutPLOT(cfile,bfile,genes_reactions_link,table_cutoff,output_folder)
        ciPLOT(cfile,bfile,genes_reactions_link,output_folder)
      } ## for all centralities calculated
    }    ## for all boosting tests
  }else{
    cat(c('No topology for component',basename(ccomp)),'\n')
  } # if component has topology 
} 


