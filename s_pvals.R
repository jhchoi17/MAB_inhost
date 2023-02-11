library(dplyr)
library(ggplot2)
library(coin)
library(FSA)

iterations=10000 ## total number of permutations tested; manually change this 

##########formatting actual distribution of mutations into gene: count: snp: indel

mutfile <- read.csv('path/to/list/of/actual/mutations') ##this file has columsn strain: type (indel/SNP): gene: count

mutfile <- na.omit(mutfile) 

mutfile$gene <- gsub('hypotheticalprotein', '\\1', mutfile$gene)
mutfile$gene <- gsub('intron', '\\1', mutfile$gene)
mutfile$gene <- gsub('unknownnode', '\\1', mutfile$gene)
mutfile <- mutfile[mutfile$gene!='',]


##optional:get rid of randoma allele numbers 
mutfile$gene <- gsub('_10','\\1',mutfile$gene)
mutfile$gene <- gsub('_11','\\1',mutfile$gene)
mutfile$gene <- gsub('_12','\\1',mutfile$gene)
mutfile$gene <- gsub('_13','\\1',mutfile$gene)
mutfile$gene <- gsub('_1','\\1',mutfile$gene)
mutfile$gene <- gsub('_2','\\1',mutfile$gene)
mutfile$gene <- gsub('_3','\\1',mutfile$gene)
mutfile$gene <- gsub('_4','\\1',mutfile$gene)
mutfile$gene <- gsub('_5','\\1',mutfile$gene)
mutfile$gene <- gsub('_6','\\1',mutfile$gene)
mutfile$gene <- gsub('_7','\\1',mutfile$gene)
mutfile$gene <- gsub('_8','\\1',mutfile$gene)
mutfile$gene <- gsub('_9','\\1',mutfile$gene)


mutfile_summ <- mutfile %>% group_by(gene) %>% summarize(n()) # how many strains per gene 

genelist <- unique(mutfile$gene) 
n <- length(genelist) 

counts <- c()
snps <- c()
indels <- c()


for (i in 1:n){
  temp <- as.character(genelist[i])
  m <- length(unique(mutfile[mutfile$gene==temp,]$strain)) #how many strains total? 
  counts <- c(counts, m)
  k <- length(unique(mutfile[mutfile$gene==temp & mutfile$type=='SNP',]$strain)) #which are snp ?
  l <- length(unique(mutfile[mutfile$gene==temp & mutfile$type=='indel',]$strain)) #which are indel?
  
  snps <- c(snps, k)
  indels <- c(indels, l)
}

realmuts <- data.frame(gene=genelist, numstrain=counts, SNP=snps,indel=indels )
actual <- realmuts ##Actual distribution data frame



#############################format permutation results 

reference <- read.csv('10000randommuts_results.csv') ##load results file from s_randommuts.py
colnames(reference) <- c('iter', 'gene','strain' )


#collapse alleles: optional (based on genome annotation method)
reference$gene <- gsub('_10','',reference$gene)
reference$gene <- gsub('_11','',reference$gene)
reference$gene <- gsub('_12','',reference$gene)
reference$gene <- gsub('_13','',reference$gene)
reference$gene <- gsub('_14','',reference$gene)
reference$gene <- gsub('_15','',reference$gene)
reference$gene <- gsub('_16','',reference$gene)
reference$gene <- gsub('_17','',reference$gene)
reference$gene <- gsub('_18','',reference$gene)

reference$gene <- gsub('_1','',reference$gene)
reference$gene <- gsub('_2','',reference$gene)
reference$gene <- gsub('_3','',reference$gene)
reference$gene <- gsub('_4','',reference$gene)
reference$gene <- gsub('_5','',reference$gene)
reference$gene <- gsub('_6','',reference$gene)
reference$gene <- gsub('_7','',reference$gene)
reference$gene <- gsub('_8','',reference$gene)
reference$gene <- gsub('_9','',reference$gene)


####generate p-value

genelist <- unique(actual$gene) 
k <- length(genelist)

actual$pval <- c('.') ## add new column for pval 

for (i in 1:k){
  target <- as.character(genelist[i]) ## target=gene 
  real <- as.numeric(as.character(actual[actual$gene==target,]$numstrain[1])) #actual number of strains mut in parallel 
  if (target %in% reference$gene){
    tempdf <- reference[reference$gene==target,]

    tempdf$strain <- as.numeric(tempdf$strain)
    tempdf <- aggregate(tempdf$strain, by=list(iter=tempdf$iter), FUN=sum) ##Aggregate separate entries for alleles so as not to double-count
    colnames(tempdf) <- c('iter', 'strain')
    n <- length(rownames(tempdf))
    diff <- iterations-n
    
    iter <- rep('.', diff)
    strains <- rep(0, diff)
    
    tempdf2 <- data.frame(iter=iter, strain=strains)
    tempdf <- rbind(tempdf, tempdf2)
    
    percentile <- ecdf(tempdf$strain)
    pvalue <- 1-percentile(real)
    
    actual[actual$gene==target,]$pval <- pvalue 
  }
}


actual$pval <- unlist(actual$pval, use.names=FALSE)

write.csv(actual, '10k_pvals_no_alleles.csv') ##edit file naame as needed

