#!/usr/bin/env Rscript
source("/home/syjoo/Script/CNV/4.GENE_sub.r")
PATH <- "/home/syjoo/Project/CNV/Batch15/SPLIT/"
mode <- c("deletion") #choose between "deletion" or "duplication" depending on your mode of analysis
path <- paste(PATH,"del_gene",sep='') # data directory for aggregation 
final(mode=mode,path=path)
