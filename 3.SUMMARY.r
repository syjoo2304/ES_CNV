#!/usr/bin/env Rscript

setwd("/home/syjoo/Script/CNV")
source('3.SUMMARY_sub.r')
pd='~/Project/CNV/Batch15/'
over.f(pd)
ind.f(pd)
comp.f(pd)
summary.f(pd)
final.f(pd)
gene_anno.f(pd,ref='/home/DATA/REF',del.dir='/~/Project/CNV/Batch15/SPLIT/del_gene',dup.dir='~/Project/CNV/Batch15/SPLIT/dup_gene')

