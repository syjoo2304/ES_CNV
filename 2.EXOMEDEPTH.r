#! /usr/bin/env Rscript
exomedepth.f <- function(path.case, path.out){
  ## Loading packages
  library(ExomeDepth)
  library(Rsamtools)
  library(GenomicAlignments)
  library(data.table)

  dir.create(path.out)
  ## getBamCounts
  target.genes <- read.table("/home/DATA/REF/Agilent_SureSelect_V5_covered.bed",header=F,stringsAsFactors=F)
  colnames(target.genes) <- c("chromosome","start","end","name")

  target.genes.GRanges <- GRanges(seqnames = target.genes$chromosome,
                                  IRanges(start=target.genes$start,end=target.genes$end),
                                  names = target.genes$name)

  bam_list.case <- list.files(path = path.case, pattern = ".bam$", full.names = TRUE)
  path.cont = "/home/DATA/YUHL/WES/CTRL/FASTQ/FASTQ_link/BAM"
  bam_list.cont <- list.files(path = path.cont, pattern = ".bam$", full.names = TRUE)
  bam_list.cont = bam_list.cont[1:10]
  bam_list <- c(bam_list.case, bam_list.cont)

  my.count <- getBamCounts(bed.frame = target.genes, bam.files= bam_list,
                           include.chr = F, referenceFasta = "/home/program/EXCAVATOR_Package_v2.2/EXCAVATOR/ucsc.hg19.fasta")

  MyCount.dafr <- as(my.count[,  colnames(my.count)], 'data.frame')

  ## remove the annoying chr letters
  MyCount.dafr$chromosome <- gsub(as.character(MyCount.dafr$chromosome), pattern = 'chr', replacement = '')
  write.table(MyCount.dafr, paste(path.out, "/MyCount_dafr.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

  ## prepare the main matrix of read count data
  MyCount.mat <- as.matrix(MyCount.dafr[, grep(names(MyCount.dafr), pattern = '*.bam')])
  nsamples <- ncol(MyCount.mat)
  write.table(MyCount.mat, paste(path.out, "/MyCount_mat.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

  ## common controls
  case.list <- colnames(MyCount.mat)[1:length(bam_list.case)]
  control.list <- colnames(MyCount.mat)[c(length(bam_list.case)+1):length(bam_list)]

  for (i in case.list){
    ## (1) Create the aggregate reference set for this sample
    ## control check
    my.choice <- select.reference.set (test.counts =  MyCount.mat[,i], reference.counts = MyCount.mat[,control.list],
                                       bin.length = (MyCount.dafr$end - MyCount.dafr$start)/1000, n.bins.reduced = 10000)
    my.reference.selected <- apply(X = MyCount.mat[, my.choice$reference.choice, drop = FALSE], MAR = 1, FUN = sum)
    write.table(cbind(MyCount.mat[,i],my.reference.selected),
                paste(path.out, "/MyCount_mat_", unlist(strsplit(i,split=".bam"))[1], ".txt",sep=""),
                col.names= c(i, paste(my.choice$reference.choice,collapse="_")),
                row.names=F,quote=F,sep="\t")
  }

  res_list <- list.files(path = path.out, pattern = "MyCount_mat_", full.names = TRUE)

  for(i in 1:length(res_list)){
    tmp.res <- read.table(res_list[i],header=T,stringsAsFactors=F)
    tmp.cont <- colnames(tmp.res)[2]
    if(i == 1){
      write.table(t(c(colnames(tmp.res[1]), unlist(strsplit(unlist(strsplit(tmp.cont,split=".bam_")),split=".bam")) )),
                  paste(path.out, "/summary_MyCount_mat.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
    }
    if(i != 1){
      write.table(t(c(colnames(tmp.res[1]), unlist(strsplit(unlist(strsplit(tmp.cont,split=".bam_")),split=".bam")) )),
                  paste(path.out, "/summary_MyCount_mat.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
    }
  }

  dat <- data.frame(fread(paste(path.out, "/summary_MyCount_mat.txt",sep=""),header=F,stringsAsFactors=F,sep="\t", fill=T))
  write.table(t(dat),paste(path.out, "/summary_t.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  dat <- read.csv(paste(path.out, "/summary_t.txt",sep=""),header=T,stringsAsFactors=F,sep="\t")
  res <- c()
  for(i in 1:ncol(dat)){
    tmp.res <- dat[which(dat[,i]!=""),i]
    res <- c(res, tmp.res)
  }

  t.res <- table(res)

  for (i in case.list){
    ## (1) Create the aggregate reference set for this sample
    ## control check
    my.choice <- select.reference.set (test.counts =  MyCount.mat[,i], reference.counts = MyCount.mat[,control.list],
                                       bin.length = (MyCount.dafr$end - MyCount.dafr$start)/1000, n.bins.reduced = 10000)
    my.reference.selected <- apply(X = MyCount.mat[, my.choice$reference.choice, drop = FALSE], MAR = 1, FUN = sum)
    
    all.exons <- new('ExomeDepth', test = MyCount.mat[,i],
                     reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')
                     
    ################ Now call the CNVs
    all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4, chromosome = MyCount.dafr$chromosome,
                          start = MyCount.dafr$start, end = MyCount.dafr$end, name = MyCount.dafr$exon)
                          
    tmp.all.exons <- try(AnnotateExtra(x = all.exons, reference.annotation = target.genes.GRanges,
                                       min.overlap = 0.0001, column.name = 'regions'),silent=T)
                                       
    if(class(tmp.all.exons)!="try-error"){
      all.exons <- tmp.all.exons
    } 
    
    output.csv <- paste(path.out, "/", unlist(strsplit(i,split=".bam"))[1], ".csv", sep = '')
    
    if(nrow(all.exons@CNV.calls)>0){
      
      for(k in 1:nrow(all.exons@CNV.calls)){
        tmp.pos <- target.genes[which(target.genes[,1]==paste("chr",all.exons@CNV.calls[k,7],sep="")),]
        tmp.gene <- unique(tmp.pos[which(as.numeric(tmp.pos[,2])>=c(as.numeric(all.exons@CNV.calls[k,5]-1)) &
                                           as.numeric(tmp.pos[,3])<=as.numeric(all.exons@CNV.calls[k,6])),4])
        all.exons@CNV.calls[k,13] <- paste(tmp.gene,collapse=",")
      } 
      write.csv(file = output.csv, x = all.exons@CNV.calls, row.names = FALSE)
    } 
  } 
} 

exomedepth.f(path.case = "/home/DATA/YUHL/WES/Batch_macrogen1/BAM", path.out = "/home/DATA/YUHL/WES/Batch_macrogen1/CNV/exomedepth")


