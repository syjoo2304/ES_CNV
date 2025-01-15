#!/usr/bin/env Rscript

# SourceTarget.txt 
## ./hg19_uniqueome.coverage.base-space.25.1.Wig ./ucsc.hg19.fasta

./TargetPerla.pl SourceTarget.txt Agilent_SureSelect_V5_covered.bed Target_region_230418


excavator.f <- function(path.case, path.out){
  
  dir.create(path.out)
  path.cont = "/home/DATA/YUHL/WES/CTRL/FASTQ/FASTQ_link/BAM"
  setwd(path.cont) 
  bam.list.cont <- dir(pattern=".bam$") #
  bam.list.cont = bam.list.cont[1:10]  
  setwd(path.case) 
  bam.list.case <- dir(pattern=".bam$")
  
  setwd("/home/program/EXCAVATOR_Package_v2.2/EXCAVATOR")
  
  # case
  k=1
  for(i in 1:length(bam.list.case)){
    SampleId = unlist(strsplit(bam.list.case[i],split=".",fixed=T))[1]
    write.table(paste("Target_region_230418 hg19 ", path.case, "/",unlist(bam.list.case)[i],
                      " ", SampleId, " T", k, sep=""),
                paste(path.out, "/ReadInput_", length(bam.list.case), "_", length(bam.list.cont), ".txt", sep=""),
                col.names=F, row.names=F, quote=F, append=T)
    k = k + 1
  }
  
 
  # control
   k=1
   for(i in 1:length(bam.list.cont)){
   ContrId = unlist(strsplit(bam.list.cont[i],split=".",fixed=T))[1]
    write.table(paste("Target_region_230418 hg19 ", path.cont, "/",unlist(bam.list.cont)[i],
                      " ", ContrId, " C", k, sep=""),
                paste(path.out, "/ReadInput_", length(bam.list.case), "_", length(bam.list.cont), ".txt", sep=""),
                col.names=F, row.names=F, quote=F, append=T)
    k = k + 1
  }

  # Running ReadPerla.pl
  ## --mode somatic: case & control= 1:1 matching
  ## --mode pooling

  system(paste("./ReadPerla.pl ", path.out, "/ReadInput_", length(bam.list.case), "_", length(bam.list.cont), ".txt ", path.out, " --mode pooling",sep=""))
}


# patient_wes: 9 samples
excavator.f(path.case = "/home/DATA/YUHL/WES/Batch_macrogen1/BAM", path.out = "/home/DATA/YUHL/WES/Batch_macrogen1/CNV")


