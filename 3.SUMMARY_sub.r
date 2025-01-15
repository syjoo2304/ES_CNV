over.f <- function(pd){
  crt.path <- getwd()
  setwd(paste(pd, "/exomedepth/", sep=""))
  res.list <- dir(pattern=".csv")
  id.list <- c()
  for(i in 1:length(res.list)){
    id.list <- c(id.list, unlist(strsplit(res.list[i],split=".csv"))[1])
  }
  setwd(pd)
  
  for(id in id.list){
    dat.ed <- read.table(paste(pd, "/exomedepth/", id, ".csv", sep=""),header=T,stringsAsFactors=F,sep=",")
    dat.ev <- read.table(paste(pd, "/excavator/Results/", id, "/FastCallResults_", id, ".txt", sep=""),header=T,stringsAsFactors=F)
    dat.ed$chromosome <- paste("chr",dat.ed$chromosome,sep="")
    dat.ev <- dat.ev[which(dat.ev$ProbCall >= 0.8),]
    all.chr <- paste("chr",c(1:22,"X","Y"),sep="")
    if(sum(dir()%in%"output_overlap")==0){
      dir.create("OVERLAP")
    }
    for(itr.chr in all.chr){
      if(sum(dat.ed$chromosome==itr.chr) * sum(dat.ev$Chromosome==itr.chr)>0){
        tmp.dat.ed <- data.frame(dat.ed[which(dat.ed$chromosome==itr.chr),])
        tmp.dat.ev <- data.frame(dat.ev[which(dat.ev$Chromosome==itr.chr),])
        for(itr in 1:nrow(tmp.dat.ed)){
          tmp.dat.ed.itr.start <- tmp.dat.ed[itr, 'start']
          tmp.dat.ed.itr.end   <- tmp.dat.ed[itr, 'end']
          tmp.dat.ed.itr.type  <- tmp.dat.ed[itr, 'type']
          case.1 <- which(tmp.dat.ev[,'Start'] == tmp.dat.ed.itr.start & tmp.dat.ev[,'End'] == tmp.dat.ed.itr.end)
          case.2 <- which(tmp.dat.ev[,'Start'] == tmp.dat.ed.itr.start & tmp.dat.ev[,'End'] < tmp.dat.ed.itr.end)
          case.3 <- which(tmp.dat.ev[,'Start'] == tmp.dat.ed.itr.start & tmp.dat.ev[,'End'] > tmp.dat.ed.itr.end)
          case.4 <- which(tmp.dat.ev[,'Start'] < tmp.dat.ed.itr.start & tmp.dat.ev[,'End'] == tmp.dat.ed.itr.end)
          case.5 <- which(tmp.dat.ev[,'Start'] < tmp.dat.ed.itr.start & tmp.dat.ev[,'End'] < tmp.dat.ed.itr.end & tmp.dat.ev[,'End'] >= tmp.dat.ed.itr.start)
          case.6 <- which(tmp.dat.ev[,'Start'] < tmp.dat.ed.itr.start & tmp.dat.ev[,'End'] > tmp.dat.ed.itr.end)
          case.7 <- which(tmp.dat.ev[,'Start'] > tmp.dat.ed.itr.start & tmp.dat.ev[,'End'] < tmp.dat.ed.itr.end)
          if(length(case.1)>0){
            tmp.det.ev.c1 <- tmp.dat.ev[case.1,c(1,2,3,8,6,7)]
            if(tmp.dat.ed.itr.type=="deletion"){
              if(length(which(tmp.det.ev.c1[,6]<0))>0){
                res.1 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c1[which(tmp.det.ev.c1[,6]<0),], tmp.dat.ed[itr, c(7,5,6)])
                write.table(res.1, paste(pd, "/OVERLAP/", id, "_case1.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(tmp.dat.ed.itr.type=="duplication"){
              if(length(which(tmp.det.ev.c1[,6]>0))>0){
                res.1 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c1[which(tmp.det.ev.c1[,6]>0),], tmp.dat.ed[itr, c(7,5,6)])
                write.table(res.1, paste(pd, "/OVERLAP/", id, "_case1.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
          if(length(case.2)>0){
            tmp.det.ev.c2 <- tmp.dat.ev[case.2,c(1,2,3,8,6,7)]
            if(tmp.dat.ed.itr.type=="deletion"){
              if(length(which(tmp.det.ev.c2[,6]<0))>0){
                res.2 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c2[which(tmp.det.ev.c2[,6]<0),], tmp.det.ev.c2[which(tmp.det.ev.c2[,6]<0),c(1,2,3)])
                write.table(res.2, paste(pd, "/OVERLAP/", id, "_case2.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(tmp.dat.ed.itr.type=="duplication"){
              if(length(which(tmp.det.ev.c2[,6]>0))>0){
                res.2 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c2[which(tmp.det.ev.c2[,6]>0),], tmp.det.ev.c2[which(tmp.det.ev.c2[,6]>0),c(1,2,3)])
                write.table(res.2, paste(pd, "/OVERLAP/", id, "_case2.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
          if(length(case.3)>0){
            tmp.det.ev.c3 <- tmp.dat.ev[case.3,c(1,2,3,8,6,7)]
            if(tmp.dat.ed.itr.type=="deletion"){
              if(length(which(tmp.det.ev.c3[,6]<0))>0){
                res.3 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c3[which(tmp.det.ev.c3[,6]<0),], tmp.dat.ed[itr, c(7,5,6)])
                write.table(res.3, paste(pd, "/OVERLAP/", id, "_case3.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(tmp.dat.ed.itr.type=="duplication"){
              if(length(which(tmp.det.ev.c3[,6]>0))>0){
                res.3 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c3[which(tmp.det.ev.c3[,6]>0),], tmp.dat.ed[itr, c(7,5,6)])
                write.table(res.3, paste(pd, "/OVERLAP/", id, "_case3.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
          if(length(case.4)>0){
            tmp.det.ev.c4 <- tmp.dat.ev[case.4,c(1,2,3,8,6,7)]
            if(tmp.dat.ed.itr.type=="deletion"){
              if(length(which(tmp.det.ev.c4[,6]<0))>0){
                res.4 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c4[which(tmp.det.ev.c4[,6]<0),], tmp.dat.ed[itr, c(7,5,6)])
                write.table(res.4, paste(pd, "/OVERLAP/", id, "_case4.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(tmp.dat.ed.itr.type=="duplication"){
              if(length(which(tmp.det.ev.c4[,6]>0))>0){
                res.4 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c4[which(tmp.det.ev.c4[,6]>0),], tmp.dat.ed[itr, c(7,5,6)])
                write.table(res.4, paste(pd, "/OVERLAP/", id, "_case4.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
          if(length(case.5)>0){
            tmp.det.ev.c5 <- tmp.dat.ev[case.5,c(1,2,3,8,6,7)]
            if(tmp.dat.ed.itr.type=="deletion"){
              if(length(which(tmp.det.ev.c5[,6]<0))>0){
                res.5 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c5[which(tmp.det.ev.c5[,6]<0),], tmp.dat.ed[itr, c(7,5)], tmp.det.ev.c5[which(tmp.det.ev.c5[,6]<0),c(3)])
                write.table(res.5, paste(pd, "/OVERLAP/", id, "_case5.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(tmp.dat.ed.itr.type=="duplication"){
              if(length(which(tmp.det.ev.c5[,6]>0))>0){
                res.5 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c5[which(tmp.det.ev.c5[,6]>0),], tmp.dat.ed[itr, c(7,5)], tmp.det.ev.c5[which(tmp.det.ev.c5[,6]>0),c(3)])
                write.table(res.5, paste(pd, "/OVERLAP/", id, "_case5.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
          if(length(case.6)>0){
            tmp.det.ev.c6 <- tmp.dat.ev[case.6,c(1,2,3,8,6,7)]
            if(tmp.dat.ed.itr.type=="deletion"){
              if(length(which(tmp.det.ev.c6[,6]<0))>0){
                res.6 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c6[which(tmp.det.ev.c6[,6]<0),], tmp.dat.ed[itr, c(7,5,6)])
                write.table(res.6, paste(pd, "/OVERLAP/", id, "_case6.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(tmp.dat.ed.itr.type=="duplication"){
              if(length(which(tmp.det.ev.c6[,6]>0))>0){
                res.6 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c6[which(tmp.det.ev.c6[,6]>0),], tmp.dat.ed[itr, c(7,5,6)])
                write.table(res.6, paste(pd, "/OVERLAP/", id, "_case6.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
          if(length(case.7)>0){
            tmp.det.ev.c7 <- tmp.dat.ev[case.7,c(1,2,3,8,6,7)]
            if(tmp.dat.ed.itr.type=="deletion"){
              if(length(which(tmp.det.ev.c7[,6]<0))>0){
                res.7 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c7[which(tmp.det.ev.c7[,6]<0),], tmp.det.ev.c7[which(tmp.det.ev.c7[,6]<0),c(1,2,3)])
                write.table(res.7, paste(pd, "/OVERLAP/", id, "_case7.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(tmp.dat.ed.itr.type=="duplication"){
              if(length(which(tmp.det.ev.c7[,6]>0))>0){
                res.7 <- cbind(tmp.dat.ed[itr,c(7,5,6,9,10,11,12,4,3)], tmp.det.ev.c7[which(tmp.det.ev.c7[,6]>0),], tmp.det.ev.c7[which(tmp.det.ev.c7[,6]>0),c(1,2,3)])
                write.table(res.7, paste(pd, "/OVERLAP/", id, "_case7.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
        }
      }
    }
  }
}

ind.f <- function(pd){
  crt.path <- getwd()
  setwd(paste(pd, "/exomedepth/", sep=""))
  res.list <- dir(pattern=".csv")
  id.list <- c()
  for(i in 1:length(res.list)){
    id.list <- c(id.list, unlist(strsplit(res.list[i],split=".csv"))[1])
  }
  setwd(pd)
  crt.dir <- getwd()
  if(sum(dir()%in%"SPLIT")==0){
    dir.create("SPLIT")
  } 
  setwd("SPLIT")
  if(sum(dir()%in%"del")==0){
    dir.create("del")
  }
  if(sum(dir()%in%"dup")==0){
    dir.create("dup")
  }
  setwd(paste(pd, "/OVERLAP/",sep=""))
  
  for(id in id.list){
    tmp.res <- dir(pattern=paste(id, "_case", sep=""))
    if(length(tmp.res)>0){
      for(i in 1:length(tmp.res)){
        tmp.output <- read.table(paste(pd, "/OVERLAP/", tmp.res[i], sep=""), header=F, stringsAsFactors=F)
        colnames(tmp.output)<- c("chromosome", "start", "end", "BF", "reads.expected", "reads.observed", 
                                 "reads.ratio", "nexons", "type", "Chromosome", "Start", "End", "ProbCall", 
                                 "CN", "Call", "F.CHR", "F.START", "F.END")
        if(i==1){
          final.output <- tmp.output
        }
        if(i!=1){
          final.output <- rbind(final.output, tmp.output)
        }
      }
      summary.output <- final.output[,c("F.CHR","F.START","F.END","Call")]
      if(length(which(summary.output[,4]<0))>0){
        write.table(summary.output[which(summary.output[,4]<0),], 
                    paste(pd, "/SPLIT/del/", id, "_deletion.txt", sep=""),
                    col.names=c("chr","start","end","type"),row.names=F,quote=F,sep="\t")
      }
      if(length(which(summary.output[,4]>0))>0){
        write.table(summary.output[which(summary.output[,4]>0),], 
                    paste(pd, "/SPLIT/dup/", id, "_duplication.txt", sep=""),
                    col.names=c("chr","start","end","type"),row.names=F,quote=F,sep="\t")
      }
    }
  }
}

comp.f = function(pd){
  crt.path <- getwd()
  setwd(paste(pd, "/exomedepth/", sep=""))
  res.list <- dir(pattern=".csv")
  id.list <- c()
  for(i in 1:length(res.list)){
    id.list <- c(id.list, unlist(strsplit(res.list[i],split=".csv"))[1])
  }
  setwd(pd)
  crt.dir <- getwd()
  if(sum(dir()%in%"OVERLAP_SAMPLE")==0){
    dir.create("OVERLAP_SAMPLE")
  }
  setwd("OVERLAP_SAMPLE")
  if(sum(dir()%in%"del")==0){
    dir.create("del")
  }
  if(sum(dir()%in%"dup")==0){
    dir.create("dup")
  }
  setwd(crt.dir)
  
  for(i in 1:(length(id.list)-1)){
    for(j in c(i+1):length(id.list)){
      id1=id.list[i]
      id2=id.list[j]
      
      # deletion
      dat.del.1 <-try(read.table(paste(pd, "/SPLIT/del/", id1, "_deletion.txt", sep=""), header=T, stringsAsFactors=F),silent=T)
      dat.del.2 <-try(read.table(paste(pd, "/SPLIT/del/", id2, "_deletion.txt", sep=""), header=T, stringsAsFactors=F),silent=T)
      
      if(class(dat.del.1)!="try-error" & class(dat.del.2)!="try-error"){
        for(itr in 1:nrow(dat.del.1)){
          if(length(which(dat.del.2[,1]==dat.del.1[itr,1]))>0){
            sub.dat.del.2 <- dat.del.2[which(dat.del.2[,1]==dat.del.1[itr,1]),]
            
            tmp.dat.del.1.chr <- dat.del.1[itr,1]
            tmp.dat.del.1.start <- dat.del.1[itr,2]
            tmp.dat.del.1.end   <- dat.del.1[itr,3]
            tmp.dat.del.1.type  <- dat.del.1[itr,4]
            
            case.1 <- which(sub.dat.del.2[,'start'] == tmp.dat.del.1.start & sub.dat.del.2[,'end'] == tmp.dat.del.1.end)
            case.2 <- which(sub.dat.del.2[,'start'] == tmp.dat.del.1.start & sub.dat.del.2[,'end'] < tmp.dat.del.1.end)
            case.3 <- which(sub.dat.del.2[,'start'] == tmp.dat.del.1.start & sub.dat.del.2[,'end'] > tmp.dat.del.1.end)
            case.4 <- which(sub.dat.del.2[,'start'] < tmp.dat.del.1.start & sub.dat.del.2[,'end'] == tmp.dat.del.1.end)
            case.5 <- which(sub.dat.del.2[,'start'] < tmp.dat.del.1.start & sub.dat.del.2[,'end'] < tmp.dat.del.1.end & sub.dat.del.2[,'end'] >= tmp.dat.del.1.start)
            case.6 <- which(sub.dat.del.2[,'start'] < tmp.dat.del.1.start & sub.dat.del.2[,'end'] > tmp.dat.del.1.end)
            case.7 <- which(sub.dat.del.2[,'start'] > tmp.dat.del.1.start & sub.dat.del.2[,'end'] < tmp.dat.del.1.end)
            
            if(length(case.1)>0){
              tmp.sub.dat.del.2.c1 <- sub.dat.del.2[case.1,]
              if(length(which(tmp.sub.dat.del.2.c1[,4]<0))>0){
                res.1 <- cbind(dat.del.1[itr,], tmp.sub.dat.del.2.c1[which(tmp.sub.dat.del.2.c1[,4]<0),], dat.del.1[itr, c(1,2,3)])
                write.table(res.1, paste(pd, "/OVERLAP_SAMPLE/del/", id1, "_", id2, "_case1.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.2)>0){
              tmp.sub.dat.del.2.c2 <- sub.dat.del.2[case.2,]
              if(length(which(tmp.sub.dat.del.2.c2[,4]<0))>0){
                res.2 <- cbind(dat.del.1[itr,], tmp.sub.dat.del.2.c2[which(tmp.sub.dat.del.2.c2[,4]<0),], tmp.sub.dat.del.2.c2[which(tmp.sub.dat.del.2.c2[,4]<0),c(1,2,3)])
                write.table(res.2, paste(pd, "/OVERLAP_SAMPLE/del/", id1, "_", id2, "_case2.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.3)>0){
              tmp.sub.dat.del.2.c3 <- sub.dat.del.2[case.3,]
              if(length(which(tmp.sub.dat.del.2.c3[,4]<0))>0){
                res.3 <- cbind(dat.del.1[itr,], tmp.sub.dat.del.2.c3[which(tmp.sub.dat.del.2.c3[,4]<0),], dat.del.1[itr, c(1,2,3)])
                write.table(res.3, paste(pd, "/OVERLAP_SAMPLE/del/", id1, "_", id2, "_case3.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.4)>0){
              tmp.sub.dat.del.2.c4 <- sub.dat.del.2[case.4,]
              if(length(which(tmp.sub.dat.del.2.c4[,4]<0))>0){
                res.4 <- cbind(dat.del.1[itr,], tmp.sub.dat.del.2.c4[which(tmp.sub.dat.del.2.c4[,4]<0),], dat.del.1[itr, c(1,2,3)])
                write.table(res.4, paste(pd, "/OVERLAP_SAMPLE/del/", id1, "_", id2, "_case4.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.5)>0){
              tmp.sub.dat.del.2.c5 <- sub.dat.del.2[case.5,]
              if(length(which(tmp.sub.dat.del.2.c5[,4]<0))>0){
                res.5 <- cbind(dat.del.1[itr,], tmp.sub.dat.del.2.c5[which(tmp.sub.dat.del.2.c5[,4]<0),], dat.del.1[itr, c(1,2)], tmp.sub.dat.del.2.c5[which(tmp.sub.dat.del.2.c5[,4]<0),c(3)])
                write.table(res.5, paste(pd, "/OVERLAP_SAMPLE/del/", id1, "_", id2, "_case5.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.6)>0){
              tmp.sub.dat.del.2.c6 <- sub.dat.del.2[case.6,]
              if(length(which(tmp.sub.dat.del.2.c6[,4]<0))>0){
                res.6 <- cbind(dat.del.1[itr,], tmp.sub.dat.del.2.c6[which(tmp.sub.dat.del.2.c6[,4]<0),], dat.del.1[itr, c(1,2,3)])
                write.table(res.6, paste(pd, "/OVERLAP_SAMPLE/del/", id1, "_", id2, "_case6.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.7)>0){
              tmp.sub.dat.del.2.c7 <- sub.dat.del.2[case.7,]
              if(length(which(tmp.sub.dat.del.2.c7[,4]<0))>0){
                res.7 <- cbind(dat.del.1[itr,], tmp.sub.dat.del.2.c7[which(tmp.sub.dat.del.2.c7[,4]<0),], tmp.sub.dat.del.2.c7[which(tmp.sub.dat.del.2.c7[,4]<0),c(1,2,3)])
                write.table(res.7, paste(pd, "/OVERLAP_SAMPLE/del/", id1, "_", id2, "_case7.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
        }
      }
      
      # duplication
      dat.dup.1 <- try(read.table(paste(pd, "/SPLIT/dup/", id1, "_duplication.txt", sep=""), header=T, stringsAsFactors=F),silent=T)
      dat.dup.2 <- try(read.table(paste(pd, "/SPLIT/dup/", id2, "_duplication.txt", sep=""), header=T, stringsAsFactors=F),silent=T)
      
      if(class(dat.dup.1)!="try-error" & class(dat.dup.2)!="try-error"){
        for(itr in 1:nrow(dat.dup.1)){
          if(length(which(dat.dup.2[,1]==dat.dup.1[itr,1]))>0){
            sub.dat.dup.2 <- dat.dup.2[which(dat.dup.2[,1]==dat.dup.1[itr,1]),]
            
            tmp.dat.dup.1.chr <- dat.dup.1[itr,1]
            tmp.dat.dup.1.start <- dat.dup.1[itr,2]
            tmp.dat.dup.1.end   <- dat.dup.1[itr,3]
            tmp.dat.dup.1.type  <- dat.dup.1[itr,4]
            
            case.dup.1 <- which(sub.dat.dup.2[,'start'] == tmp.dat.dup.1.start & sub.dat.dup.2[,'end'] == tmp.dat.dup.1.end)
            case.dup.2 <- which(sub.dat.dup.2[,'start'] == tmp.dat.dup.1.start & sub.dat.dup.2[,'end'] < tmp.dat.dup.1.end)
            case.dup.3 <- which(sub.dat.dup.2[,'start'] == tmp.dat.dup.1.start & sub.dat.dup.2[,'end'] > tmp.dat.dup.1.end)
            case.dup.4 <- which(sub.dat.dup.2[,'start'] < tmp.dat.dup.1.start & sub.dat.dup.2[,'end'] == tmp.dat.dup.1.end)
            case.dup.5 <- which(sub.dat.dup.2[,'start'] < tmp.dat.dup.1.start & sub.dat.dup.2[,'end'] < tmp.dat.dup.1.end & sub.dat.dup.2[,'end'] >= tmp.dat.dup.1.start)
            case.dup.6 <- which(sub.dat.dup.2[,'start'] < tmp.dat.dup.1.start & sub.dat.dup.2[,'end'] > tmp.dat.dup.1.end)
            case.dup.7 <- which(sub.dat.dup.2[,'start'] > tmp.dat.dup.1.start & sub.dat.dup.2[,'end'] < tmp.dat.dup.1.end)
            
            if(length(case.dup.1)>0){
              tmp.sub.dat.dup.2.c1 <- sub.dat.dup.2[case.dup.1,]
              if(length(which(tmp.sub.dat.dup.2.c1[,4]>0))>0){
                res.1 <- cbind(dat.dup.1[itr,], tmp.sub.dat.dup.2.c1[which(tmp.sub.dat.dup.2.c1[,4]>0),], dat.dup.1[itr, c(1,2,3)])
                write.table(res.1, paste(pd, "/OVERLAP_SAMPLE/dup/", id1, "_", id2, "_case1.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.dup.2)>0){
              tmp.sub.dat.dup.2.c2 <- sub.dat.dup.2[case.dup.2,]
              if(length(which(tmp.sub.dat.dup.2.c2[,4]>0))>0){
                res.2 <- cbind(dat.dup.1[itr,], tmp.sub.dat.dup.2.c2[which(tmp.sub.dat.dup.2.c2[,4]>0),], tmp.sub.dat.dup.2.c2[which(tmp.sub.dat.dup.2.c2[,4]>0),c(1,2,3)])
                write.table(res.2, paste(pd, "/OVERLAP_SAMPLE/dup/", id1, "_", id2, "_case2.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.dup.3)>0){
              tmp.sub.dat.dup.2.c3 <- sub.dat.dup.2[case.dup.3,]
              if(length(which(tmp.sub.dat.dup.2.c3[,4]>0))>0){
                res.3 <- cbind(dat.dup.1[itr,], tmp.sub.dat.dup.2.c3[which(tmp.sub.dat.dup.2.c3[,4]>0),], dat.dup.1[itr, c(1,2,3)])
                write.table(res.3, paste(pd, "/OVERLAP_SAMPLE/dup/", id1, "_", id2, "_case3.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.dup.4)>0){
              tmp.sub.dat.dup.2.c4 <- sub.dat.dup.2[case.dup.4,]
              if(length(which(tmp.sub.dat.dup.2.c4[,4]>0))>0){
                res.4 <- cbind(dat.dup.1[itr,], tmp.sub.dat.dup.2.c4[which(tmp.sub.dat.dup.2.c4[,4]>0),], dat.dup.1[itr, c(1,2,3)])
                write.table(res.4, paste(pd, "/OVERLAP_SAMPLE/dup/", id1, "_", id2, "_case4.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.dup.5)>0){
              tmp.sub.dat.dup.2.c5 <- sub.dat.dup.2[case.dup.5,]
              if(length(which(tmp.sub.dat.dup.2.c5[,4]>0))>0){
                res.5 <- cbind(dat.dup.1[itr,], tmp.sub.dat.dup.2.c5[which(tmp.sub.dat.dup.2.c5[,4]>0),], dat.dup.1[itr, c(1,2)], tmp.sub.dat.dup.2.c5[which(tmp.sub.dat.dup.2.c5[,4]>0),c(3)])
                write.table(res.5, paste(pd, "/OVERLAP_SAMPLE/dup/", id1, "_", id2, "_case5.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.dup.6)>0){
              tmp.sub.dat.dup.2.c6 <- sub.dat.dup.2[case.dup.6,]
              if(length(which(tmp.sub.dat.dup.2.c6[,4]>0))>0){
                res.6 <- cbind(dat.dup.1[itr,], tmp.sub.dat.dup.2.c6[which(tmp.sub.dat.dup.2.c6[,4]>0),], dat.dup.1[itr, c(1,2,3)])
                write.table(res.6, paste(pd, "/OVERLAP_SAMPLE/dup/", id1, "_", id2, "_case6.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
            if(length(case.dup.7)>0){
              tmp.sub.dat.dup.2.c7 <- sub.dat.dup.2[case.dup.7,]
              if(length(which(tmp.sub.dat.dup.2.c7[,4]>0))>0){
                res.7 <- cbind(dat.dup.1[itr,], tmp.sub.dat.dup.2.c7[which(tmp.sub.dat.dup.2.c7[,4]>0),], tmp.sub.dat.dup.2.c7[which(tmp.sub.dat.dup.2.c7[,4]>0),c(1,2,3)])
                write.table(res.7, paste(pd, "/OVERLAP_SAMPLE/dup/", id1, "_", id2, "_case7.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
              }
            }
          }
        }
      }
    }
  }
}

summary.f <- function(pd){
  ## deletion
  setwd(paste(pd, "/OVERLAP_SAMPLE/del",sep=""))
  res.list <- dir(pattern=".txt")
  setwd(pd)
  
  if(length(res.list)>0){
    for(i in 1:length(res.list)){
      tmp.res <- read.table(paste(pd, "/OVERLAP_SAMPLE/del/", res.list[i],sep=""),header=F,stringsAsFactors=F)
      tmp.ids <- unlist(strsplit(res.list[i],split="_case"))[1]
      id1 <- paste(unlist(strsplit(tmp.ids,split="[_]"))[c(1)],collapse="_")
      id2 <- paste(unlist(strsplit(tmp.ids,split="[_]"))[c(2)],collapse="_")
      tmp.res <- cbind(tmp.res, id1, id2)
      colnames(tmp.res)<- c("chr.1", "start.1", "end.1", "CNV.1", 
                            "chr.2", "start.2", "end.2", "CNV.2", 
                            "F.CHR", "F.START", "F.END","id1","id2")
      output <- tmp.res[,c(9:11,12,1:4,13,5:8)]
      if(i==1){
        write.table(output, paste(pd, "/OVERLAP_SAMPLE/summary_del.txt", sep=""),col.names=T,row.name=F,quote=F,sep="\t")
      }
      if(i!=1){
        write.table(output, paste(pd, "/OVERLAP_SAMPLE/summary_del.txt", sep=""),col.names=F,row.name=F,quote=F,sep="\t",append=T)
      }
    }
  }
  
  ## duplication
  setwd(paste(pd, "/OVERLAP_SAMPLE/dup",sep=""))
  res.list <- dir(pattern=".txt")
  setwd(pd)
  
  if(length(res.list)>0){
    for(i in 1:length(res.list)){
      tmp.res <- read.table(paste(pd, "/OVERLAP_SAMPLE/dup/", res.list[i],sep=""),header=F,stringsAsFactors=F)
      tmp.ids <- unlist(strsplit(res.list[i],split="_case"))[1]
      id1 <- paste(unlist(strsplit(tmp.ids,split="[_]"))[c(1)],collapse="_")
      id2 <- paste(unlist(strsplit(tmp.ids,split="[_]"))[c(2)],collapse="_")
      tmp.res <- cbind(tmp.res, id1, id2)
      colnames(tmp.res)<- c("chr.1", "start.1", "end.1", "CNV.1", 
                            "chr.2", "start.2", "end.2", "CNV.2", 
                            "F.CHR", "F.START", "F.END","id1","id2")
      output <- tmp.res[,c(9:11,12,1:4,13,5:8)]
      if(i==1){
        write.table(output, paste(pd, "/OVERLAP_SAMPLE/summary_dup.txt",sep=""),col.names=T,row.name=F,quote=F,sep="\t")
      }
      if(i!=1){
        write.table(output, paste(pd, "/OVERLAP_SAMPLE/summary_dup.txt",sep=""),col.names=F,row.name=F,quote=F,sep="\t",append=T)
      }
    }
  }
}

final.f <- function(pd){
  for(in.dat in c("summary_del.txt", "summary_dup.txt")){
    out.dat = paste("final_", in.dat, sep="")
    dat <- try(read.table(paste(pd, "/OVERLAP_SAMPLE/", in.dat, sep=""), header=T, stringsAsFactors=F),silent=T)
    
    if(class(dat)!="try-error"){
      write.table(t(c("CHR","START","END","GENE","UCSC","IDs","Info")), 
                  paste(pd, "/OVERLAP_SAMPLE/", out.dat, sep=""), col.names=F, row.names=F, quote=F, sep="\t")
      bed <- read.table("/home/DATA/REF/Agilent_SureSelect_V5_covered.bed", header=F, stringsAsFactors=F)
      
      # ordering data set
      all.chr <- paste("chr",c(1:22,"X","Y"),sep="")
      itr=1
      for(i in all.chr){
        if(length(which(dat[,1]==i))>0){
          tmp.dat <- dat[which(dat[,1]==i),]
          o.tmp.dat <- tmp.dat[order(tmp.dat[,2],tmp.dat[,3],tmp.dat[,4]),]
          if(itr==1){
            res.dat <- o.tmp.dat	
          }
          else{
            res.dat <- rbind(res.dat, o.tmp.dat)
          }
          itr = itr + 1
        }
      }
      
      # checking overlap region
      for(i in all.chr){
        if(length(which(res.dat[,1]==i))>0){
          tmp.dat <- res.dat[which(res.dat[,1]==i),]
          itr.str <- unique(tmp.dat[,2])
          for(j in 1:length(itr.str)){
            tmp.dat.j <- tmp.dat[which(tmp.dat[,2]==itr.str[j]),]
            itr.end <- unique(tmp.dat.j[,3])
            for(k in 1:length(itr.end)){
              tmp.dat.j.k <- tmp.dat.j[which(tmp.dat.j[,3]==itr.end[k]),]
              for(l in 1:nrow(tmp.dat.j.k)){
                if(l==1){
                  tmp.dat.j.k.l.chr <- tmp.dat.j.k[l,1]
                  tmp.dat.j.k.l.str <- tmp.dat.j.k[l,2]
                  tmp.dat.j.k.l.end <- tmp.dat.j.k[l,3]
                  tmp.dat.j.k.l.id <- c(tmp.dat.j.k[l,4], tmp.dat.j.k[l,9])
                  tmp.dat.j.k.l.info <- paste(tmp.dat.j.k[l,4],":",tmp.dat.j.k[l,5],":",tmp.dat.j.k[l,6],"-",
                                              tmp.dat.j.k[l,7],":",tmp.dat.j.k[l,8],sep="")
                  tmp.dat.j.k.l.info <- c(tmp.dat.j.k.l.info, paste(tmp.dat.j.k[l,9],":",tmp.dat.j.k[l,10],":",
                                                                    tmp.dat.j.k[l,11],"-",tmp.dat.j.k[l,12],":",tmp.dat.j.k[l,13],sep=""))
                }
                tmp.dat.j.k.l.id <- unique(c(tmp.dat.j.k.l.id, tmp.dat.j.k[l,4], tmp.dat.j.k[l,9]))
                tmp.dat.j.k.l.info <- unique(c(tmp.dat.j.k.l.info, paste(tmp.dat.j.k[l,4],":",tmp.dat.j.k[l,5],":",
                                                                         tmp.dat.j.k[l,6],"-",tmp.dat.j.k[l,7],":",tmp.dat.j.k[l,8],sep=""),
                                               paste(tmp.dat.j.k[l,9],":",tmp.dat.j.k[l,10],":",tmp.dat.j.k[l,11],
                                                     "-",tmp.dat.j.k[l,12],":",tmp.dat.j.k[l,13],sep="")))
              }
              tmp.dat.j.k.l.ucsc <- paste("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=",
                                          tmp.dat.j.k.l.chr,":",tmp.dat.j.k.l.str,"-",tmp.dat.j.k.l.end,sep="")
              
              tmp.bed <- bed[which(bed[,1]==tmp.dat.j.k.l.chr),]
              
              case.1 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] == tmp.dat.j.k.l.end)
              case.2 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end)
              case.3 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] > tmp.dat.j.k.l.end)
              case.4 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] == tmp.dat.j.k.l.end)
              case.5 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end & tmp.bed[,3] >= tmp.dat.j.k.l.str)
              case.6 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] > tmp.dat.j.k.l.end)
              case.7 <- which(tmp.bed[,2] > tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end)
              case.all <- unique(c(case.1, case.2, case.3, case.4, case.5, case.6, case.7))
              
              if(length(case.all)==0){
                tmp.dat.j.k.gene <- "."
              }
              if(length(case.all)!=0){
                tmp.bed.all <- tmp.bed[case.all,]
                tmp.dat.j.k.gene <- c()
                for(m in 1:nrow(tmp.bed.all)){
                  tmp.bed.all.m <- unlist(strsplit(unlist(strsplit(tmp.bed.all[m,4],split="[,]"))[1],split="[|]"))
                  if(tmp.bed.all.m[1] == "ref"){
                    tmp.dat.j.k.m.gene <- tmp.bed.all.m[2]
                    tmp.dat.j.k.gene <- c(tmp.dat.j.k.gene, tmp.dat.j.k.m.gene)
                  }
                }
              }
              if(is.null(tmp.dat.j.k.gene)){
                tmp.dat.j.k.gene <- "."
              }
              write.table(t(c(tmp.dat.j.k.l.chr, tmp.dat.j.k.l.str, tmp.dat.j.k.l.end, paste(tmp.dat.j.k.gene,collapse=";"), 
                              tmp.dat.j.k.l.ucsc, paste(tmp.dat.j.k.l.id,collapse=";"), paste(tmp.dat.j.k.l.info,collapse=";"))), 
                          paste(pd, "/OVERLAP_SAMPLE/", out.dat, sep=""), col.names=F, row.names=F, quote=F, sep="\t", append=T)
            }
          }
        }
      }
    }
  }
}


gene_anno.f <- function(pd="",ref="",del.dir="",dup.dir=""){

#ref = "/home/DATA/REF"
dir.create(del.dir)
dir.create(dup.dir)
setwd(pd)

# id check
setwd("./SPLIT/del")
res.list <- dir(pattern=".txt")

id.list <- c()
for(i in 1:length(res.list)){
	id.list <- c(id.list, unlist(strsplit(res.list[i],split="_deletion.txt"))[1])
}

setwd(ref)

bed <- read.table("Agilent_SureSelect_V5_covered.bed", header=F, stringsAsFactors=F)

setwd(pd)

for(i in 1:length(id.list)){
	id=id.list[i]

# deletion
	dat.del <- try(read.table(paste(pd, "/SPLIT/del/", id, "_deletion.txt", sep=""), header=T, stringsAsFactors=F),silent=T)

	if(class(dat.del)!="try-error"){
	# ordering data set
		all.chr <- paste("chr",c(1:22,"X","Y"),sep="")
		itr=1
		for(j in all.chr){
			if(length(which(dat.del[,1]==j))>0){
				tmp.dat <- dat.del[which(dat.del[,1]==j),]
				o.tmp.dat <- tmp.dat[order(tmp.dat[,2],tmp.dat[,3]),]
				if(itr==1){
					res.dat <- o.tmp.dat
				}
				else{
					res.dat <- rbind(res.dat, o.tmp.dat)
				}
				itr = itr + 1
			}
		}

	# annotation
		for(j in 1:nrow(res.dat)){
			tmp.dat.j.k.l.chr <- res.dat[j,1]
			tmp.dat.j.k.l.str <- res.dat[j,2]
			tmp.dat.j.k.l.end <- res.dat[j,3]
			tmp.dat.j.k.l.type <- res.dat[j,4]
			tmp.bed <- bed[which(bed[,1]==tmp.dat.j.k.l.chr),]

			case.1 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] == tmp.dat.j.k.l.end)
			case.2 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end)
			case.3 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] > tmp.dat.j.k.l.end)
			case.4 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] == tmp.dat.j.k.l.end)
			case.5 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end & tmp.bed[,3] >= tmp.dat.j.k.l.str)
			case.6 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] > tmp.dat.j.k.l.end)
			case.7 <- which(tmp.bed[,2] > tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end)
			case.all <- unique(c(case.1, case.2, case.3, case.4, case.5, case.6, case.7))

			if(length(case.all)==0){
				tmp.dat.j.k.gene <- "."
			}

			if(length(case.all)!=0){
				tmp.bed.all <- tmp.bed[case.all,]
				tmp.dat.j.k.gene <- c()
				for(m in 1:nrow(tmp.bed.all)){
					tmp.bed.all.m <- unlist(strsplit(unlist(strsplit(tmp.bed.all[m,4],split="[,]"))[1],split="[|]"))
					if(tmp.bed.all.m[1] == "ref"){
						tmp.dat.j.k.m.gene <- tmp.bed.all.m[2]
						tmp.dat.j.k.gene <- c(tmp.dat.j.k.gene, tmp.dat.j.k.m.gene)
						tmp.dat.j.k.gene <- unique(tmp.dat.j.k.gene)
					}
				}
			}

			if(is.null(tmp.dat.j.k.gene)){
				tmp.dat.j.k.gene <- "."
			}

			if(j==1){
				write.table(t(c(tmp.dat.j.k.l.chr, tmp.dat.j.k.l.str, tmp.dat.j.k.l.end, tmp.dat.j.k.l.type, 
					paste(tmp.dat.j.k.gene,collapse=";"))), paste(pd, "/SPLIT/del_gene/", id, "_deletion.txt", sep=""), 
					col.names=c("chr","start","end","type","GENE"), row.names=F, quote=F, sep="\t")
			}
			if(j>1){
				write.table(t(c(tmp.dat.j.k.l.chr, tmp.dat.j.k.l.str, tmp.dat.j.k.l.end, tmp.dat.j.k.l.type, 
				paste(tmp.dat.j.k.gene,collapse=";"))), paste(pd, "/SPLIT/del_gene/", id, "_deletion.txt", sep=""), 
				col.names=F, row.names=F, quote=F, sep="\t", append=T)
			}
		}
	}

# duplication
	dat.dup <- try(read.table(paste(pd, "/SPLIT/dup/", id, "_duplication.txt", sep=""), header=T, stringsAsFactors=F),silent=T)


	if(class(dat.dup)!="try-error"){
	# ordering data set
		all.chr <- paste("chr",c(1:22,"X","Y"),sep="")
		itr=1
		for(j in all.chr){
			if(length(which(dat.dup[,1]==j))>0){
				tmp.dat <- dat.dup[which(dat.dup[,1]==j),]
				o.tmp.dat <- tmp.dat[order(tmp.dat[,2],tmp.dat[,3]),]
				if(itr==1){
					res.dat <- o.tmp.dat
				}
				else{
					res.dat <- rbind(res.dat, o.tmp.dat)
				}
				itr = itr + 1
			}
		}

	# annotation
		for(j in 1:nrow(res.dat)){
			tmp.dat.j.k.l.chr <- res.dat[j,1]
			tmp.dat.j.k.l.str <- res.dat[j,2]
			tmp.dat.j.k.l.end <- res.dat[j,3]
			tmp.dat.j.k.l.type <- res.dat[j,4]
			tmp.bed <- bed[which(bed[,1]==tmp.dat.j.k.l.chr),]

			case.1 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] == tmp.dat.j.k.l.end)
			case.2 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end)
			case.3 <- which(tmp.bed[,2] == tmp.dat.j.k.l.str & tmp.bed[,3] > tmp.dat.j.k.l.end)
			case.4 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] == tmp.dat.j.k.l.end)
			case.5 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end & tmp.bed[,3] >= tmp.dat.j.k.l.str)
			case.6 <- which(tmp.bed[,2] < tmp.dat.j.k.l.str & tmp.bed[,3] > tmp.dat.j.k.l.end)
			case.7 <- which(tmp.bed[,2] > tmp.dat.j.k.l.str & tmp.bed[,3] < tmp.dat.j.k.l.end)
			case.all <- unique(c(case.1, case.2, case.3, case.4, case.5, case.6, case.7))

			if(length(case.all)==0){
				tmp.dat.j.k.gene <- "."
			}

			if(length(case.all)!=0){
				tmp.bed.all <- tmp.bed[case.all,]
				tmp.dat.j.k.gene <- c()
				for(m in 1:nrow(tmp.bed.all)){
					tmp.bed.all.m <- unlist(strsplit(unlist(strsplit(tmp.bed.all[m,4],split="[,]"))[1],split="[|]"))
					if(tmp.bed.all.m[1] == "ref"){
						tmp.dat.j.k.m.gene <- tmp.bed.all.m[2]
						tmp.dat.j.k.gene <- c(tmp.dat.j.k.gene, tmp.dat.j.k.m.gene)
						tmp.dat.j.k.gene <- unique(tmp.dat.j.k.gene)
					}
				}
			}

			if(is.null(tmp.dat.j.k.gene)){
				tmp.dat.j.k.gene <- "."
			}

			if(j==1){
				write.table(t(c(tmp.dat.j.k.l.chr, tmp.dat.j.k.l.str, tmp.dat.j.k.l.end, tmp.dat.j.k.l.type, 
					paste(tmp.dat.j.k.gene,collapse=";"))), paste(pd, "/SPLIT/dup_gene/", id, "_duplication.txt", sep=""), 
					col.names=c("chr","start","end","type","GENE"), row.names=F, quote=F, sep="\t")
			}
			if(j>1){
				write.table(t(c(tmp.dat.j.k.l.chr, tmp.dat.j.k.l.str, tmp.dat.j.k.l.end, tmp.dat.j.k.l.type, 
				paste(tmp.dat.j.k.gene,collapse=";"))), paste(pd, "/SPLIT/dup_gene/", id, "_duplication.txt", sep=""), 
				col.names=F, row.names=F, quote=F, sep="\t", append=T)
			}
		}
	}
}
}
