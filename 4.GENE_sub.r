
final <- function(mode="",path=""){

	#hearing_loss_geneset
	HL <- read.table("/home/syjoo/REF/HL_KNOWN/HL-KNOWN-GENE_v201030.txt",sep='\t',stringsAsFactors=F,header=T)
	HL.gene <- HL$Name

	setwd(path)
	x <- paste('list.',mode,sep='')
	x <- dir(pattern=".txt")	
	data2 <- data.frame()
	for (i in 1:length(x)){
		Sample.ID <- unlist(strsplit(x[i],split= paste('_dedup_',mode,'.txt',sep='')))[1]
		test2 <- read.table(x[i],sep='\t',stringsAsFactors=F,header=T)
		test2 <- test2[which(test2$GENE != '.'),]
		test2 <- unique(test2)
		test2$"Sample" <- Sample.ID
		data2 <- rbind(data2, test2,stringsAsFactors=F)
	}

	data2 <- data2[order(data2$GENE),]
	write.table(data2, paste("combined.",mode,".txt",sep=''),quote=F,sep='\t',row.names=F)

	Gene2 <- c()
	for (i in 1:nrow(data2)){
	gene_temp <- unlist(strsplit(data2[i,5],split=';'))
	        if (gene_temp != '.'){
	                Gene2 <- c(Gene2,gene_temp)
	                Gene2 <- unique(Gene2)
	        }
	}

	CNV2 <- intersect(HL.gene,Gene2)

	dat <- data.frame()
	for (i in 1:length(CNV2)){
	        dat.temp <- data2[grepl(CNV2[i],data2$GENE),]
		dat.temp$"HLgene" <- CNV2[i]
		dat.temp <- dat.temp[,c("chr","start","end","type","HLgene","GENE","Sample")]
	        dat <- rbind(dat,dat.temp,stringsAsFactors=F)
	}

	write.table(dat, paste("HL.",mode,".txt",sep=''),quote=F,sep='\t', row.names=F)

}
