####################################################
##
## Primary Researcher: Ben and Carlos
##		Washington State University
## 		Skinner labratory
##  Nimblegen Promotor array CHiP data Dataset
##  January 25 2008
##
## Bioinformatician: Matt Settles
##	Washington State University
##
## file = ComputeDiff.R
##  Modified by Md Haque
########## PREAMBLE ############

base <- ""

setwd(base)
rawPath <- file.path(getwd(),"Data")
figurePath <- file.path(getwd(), "Figures")
tablePath <- file.path(getwd(), "Tables")
pairPath <- file.path(getwd(),"Raw_Data_Files")
designPath <- file.path( base,"Design_Information")
probeAnnoPath <- file.path( base,"Probe_Anno")
library(limma)
library(Biobase)

designids <- c(designID)

allChr <- paste("chr",c(1:20,"1_random","2_random","3_random","4_random",
				"5_random","6_random","7_random","8_random",
				"9_random","11_random","12_random","15_random",
				"16_random","19_random","20_random","Un","X","X_random"),sep="") 

maxp <- 1e-5

for (i in 1:length(designids)){
	library(limma)
	design_id <- designids[i] 
	load(file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
	load(file=file.path(rawPath,paste("GCNormSmoothed",design_id,"RData",sep=".")))
	
	rowM <- rowMedians(smooth300.gc$M)
	rowA <- rowMedians(smooth300.gc$A)
	sdA <- apply(smooth300.gc$A,1,sd,na.rm=TRUE)
	meanM <- mean(rowM,na.rm=TRUE);sdM <- sd(rowM,na.rm=TRUE)
	rawp <- 2*pnorm(abs(rowM-meanM),mean=0,sd=sdM,lower.tail=FALSE)

	res <- data.frame(M=rowM,A=rowA,sdA = sdA,t=NA,rawp=rawp,row.names=rownames(smooth300.gc))
	write.table(res,file=file.path(tablePath,paste("tresults",design_id,"csv",sep=".")),sep=",",row.names=TRUE,col.names=TRUE)

	probes <- res[which(res$rawp < maxp),]	
	probes <- data.frame(probes,POS[match(rownames(probes),POS$PROBE_ID),c("CHROMOSOME","POSITION")])
	probes$POSITION <- probes$POSITION + ceiling(POS[match(rownames(probes),POS$PROBE_ID),c("LENGTH")]/2)
	
	probes$Dir <-  cut(probes$M, c(-10000,0 ,10000), labels = c("Cont","Vinc"))
	probes$rank= floor(-log10(probes$rawp))
	probes$cluster <- 1
	probes <- probes[order(probes$CHROMOSOME,as.numeric(as.character(probes$POSITION))),]
	for (i in 2:nrow(probes)){
		if ( probes[i,"CHROMOSOME"] == probes[i-1,"CHROMOSOME"]){
			if( (as.numeric(as.character(probes[i,"POSITION"])) - 300) < (as.numeric(as.character(probes[i-1,"POSITION"])) + 300) & (probes[i,"Dir"] == probes[i-1,"Dir"])){
				probes[i,"cluster"] <- probes[i-1,"cluster"]
			} else {
				probes[i,"cluster"] <- probes[i-1,"cluster"] + 1
			}
		} else {
			probes[i,"cluster"] <- probes[(i-1),"cluster"] + 1
		}       
	}
	clusterID <- tapply(1:length(probes$cluster),probes$cluster,function(x) paste(probes[["CHROMOSOME"]][x][1],":",min(probes[["POSITION"]][x])-300,"-",
	max(probes[["POSITION"]][x]+300),sep=""))
	probes$clusterloc <- clusterID[probes$cluster]
	
	#################### Clusters #####################################
	library(BSgenome)
	library("BSgenome.Rnorvegicus.UCSC.rn4")
	
	sregions <- split(probes,probes$cluster)
	regions <- sapply(sregions, function(x) {
				z <- c(ClusterID=x[1,"clusterloc"],Chromosome = x[1,"CHROMOSOME"],cSTART = min(x$POSITION)-300, cSTOP=max(x$POSITION)+300, meanM = mean(x$M),meanA = mean(x$A),minP = min(x$rawp),nprobes = nrow(x))
				a <- c(Score = max(x$rank))
				t(c(z,a))
			})
	regions <- as.data.frame(t(regions))
	colnames(regions) <- c("ClusterID","Chromosome","cSTART","cSTOP","meanM","meanA","minP","nProbes","Score")
	Seq <- apply(regions[,c("Chromosome","cSTART","cSTOP")],1,function(x) as.character(getSeq(Rnorvegicus, x[1], start=as.integer(x[2]),
	 end=as.integer(x[3]),as.character=FALSE)))

	base <- alphabetFrequency(DNAStringSet(Seq),baseOnly=TRUE)
	regions$Length <- nchar(Seq)
	regions$GC <- rowSums(base[,2:3])
	regions$CpG <- dinucleotideFrequency(DNAStringSet(Seq))[,"CG"]
	regions$CpGdensity <- (as.numeric(as.character(regions$CpG)) / ( as.numeric(as.character(regions$cSTOP)) - as.numeric(as.character(regions$cSTART)) ))
	regions$Seq <- Seq
	regions <- regions[which(regions$CpGdensity >= 0.01),]
	regions <- regions[order(regions$meanA),]
	regions <- regions[which(as.numeric(as.character(regions$meanA)) > 9.5),]
	write.table(probes,file=file.path(tablePath,paste("ProbeResults",design_id,"csv",sep=".")),sep=",",row.names=TRUE,col.names=TRUE)
	write.table(regions,file=file.path(tablePath,paste("RegionsResults",design_id,"csv",sep=".")),sep=",",row.names=TRUE,col.names=TRUE,quote=FALSE)

}

for (i in 1:length(designids)){
	design_id <- designids[i] 
	
	tmp <- read.table(file=file.path(tablePath,paste("RegionsResults",design_id,"csv",sep=".")),sep=",")
	if (i == 1)
		cresults <- tmp
	else
		cresults <- rbind(cresults,tmp)
}
write.table(cresults,file=file.path(tablePath,paste("RegionsCombined","csv",sep=".")),sep=",",row.names=TRUE,col.names=TRUE,quote=FALSE)


