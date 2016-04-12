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
##  Normalize by Adjusted GC content
##
## file = normalizeLoessGC.R
########
########## PREAMBLE ############

base <- ""
experimentName <- "ExpName"


setwd(base)
codePath <- "http://bioinfo-mite.crb.wsu.edu/Rcode"
rawPath <- file.path(getwd(),"Data")
figurePath <- file.path(getwd(), "Figures")
tablePath <- file.path(getwd(), "Tables")
pairPath <- file.path(getwd(),"Raw_Data_Files")
designPath <- file.path( base,"Design_Information")
probeAnnoPath <- file.path( base,"Probe_Anno")
library(limma)

targetsFile <- "TargetFile.txt"
designids <- c(designID)
ext <- ".pair"
################################################################################
### LOAD DATA FILES
#############################
### Load Raw Data ###
#############################

for (i in 1:length(designids)){
	design_id <- designids[i] 
	
	if (file.exists(file=file.path(rawPath,paste("UnmodifiedRawData",design_id,"RData",sep=".")))){
		load(file=file.path(rawPath,paste("UnmodifiedRawData",design_id,"RData",sep=".")))      
		#############################
		### Read ndf and Pos files
		#############################
		## load annotation information
		if(file.exists(file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))){
			load(file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
		} else {
			designName <- RG$targets$DESIGN_NAME[1]
			NDF <- read.table( file.path(designPath,paste(designName,"ndf",sep=".")) , 
					sep="\t",header=TRUE,as.is=TRUE)
			POS <- read.table( file.path(designPath,paste(designName,"pos",sep=".")) , 
					sep="\t",header=TRUE,as.is=TRUE)
			NDF <- NDF[match(RG$genes$PROBE_ID,NDF$PROBE_ID),]
			NDF$CHROMOSOME <- "RANDOM"
			NDF[match(POS$PROBE_ID,NDF$PROBE_ID),"CHROMOSOME"] <- POS$CHROMOSOME
			NDF$LENGTH <- nchar(NDF$PROBE_SEQUENCE)
			NDF$COUNT <- 0
			NDF[match(POS$PROBE_ID,NDF$PROBE_ID),"COUNT"] <- POS$COUNT
			library(Biostrings)
			NDF$GC <-sapply(gregexpr("[CG]",NDF$PROBE_SEQUENCE) ,length)
			NDF$CG <-sapply(gregexpr("CG",NDF$PROBE_SEQUENCE) ,length)
			NDF <- NDF[order(NDF$CHROMOSOME,NDF$POSITION),]
			POS <- NDF
			rm(NDF)
			save(POS,file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
		}
	} else {
		
		hybes  <- readTargets(targetsFile)
		hybes <- hybes[order(hybes$DESIGN_ID,hybes$ORD_ID,hybes$DYE),]
		hybes <- hybes[which(hybes$DESIGN_ID == design_id),]
		
		
		source(file.path(base,"readNimblegen.R"))
		design_name <- hybes[,"DESIGN_NAME"][1]
		
		narray <- length(unique(hybes[,"CHIP_ID"]))
		cy3FileName <- grep("[Cc][Yy]3", hybes$DYE)
		cy5FileName <- grep("[Cc][Yy]5",  hybes$DYE)
		files <- paste(hybes[c(cy3FileName,cy5FileName),"CHIP_ID",drop=TRUE],"_",rep(c(532,635),each=narray),ext,sep="")
		files <- matrix(files, ncol = 2, byrow = FALSE)
		RG <- readPairFiles(files, verbose=TRUE, path=pairPath)
		RG$Rb <- RG$Gb <- NULL 
		spottypes <- readSpotTypes("spottypes.txt")
		RG$genes$Status <- controlStatus(spottypes, RG$genes)
		RG$genes <- RG$genes[, match(c("PROBE_ID","Status"),colnames(RG$genes))]
		RG$targets <- cbind(hybes[cy3FileName,c("ORD_ID","CHIP_ID","DESIGN_NAME","DESIGN_ID","SAMPLE_SPECIES")],
				matrix(hybes$SAMPLE_DESCRIPTION[c(cy3FileName,cy5FileName)],ncol=2,byrow=FALSE,dimnames=list(NULL,c("Cy3","Cy5"))))
		if (!all(!duplicated(RG$genes$PROBE_ID)))
			RG <- RG[-which(duplicated(RG$genes$PROBE_ID)==TRUE),]
		
		
		#############################
		### Read ndf and Pos files
		#############################
		## load annotation information
		if(file.exists(file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))){
			load(file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
		} else {
			designName <- RG$targets$DESIGN_NAME[1]
			NDF <- read.table( file.path(designPath,paste(designName,"ndf",sep=".")) , 
					sep="\t",header=TRUE,as.is=TRUE)
			POS <- read.table( file.path(designPath,paste(designName,"pos",sep=".")) , 
					sep="\t",header=TRUE,as.is=TRUE)
			NDF <- NDF[match(RG$genes$PROBE_ID,NDF$PROBE_ID),]
			NDF$CHROMOSOME <- "RANDOM"
			NDF[match(POS$PROBE_ID,NDF$PROBE_ID),"CHROMOSOME"] <- POS$CHROMOSOME
			NDF$LENGTH <- nchar(NDF$PROBE_SEQUENCE)
			NDF$COUNT <- 0
			NDF[match(POS$PROBE_ID,NDF$PROBE_ID),"COUNT"] <- POS$COUNT
			library(Biostrings)
			NDF$GC <-sapply(gregexpr("[CG]",NDF$PROBE_SEQUENCE) ,length)
			NDF$CG <-sapply(gregexpr("CG",NDF$PROBE_SEQUENCE) ,length)
			NDF <- NDF[order(NDF$CHROMOSOME,NDF$POSITION),]
			POS <- NDF
			rm(NDF)
			save(POS,file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
		}
		RG <- RG[match(POS$PROBE_ID,RG$genes$PROBE_ID),]
		save(RG,file=file.path(rawPath,paste("UnmodifiedRawData",design_id,"RData",sep=".")))
		
		########### produce MA-plots and Histogram of RawData
		source(file.path(codePath,"maplot.R"))
		library(affy)
		narray <- ncol(RG$R)
		pdf(file=file.path(figurePath,paste(experimentName,design_id,"RawHist","pdf",sep=".")),width=7,height=10,pointsize=10)
		par(mfrow=c(2,1))
		plotDensity(cbind(log2(RG$R),log2(RG$G)),col=rep(c("red","green"),each=narray),
				lty=c(1:narray,1:narray),xlab="Log-intensity",ylab="Density",lwd=2,
				main="Channel bias in log-intensity",sub=paste("ChipSet",design_id))
		plotDensity(log2(RG$R)-log2(RG$G), lwd=2,
				xlab="Log-fold change",ylab="Density",
				main="Raw data in log-fold change, [log2(R) - log2(G)]",
				sub=paste("ChipSet",design_id))
		dev.off()
		
		for(i in 1:narray){
			png(filename = file.path(figurePath,paste(experimentName,design_id,"RawMA","array", RG$targets$CHIP_ID[i],"png",sep=".")),
					width = 7, height = 7, units = "in",res=300,
					pointsize = 14)
			maplot(RG,array=i,show.statistics=TRUE,ylim=c(-5,5),
					main=paste("Array",  RG$targets$CHIP_ID[i],"\nChipset",design_id)) ## can't use sub
			dev.off()
		}
		
	
		minE <-  floor(min(log2(cbind(RG$R,RG$G))))
		maxE <-  ceiling(max(log2(cbind(RG$R,RG$G))))
		#### Boxplots by column
		png(file=file.path(figurePath,paste(experimentName,design_id,"qcBP","png",sep=".")),
				height=10,width=7,pointsize=8,res=300,units="in")
		par(mfrow=c(3,1))
		boxplot(log2(RG$R) ~ col(RG$R),ylim=c(minE,maxE),col="red", 
				main=paste("all arrays Chipset",design_id,"Red"))
		boxplot(log2(RG$G) ~ col(RG$G),ylim=c(minE,maxE),col="green", 
				main=paste("all arrays Chipset",design_id,"Green"))
		boxplot(log2(RG$R) - log2(RG$G) ~ col(RG$R),col="blue", 
				main=paste("log-FC all arrays Chipset",design_id))
		dev.off()
		
		#### Boxplots by probelengths
		png(file=file.path(figurePath,paste(experimentName,design_id,"qcBPprobelengths","png",sep=".")),
				height=10,width=7,pointsize=8,res=300,units="in")
		par(mfrow=c(3,1))
		boxplot(log2(RG$R) ~ POS$LENGTH,ylim=c(minE,maxE),col="red", 
				main=paste("probe length all arrays Chipset",design_id,"Red"))
		boxplot(log2(RG$G) ~ POS$LENGTH,ylim=c(minE,maxE),col="green", 
				main=paste("probe length all arrays Chipset",design_id,"Green"))
		boxplot(log2(RG$R) - log2(RG$G) ~ POS$LENGTH,col="blue", 
				main=paste("probe length vs log-FC all arrays Chipset",design_id))
		dev.off()
		
		#### Boxplots by probe copy number
		png(file=file.path(figurePath,paste(experimentName,design_id,"qcBPprobecopy","png",sep=".")),
				height=10,width=7,pointsize=8,res=300,units="in")
		par(mfrow=c(3,1))
		boxplot(log2(RG$R) ~ POS$COUNT,ylim=c(minE,maxE),col="red", 
				main=paste("copy count all arrays Chipset",design_id,"Red"))
		boxplot(log2(RG$G) ~ POS$COUNT,ylim=c(minE,maxE),col="green", 
				main=paste("copy count all arrays Chipset",design_id,"Green"))
		boxplot(log2(RG$R) - log2(RG$G) ~ POS$COUNT,col="blue", 
				main=paste("probe count vs log-FC all arrays Chipset",design_id))
		dev.off()
		
		## Intensity by CpG
		png(file=file.path(figurePath,paste(experimentName,design_id,"qcBPCpG","png",sep=".")),
				height=10,width=7,pointsize=8,res=300,units="in")
		par(mfrow=c(3,1))
		boxplot(log2(RG$R) ~ POS$CG,ylim=c(minE,maxE),col="red",
				main=paste("Intensity vs GC count all arrays Chipset",design_id,"Red"))
		boxplot(log2(RG$G) ~ POS$CG,ylim=c(minE,maxE),col="green", 
				main=paste("Intensity vs GC arrays Chipset",design_id,"Green"))
		boxplot(log2(RG$R)-log2(RG$G) ~ POS$CG,col="blue", 
				main=paste("log-FC vs GC count arrays Chipset",design_id))
		dev.off()
		
		## Intensity by GC    
		png(file=file.path(figurePath,paste(experimentName,design_id,"qcBPGC","png",sep=".")),
				height=10,width=7,pointsize=8,res=300,units="in")
		par(mfrow=c(3,1))
		boxplot(log2(RG$R) ~ POS$GC,ylim=c(minE,maxE),col="red",
				main=paste("Intensity vs GC count all arrays Chipset",design_id,"Red"))
		boxplot(log2(RG$G) ~ POS$GC,ylim=c(minE,maxE),col="green", 
				main=paste("Intensity vs GC arrays Chipset",design_id,"Green"))
		boxplot(log2(RG$R)-log2(RG$G) ~ POS$GC,col="blue", 
				main=paste("log-FC vs GC count arrays Chipset",design_id))
		dev.off()
		
		corc1 <- sapply(names(table(POS$GC)),function(x) cor(as.numeric(log2(RG$R[which(POS$GC == x),])),as.numeric(log2(RG$G[which(POS$GC == x),]))))
		png(file=file.path(figurePath,paste(experimentName,design_id,"qcScatter","png",sep=".")),
				height=10,width=7,pointsize=8,res=300,units="in")
		par(mfrow=c(2,1))
		plot(log2(RG$R[which(POS$GC == 11),]),log2(RG$G[which(POS$GC == 11),]), 
				main=paste("Log Intensity for GC=11, Chipset",design_id),
				ylim=c(minE,maxE),xlim=c(minE,maxE),
				xlab="Cy5 log Intensity",ylab="Cy3 log Intensity")
		text(maxE-1,minE+1,paste("Correlation =",signif(corc1["11"],3)))
		plot(log2(RG$R[which(POS$GC == 39),]),log2(RG$G[which(POS$GC == 39),]), 
				main=paste("Log Intensity for GC=39, Chipset",design_id),
				ylim=c(minE,maxE),xlim=c(minE,maxE),
				xlab="Cy5 log Intensity",ylab="Cy3 log Intensity")
		text(maxE-1,minE+1,paste("Correlation =",signif(corc1["39"],3)))
		dev.off()
		## qc cor
		png(file=file.path(figurePath,paste(experimentName,design_id,"qcCorr","png",sep=".")),
				height=10,width=7,pointsize=8,res=300,units="in")
		plot(as.numeric(names(table(POS$GC))),corc1,type="l",main=paste("Correlation by GC count Chipset",design_id),xlab="GC count",ylab="Pearson correlation")
		dev.off()
		## qc hist
		png(file=file.path(figurePath,paste(experimentName,design_id,"qcHist","png",sep=".")),
				height=10,width=7,pointsize=8,res=300,units="in")
		par(mfrow=c(2,1))
		hist(log2(RG$G[,1]),breaks=100,
				main= paste("Green Channel, G-C bias in log-intensity (GC = 20), Chipset",design_id,"Array:",RG$targets$CHIP_ID[1]),
				xlim=c(minE,maxE),xlab="Log-intensity",ylab="Frequency")
		hist(log2(RG$G[which(POS$GC <=20),1]),
				col="green",breaks=100,add=TRUE)
		hist(log2(RG$R[,1]),breaks=100,
				main= paste("Red Channel, G-C bias in log-intensity (GC = 20), Chipset",design_id,"Array:",RG$targets$CHIP_ID[1]),
				xlim=c(minE,maxE),xlab="Log-intensity",ylab="Frequency")
		hist(log2(RG$R[which(POS$GC <=20),1]),col="red",breaks=100,add=TRUE)
		dev.off()
		
	}
	#######################################################
	## plot quantiles
	Rq <- as.data.frame(apply(log2(RG$R),2,quantile,probs = seq(0, 1, 0.001),na.rm=TRUE))
	Gq <- as.data.frame(apply(log2(RG$G),2,quantile,probs = seq(0, 1, 0.001),na.rm=TRUE))
	
	png(file=file.path(figurePath,paste(experimentName,design_id,"RvGquantiles","png",sep=".")),height=7,width=7,pointsize=8,res=300,units="in")
	matplot(Rq,Gq,type="p",col=1:ncol(RG$R),cex=0.5)
	abline(a=0,b=1)
	dev.off()
	
	### Remove Saturated Probes
	sat <- 2^(15.5)
	
	removeA1 <- c(apply(RG$R,2,function(x) which (x > sat)),
			apply(RG$G,2,function(x) which (x > sat)))
	matrix(lapply(removeA1,length),ncol=2,dimnames=list(RG$targets$CHIP_ID,c("Red","Green")))

	length(unique(unlist(removeA1)))
	duplicated <- which(POS$COUNT > 1)
	length(duplicated)

	RG$R[as.numeric((unique(c(unlist(removeA1),duplicated)))),] <- NA
	RG$G[as.numeric((unique(c(unlist(removeA1),duplicated)))),] <- NA
	## removed saturated and duplicated probes
	save(RG,file=file.path(rawPath,paste("PreprocRawData",design_id,"RData",sep=".")))
		
	
	################################################33
	# Log Values
	# MA2C
	logR <- log2(RG$R)
	logG <- log2(RG$G)
	## Adjust and GC content
	narray <- ncol(RG$R)

	min <- 2
	GCcount = table(as.factor(POS$GC))
	Rmean = apply(logR,2, function (x) tapply(x,POS$GC, mean,na.rm=TRUE))
	Gmean = apply(logG,2, function (x) tapply(x,POS$GC, mean,na.rm=TRUE))
	Rvar = apply(logR,2, function(x) tapply(x,POS$GC,var,na.rm=TRUE))
	Gvar = apply(logG,2, function(x) tapply(x,POS$GC,var,na.rm=TRUE))
	gccov = matrix(NA,nrow=length(GCcount),ncol=ncol(logR))
	for(i in 1:ncol(logR)) 
		gccov[,i] <- tapply(1:nrow(logR),POS$GC,
				function(x) cov(cbind(logR[x,i],logG[x,i]),use="pairwise.complete.obs")[1,2])
	
	GCsum <- list(GCcount=GCcount,Rmean=Rmean,Gmean=Gmean,Rvar=Rvar,Gvar=Gvar,gccov=gccov)
	
	png(file=file.path(figurePath,paste(experimentName,design_id,"meanIntvsGC","png",sep=".")),height=7,width=7,pointsize=8,res=300,units="in")
	matplot(names(GCcount),cbind(Rmean,Gmean),type = "l",col=rep(c("red","green"),each=narray),lty=rep(1:narray,times=2),main=paste("Mean expression within GC Chip",design_id))
	dev.off()
	
	## no background correction
	object <- MA.RG(RG, bc.method = "none", offset = 0)
	if (is.vector(object$M)) 
		object$M <- as.matrix(object$M)
	if (is.vector(object$A)) 
		object$A <- as.matrix(object$A)
	
	## GC correction 
	plotit=TRUE
	pdf(file=file.path(figurePath,paste(experimentName,design_id,"loesscurvesbyGC","pdf",
							sep=".")),height=7,width=7,pointsize=8)
	for (j in 1:narray) {
		if (plotit){
			plot(object$A[,j],object$M[,j],type="n",cex = 0.5,ylim=c(-5,5),xlab="A",ylab="M",main=paste("Loess Curves by GC, Chipset", design_id,", Array:",RG$targets$CHIP_ID[1]))
			abline(0, 0, col = "blue")
		}
		cat("array", j,"\tGC:")
		# look into each GC levels, may vary from 1 to 50
		for (GCp in levels(factor(POS$GC)) ) {
			cat(GCp,":")
			# pick up the probes with GC level N (N=1..50)
			spots <- which(POS$GC == as.numeric(GCp) )
			cat(length(spots),"..")
			# Remove probes with NA and see we have more than 1 probe
			if (length(na.exclude(object$M[spots,j])) > 1){
				y <- object$M[spots, j]
				x <- object$A[spots, j]
				#w <- weights[spots, j]
				lfit <- loessFit(y, x, span = 0.3, iterations = 4)
				# lfit has two list (fitted and residuals)
				# replace spot intensity M with residuals
				object$M[spots, j] <- lfit$residuals
				if(plotit){
					o <- order(x)
					# order A
					A <- x[o]
					# order M with same order
					M <- lfit$fitted[o]
					# get the probes with unique values (A)
					o <- which(!duplicated(x))
					lines(approx(A[o], M[o]), col = "red", lwd = 1, lty = 1)
				}
			} else {
				# if only one spot with the GC level, mark the probe as NA
				object$M[spots,j] <- NA
			}
		}
		cat("\n")
	}
	dev.off()
	
	
	object <- normalizeBetweenArrays(object) # A-quantile , quantile norm the A values
	
	
	#####################
	tmp <- RG.MA(object)
	## plot quantiles
	Rq <- as.data.frame(apply(log2(tmp$R),2,quantile,probs = seq(0, 1, 0.001),na.rm=TRUE))
	Gq <- as.data.frame(apply(log2(tmp$G),2,quantile,probs = seq(0, 1, 0.001),na.rm=TRUE))
	
	png(file=file.path(figurePath,paste(experimentName,design_id,"RvGpostGCquantiles","png",sep=".")),height=7,width=7,pointsize=8,res=300,units="in")
	matplot(Rq,Gq,type="p",col=1:ncol(tmp$R),cex=0.5)
	abline(a=0,b=1)
	dev.off()
	
	library(affy)
	narray <- ncol(tmp$R)
	pdf(file=file.path(figurePath,paste(experimentName,design_id,"postGCHist","pdf",sep=".")),width=7,height=10,pointsize=10)
	par(mfrow=c(2,1))
	plotDensity(cbind(log2(tmp$R),log2(tmp$G)),col=rep(c("red","green"),each=narray),
			lty=c(1:narray,1:narray),xlab="Log-intensity",ylab="Density",lwd=2,
			main="Channel bias in log-intensity",sub=paste("ChipSet",design_id))
	plotDensity(log2(tmp$R)-log2(tmp$G), lwd=2,
			xlab="Log-fold change",ylab="Density",
			main="Raw data in log-fold change, [log2(R) - log2(G)]",
			sub=paste("ChipSet",design_id))
	dev.off()
	#####################
	
	############## MA2C
	Rmean = apply(logR,2, function (x) tapply(x,POS$GC, mean,na.rm=TRUE))
	numer <- (logR - Rmean[as.character(POS$GC),]) - (logG - Gmean[as.character(POS$GC),])
	denom <- sqrt(Gvar + Rvar - 2*gccov)[as.character(POS$GC),]
	MA2C <- numer/denom  ## adjusted further for correlation
	MA2C <- scale(MA2C,center=FALSE)
	MA2CA <-  (logR + logG)/2

	save(object,MA2C,MA2CA,GCsum,file=file.path(rawPath,paste("GCNormalized",design_id,"RData",sep=".")))
	
	#maplotNormalized
	source(file.path(codePath,"maplot.R"))
	for(i in 1:narray){
		png(filename=file.path(figurePath,paste(experimentName,design_id,"GCloessMA","array",RG$targets$CHIP_ID[i],"png",sep=".")),
				width = 7, height = 7, units = "in",res=300,
				pointsize = 14)
		maplot(object,array=i,show.statistics=TRUE,ylim=c(-5,5),
				main=paste("Array ", RG$targets$CHIP_ID[i],"\nChipset",design_id)) ## can't use sub
		dev.off()
	}
	#maplotNormalized
	source(file.path(codePath,"maplot.R"))
	for(i in 1:narray){
		png(filename=file.path(figurePath,paste(experimentName,design_id,"MA2CMA","array",RG$targets$CHIP_ID[i],"png",sep=".")),
				width = 7, height = 7, units = "in",res=300,
				pointsize = 14)
		maplot(M=MA2C,object=MA2CA,array=i,show.statistics=TRUE,ylim=c(-5,5),
				main=paste("Array ", RG$targets$CHIP_ID[i],"\nChipset",design_id)) ## can't use sub
		dev.off()
	}
}
