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
## 	Compute Running Medians
##
## file = Median300Bp.R
########
########## PREAMBLE ############

base <- ""
experimentName <- "ExpName"

setwd(base)
rawPath <- file.path(getwd(),"Data")
figurePath <- file.path(getwd(), "Figures")
tablePath <- file.path(getwd(), "Tables")
pairPath <- file.path(getwd(),"Raw_Data_Files")
designPath <- file.path( base,"Design_Information")
probeAnnoPath <- file.path(base,"Probe_Anno")
library(limma)

designids <- c(designID)

compChr <- paste("chr",c(1:20,"1_random","2_random","3_random","4_random",
		"5_random","6_random","7_random","8_random",
		"9_random","11_random","12_random","15_random",
		"16_random","19_random","20_random","Un","X","X_random"),sep="") 

design <- c(1,1)
groupings <- c(1,2)


########## compute running medians
source("computeRunningMedians.R")

for (i in 1:length(designids)){
#	i <- 1
	design_id <- designids[i] 
	load(file=file.path(rawPath,paste("GCNormalized",design_id,"RData",sep=".")))
	load(file=file.path(probeAnnoPath,paste("ndfpos",design_id,"RData",sep=".")))
	POS <- POS[match(object$genes$PROBE_ID,POS$PROBE_ID),]
	ord <- order(POS$CHROMOSOME,POS$POSITION)
	POS <- POS[ord,]
	object<- object[ord,]
	
	chr <- intersect(unique(POS$CHROMOSOME),compChr)
	
	smooth300.gc <- computeRM(
			object = object ,
			probeAnno = POS[,c("PROBE_ID","CHROMOSOME","POSITION","LENGTH","COUNT")] ,
			design=design ,
			groupings = groupings ,
			allChr = chr ,
			winHalfSize = 300 ,
			min.probes = 4 ,
			quant=0.5 ,
			combineReplicates = FALSE , 
			checkUnique = FALSE ,
			uniqueCodes = c(1) ,
			verbose = TRUE 
	)
	save(smooth300.gc,file=file.path(rawPath,paste("GCNormSmoothed",design_id,"RData",sep="."))) 
	
	nas <- which(is.na(smooth300.gc$M[,1]))
	smooth300.gc <- smooth300.gc[-nas,]
	POS <- POS[match(smooth300.gc$genes$PROBE_ID,POS$PROBE_ID),]
	pmid <- round(POS$POSITION + (POS$LENGTH/2))
	data <- data.frame(
			seqid = POS$CHROMOSOME,
			start = pmid,
			end = pmid, 
			attributes = paste(POS$PROBE_ID,POS$SEQ_ID,sep=";"))
	
}
