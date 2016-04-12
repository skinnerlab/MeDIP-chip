library(limma)
### Read in 2-channel nimblegen array data
readNblgn <- function(targetsFile, spotTypesFile, pairpath=getwd(), verbose=TRUE, ...)
{
  # 0. check arguments:
  stopifnot(file.exists(targetsFile), file.exists(spotTypesFile))

  # 1. read in raw intensities:
  if (verbose) cat("Reading targets file...\n")
  hybes  <- readTargets(targetsFile, ...)
  cy3FileColumn <- grep("^[Ff]ile.*[Cc][Yy]3", colnames(hybes))
  cy5FileColumn <- grep("^[Ff]ile.*[Cc][Yy]5", colnames(hybes))
  if (length(cy3FileColumn)==0 | length(cy5FileColumn)==0){
    cat("Unexpected File Header of targets file:\n")
    print(hybes)
    #browser()
    while (!all(c(cy3FileColumn, cy5FileColumn) %in% 1:ncol(hybes))) {
      cy3FileColumn <- as.numeric(readline("Which column holds the Cy3 file names? "))
      cy5FileColumn <- as.numeric(readline("Which column holds the Cy5 file names? "))
    }# while
  } # if (length(cy3FileColum)==0 | length(cy5FileColumn)==0)  
  if (verbose) cat("Reading raw intensities...\n")
  RG <- readPairFiles(hybes[,c(cy3FileColumn,cy5FileColumn),drop=FALSE], verbose=verbose, path=pairpath,...)
  if (verbose) cat("Determining probe categories...\n")
  spottypes <- readSpotTypes(spotTypesFile)
  RG$genes$Status <- controlStatus(spottypes, RG$genes)
  if (is.null(RG$genes$ID)){
    RG$genes$ID <- as.character(RG$genes$PROBE_ID)}
  RG$targets <- hybes
  return(RG)
}# readNimblegen
  

readPairFiles <- function(files, path=NULL, ext=NULL, names=NULL, columns=NULL, wt.fun=NULL, verbose=TRUE, sep="\t", quote="\"",...)
{
    if (is.null(dim(files))) {
        if (length(files)%%2 == 0)
            files <- matrix(files, ncol = 2, byrow = TRUE)
        else stop("Odd number of files: should be two data files for each array" )
    }
    else {
        files <- as.matrix(files)
        if (ncol(files) != 2)
            stop("Need a two column matrix of file names")
    }
    if (!is.null(ext))
        files <- array(paste(files, ext, sep = "."), dim(files))
    narrays <- nrow(files)
    if (is.null(columns))
        columns <- list(f = "PM", b = "MM")
    if (is.null(columns$f) || is.null(columns$b))
        stop("'columns' should have components 'f' and 'b'")
    fullname <- files[1, 1]
    if (!is.null(path))
        fullname <- file.path(path, fullname)
    headers <- readNblgnHeader(fullname)
    if (verbose)
        cat("Read header information\n")
    skip <- headers$NHeaderRecords
    obj <- read.table(fullname, skip = skip, header = TRUE,
            sep = sep, quote = quote, as.is=TRUE, check.names = FALSE, comment.char = "",
            fill = TRUE, ...)
    nspots <- nrow(obj)
    YR <- YG  <- matrix(0, nspots, narrays)
    if (is.null(names)) {
        colnames(YG) <- removeExt(files[, 1])
        colnames(YR) <- removeExt(files[, 2])
        }
    else {
        colnames(YR) <- colnames(YG) <- names
        }
    RG <- list(R = YR, G = YG, Rb = YR, Gb = YG)
    if (!is.null(wt.fun))
        RG$weights <- YR
    for (i in 1:narrays) {
        fullname <- files[i, 1]
      if (!is.null(path))
            fullname <- file.path(path, fullname)
        if (i > 1) {
            headers <- readNblgnHeader(fullname)
        }
        obj <- read.table(fullname, skip = skip, header = TRUE,
            sep = sep, quote = quote, check.names = FALSE, comment.char = "",
            fill = TRUE, nrows = nspots, ...)
        if (verbose)
            cat(paste("Read", fullname, "\n"))
        if (i == 1)
            RG$genes <- obj[, c("GENE_EXPR_OPTION", "SEQ_ID", "PROBE_ID", "POSITION", "X",  "Y")]
        RG$G[, i] <- obj[, columns$f]
        RG$Gb[, i] <- obj[, columns$b]
        fullname <- files[i, 2]
        if (!is.null(path))
            fullname <- file.path(path, fullname)
        obj <- read.table(fullname, skip = skip, header = TRUE,
            sep = sep, quote = quote, check.names = FALSE, comment.char = "",
            fill = TRUE, nrows = nspots, ...)
        if (verbose)
            cat(paste("Read", fullname, "\n"))
        RG$R[, i] <- obj[, columns$f]
        RG$Rb[, i] <- obj[, columns$b]
        if (!is.null(wt.fun))
            RG$weights[, i] <- wt.fun(obj)
    }
    new("RGList", RG)
}#readNgIntensitiesTxt

readNblgnHeader <- function (file)
{
    firstfield <- scan(file, what = "", sep = "\t", quote = "\"",
        nlines = 100, flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE,
        multi.line = FALSE, allowEscape = FALSE)
    NHeaderRecords <- grep("# software=", firstfield)
    txt <- scan(file, what = "", sep = "\t", quote = "\"", nlines = NHeaderRecords -
        1, quiet = TRUE, allowEscape = FALSE)
    out <- list(NHeaderRecords = NHeaderRecords, BeginRawData = NHeaderRecords)
    out$Version <- txt[grep("version=", txt) + 1]
    out$Date <- txt[grep("date=", txt) + 1]
    out$ImageFile <- txt[grep("^Image File$", txt) + 1]
    out
}#readNimblegenHeader

## merge by column bind
merge.RGList <- function (x, y, ...)
{
    if (!is(y, "RGList"))
        stop("both x and y must be RGList objects")
    genes1 <- rownames(x$R)
    if (is.null(genes1))
        genes1 <- rownames(x$G)
    if (is.null(genes1))
        genes1 <- x$genes$ID
    genes2 <- rownames(y$R)
    if (is.null(genes2))
        genes2 <- rownames(y$G)
    if (is.null(genes2))
        genes2 <- y$genes$ID
    if (is.null(genes1) || is.null(genes2))
        stop("Need row names to align on")
    fields1 <- names(x)
    fields2 <- names(y)
    if (!identical(fields1, fields2))
        stop("The two RGLists have different components")
    ord2 <- match(makeUnique(genes1), makeUnique(genes2))
    #browser()
    data.fields <- grep("^[RG]b?$", fields1)
    # merge intensity information by 'cbind'
    for (i in data.fields) x[[i]] <- cbind(x[[i]], y[[i]][ord2, ])
    if ("targets" %in% fields1)
      x$targets <- rbind(x$targets, y$targets)
    return(x)
}#merge.RGList


joinNg <- 
function(RG1,RG2){
	YR <- YG <- matrix(0, dim(RG1)[[1]]+dim(RG2)[[1]], dim(RG1)[[2]])
	RG <- list(R = YR, G = YG, Rb = YR, Gb = YG)
	RG$R <- rbind( RG1$R[,RG1$targets$Array],RG2$R[,RG2$targets$Array] )
	RG$G <- rbind( RG1$G[,RG1$targets$Array],RG2$G[,RG2$targets$Array] )
	RG$Rb <- rbind( RG1$Rb[,RG1$targets$Array],RG2$Rb[,RG2$targets$Array] )
	RG$Gb <- rbind( RG1$Gb[,RG1$targets$Array],RG2$Gb[,RG2$targets$Array] )
	RG$genes <- rbind( RG1$genes[,RG1$targets$Array],RG2$genes[,RG2$targets$Array] )
	RG$targets <- RG1$targets
	new("RGList",RG)
}

preprocessNG <-
function (myRG, method = "vsn", returnMAList = FALSE, verbose = TRUE, 
    ...) 
{
    stopifnot(inherits(myRG, "RGList"))
    method <- match.arg(method, choices = c("vsn", "loess", "median", 
        "Gquantile", "nimblegen", "none"))
    if (method %in% c("loess", "median")) {
        if (verbose) 
            cat("Background correction...\n")
        myRG <- backgroundCorrect(myRG, method = "normexp", offset = 50)
    }
    if (verbose) 
        cat("Normalizing...\n")
    myMA <- switch(method, loess = normalizeWithinArrays(myRG, 
        method = "loess", ...), median = normalizeWithinArrays(myRG, 
        method = "median", ...), vsn = normalizeBetweenArrays(myRG, 
        method = "vsn", ...), Gquantile = normalizeBetweenArrays(myRG, 
        method = "Gquantile", ...), nimblegen = nimblegenNorm(myRG, 
        ...), none = normalizeWithinArrays(myRG, method = "none", 
        ...))
    if (returnMAList) {
        return(myMA)
    }
    else {
        return(Ringo:::asExprSet(myMA))
    }
}
