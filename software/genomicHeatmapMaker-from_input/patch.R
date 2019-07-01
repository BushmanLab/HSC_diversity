makeGRanges <- function (x, freeze = NULL, positionsOnly = FALSE, soloStart = FALSE, 
          chromCol = NULL, strandCol = NULL, startCol = NULL, stopCol = NULL, 
          keepFactors = FALSE) 
{
  if (is.null(chromCol)) {
    colIndex <- getRelevantCol(names(x), c("chr", "chromosome", 
                                           "tname", "space", "chrom", "contig", "seqnames"), 
                               "seqnames")
    chromCol <- names(x)[colIndex]
  }
  x$seqnames <- x[, chromCol]
  if (is.null(strandCol)) {
    colIndex <- getRelevantCol(names(x), c("ort", "orientation", 
                                           "strand"), "strand")
    strandCol <- names(x)[colIndex]
  }
  x$strand <- x[, strandCol]
  if (is.null(startCol)) {
    startCol <- getRelevantCol(names(x), c("position", "intsite", 
                                           "txstart", "start", "chromstart"), "start", multiple.ok = TRUE)
    startCol <- names(x)[startCol[1]]
  }
  x$start <- x[, startCol]
  if (!as.logical(soloStart) & is.null(stopCol)) {
    stopCol <- getRelevantCol(names(x), c("txend", "end", 
                                          "stop", "chromend"), "end", multiple.ok = TRUE)
    stopCol <- names(x)[stopCol[1]]
  }
  x$end <- x[, stopCol]
  if (any(is.na(x$start))) {
    stop("NAs found in column containing start positions")
  }
  if (any(is.na(x$seqnames))) {
    stop("NAs found in column containing chromosome information")
  }
  if (!keepFactors) {
    factorCols <- sapply(x, is.factor)
    if (any(factorCols)) {
      for (y in names(which(factorCols))) {
        x[, y] <- as.character(x[, y])
        if (!any(is.na(suppressWarnings(as.numeric(x[, 
                                                     y]))))) {
          x[, y] <- as.numeric(x[, y])
        }
      }
    }
  }
  if (length(startCol) > 0 & length(stopCol) > 0) {
    if (any(is.na(x$end))) {
      stop("NAs found in column containing end positions")
    }
  }
  else {
    x$end <- x$start
  }
  if (as.logical(positionsOnly)) {
    x <- x[, c("seqnames", "start", "end", "strand")]
  }
  metadataCols <- setdiff(names(x), c("seqnames", "start", 
                                      "end", "strand", "width"))
  metadataCols <- metadataCols[!is.na(metadataCols)]
  sites.gr <- GRanges(seqnames = x$seqnames, IRanges(start = x$start, 
                                                     end = x$end), strand = x$strand)
  for (f in metadataCols) {
    mcols(sites.gr)[[f]] <- x[, f]
  }
  if (!is.null(freeze)) {
    genomeLib <- grep(freeze, BSgenome::installed.genomes(), 
                      value = TRUE)
    if (length(genomeLib) != 0) {
      bsGenomeObject <- strsplit(genomeLib, "\\.")[[1]][2]
      chrom.info <- seqlengths(do.call(`:::`, list(genomeLib, 
                                                   bsGenomeObject)))
    }
    else {
      z <- gzcon(url(paste0("http://hgdownload.cse.ucsc.edu/goldenPath/", 
                            freeze, "/database/chromInfo.txt.gz")))
      zlines <- try(readLines(z))
      close(z)
      if (class(zlines) == "try-error") 
        stop("Could not get thru to UCSC server -\n               try later or drop the freeze parameter!")
      raw.data <- textConnection(zlines)
      chrom.info <- read.delim(raw.data, header = FALSE, 
                               stringsAsFactors = FALSE)[, 1:2]
      chrom.info <- structure(chrom.info$V2, names = chrom.info$V1)
      close(raw.data)
    }
    
    sites.gr <- subset(sites.gr, seqnames(sites.gr) %in% names(chrom.info))
    
    seqlevels(sites.gr) <- sortSeqlevels(seqlevelsInUse(sites.gr))
    seqlengths(sites.gr) <- chrom.info[seqlevels(sites.gr)]
  }
  sites.gr <- sort(sites.gr, ignore.strand = TRUE)
  sites.gr
}
