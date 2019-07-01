getIntSiteData <- function(sampleDB.group, intSiteDB.group, patients=NULL,
                           samples=NULL, roundTimePoints=TRUE, runIDs=NULL){
  
  options(useFancyQuotes = FALSE)
  
  if(length(patients)>0 && length(samples)>0)
    stop('Both patients and samples can not be defined.')
  
  if(length(patients)==0 && length(samples)==0)
    stop('Either patients or samples must be defined.')
  
  dbConn1  <- DBI::dbConnect(RMySQL::MySQL(), group=sampleDB.group)
  dbConn2  <- DBI::dbConnect(RMySQL::MySQL(), group=intSiteDB.group)
  
  if (length(patients) > 0){
    selectString <- "select SpecimenAccNum from gtsp where Patient='%s'"
    samples <- unname(unlist(sapply(patients, function(p){
      unlist(DBI::dbGetQuery(dbConn1, sprintf(selectString, p)))
    })))
  }
  
  ## message('Specimen database samples: ', paste0(samples, collapse = ', '))
  
  if(length(samples) == 0) stop('Error: no samples have been selected.')
  
  intSiteSamples <- DBI::dbGetQuery(dbConn2, 'select * from samples')
  intSiteSamples$GTSP <- gsub('\\-\\d+$', '', intSiteSamples$sampleName)
  sampleIDs <- unique(base::subset(intSiteSamples, GTSP %in% samples)$sampleID)
  
  ## message('intSite samples: ', paste0(sampleIDs, collapse = ', '))
  
  if(length(sampleIDs) == 0) stop('Error: no intSite sample ids have been selected.')
  
  replicateQuery <- paste('samples.sampleID in (', paste0(sampleIDs, collapse = ','), ')')
  
  if(!is.null(runIDs)){
    runIDs <- paste0('(miseqid in (', paste0(sQuote(runIDs), collapse = ','), ')) and ')
  } else {
    runIDs = ''
  }
  
  q <- sprintf("select position, chr, strand, breakpoint, count,
               sampleName from sites left join samples on
               sites.sampleID = samples.sampleID
               left join pcrbreakpoints on
               pcrbreakpoints.siteID = sites.siteID
               where %s (%s)", runIDs, replicateQuery)
  
  options(warn=-1)
  sampleData <- DBI::dbGetQuery(dbConn1, "select * from gtsp")
  sites <- dbGetQuery(dbConn2, q)
  options(warn=0)
  
  if(nrow(sites) == 0) return(GRanges())
  
  sites$GTSP      <- as.character(sub('\\-\\d+', '', sites$sampleName))
  sites$patient   <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 'Patient']
  sites$timePoint <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 'Timepoint']
  sites$cellType  <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 'CellType']
  
  sites$timePoint <- toupper(sites$timePoint)
  sites$timePoint <- gsub('_', '.', sites$timePoint)
  
  sites$timePointType<- stringr::str_match(sites$timePoint, '[DMY]')
  sites$timePointType[which(is.na(sites$timePointType))] <- 'X'
  
  sites <- do.call(rbind, lapply(split(sites, sites$timePointType), function(x){
    n <- as.numeric(stringr::str_match(x$timePoint, '[\\d\\.]+'))
    
    if(x$timePointType[1] == 'D'){
      x$timePointMonths <- n / 30.4167
      x$timePointDays   <- n
    } else if(x$timePointType[1] == 'M'){
      x$timePointMonths <- n
      x$timePointDays   <- n * 30.4167
    } else if(x$timePointType[1] == 'Y'){
      x$timePointMonths <- n * 12
      x$timePointDays   <- n * 365
    } else {
      message('Warning - could not determine date unit for: ', x$timePointType[1])
      x$timePointMonths <- n
      x$timePointDays   <- n 
    }
    x
  }))
  sites$timePointType <- NULL
  
  if(roundTimePoints) sites$timePointMonths <- base::round(sites$timePointMonths)
  
  GenomicRanges::GRanges(seqnames=S4Vectors::Rle(sites$chr),
                         ranges=IRanges::IRanges(start=pmin(sites$position, sites$breakpoint),
                                                 end=pmax(sites$position, sites$breakpoint)),
                         strand=S4Vectors::Rle(sites$strand),
                         reads=sites$count,
                         patient=sites$patient,
                         sampleName=sites$sampleName,
                         GTSP=sites$GTSP,
                         cellType=sites$cellType,
                         timePoint=sites$timePoint,
                         timePointDays=sites$timePointDays,
                         timePointMonths=sites$timePointMonths)
}



nearestGenomicFeature <- function(query, subject = NULL, subject.exons = NULL, side='either', geneList=NULL){
  
  if(is.null(subject))       stop('subject parameter can not be NULL.')
  if(is.null(subject.exons)) stop('subject.exons parameter can not be NULL.')
  if(! is.null(geneList))    subject <- GenomicRanges::subset(subject, toupper(name2) %in% toupper(geneList))
  
  # If side is not set to either, collapse the subject ranges to single positions
  if (side %in% c("5p", "3p", "midpoint")) {
    options(warn=-1)
    if (side == "5p") subject <- GenomicRanges::flank(subject, width = -1)
    if (side == "3p") subject <- GenomicRanges::flank(subject, width = -1, start = FALSE)
    if (side == "midpoint") ranges(subject) <- IRanges(mid(ranges(subject)), width = 1)
    options(warn=0)
  }
  
  options(stringsAsFactors = FALSE)
  
  query.df  <- GenomicRanges::as.data.frame(query)
  subject.df <- GenomicRanges::as.data.frame(subject)
  
  query.df$strand <- as.character(query.df$strand)
  subject.df$strand <- as.character(subject.df$strand)
  
  subject.df$name2 <- as.character(subject.df$name2)
  
  
  subject.exons.df <- GenomicRanges::as.data.frame(subject.exons)
  query.df$inFeature            <- FALSE
  query.df$nearestFeature       <- 'None.found'
  query.df$nearestFeatureStrand <- 'None.found'
  query.df$inFeatureExon        <- FALSE
  query.df$inFeatureSameOrt     <- FALSE
  query.df$nearestFeatureStart  <- Inf
  query.df$nearestFeatureEnd    <- Inf
  query.df$nearestFeatureDist   <- Inf
  
  o  <- suppressWarnings(GenomicRanges::nearest(query, subject, select='all', ignore.strand=TRUE))
  
  if(length(o) > 0){
    
    createCol <- function(a, b, n){
      paste0(unique(cbind(a, b))[,n], collapse=',')
    }
    
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = createCol(subject.df[subjectHits,]$name2, subject.df[subjectHits,][['strand']], 1),
                    strand = createCol(subject.df[subjectHits,]$name2, subject.df[subjectHits,][['strand']], 2),
                    hitStart = min(subject.df[subjectHits,][['start']]),
                    hitEnd   = max(subject.df[subjectHits,][['end']])) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene, strand, hitStart, hitEnd) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    query.df[a$queryHits,]$nearestFeature       <- a$gene
    query.df[a$queryHits,]$nearestFeatureStrand <- a$strand
    query.df[a$queryHits,]$nearestFeatureStart  <- a$hitStart
    query.df[a$queryHits,]$nearestFeatureEnd    <- a$hitEnd
  }
  
  o <- suppressWarnings(GenomicRanges::findOverlaps(query, subject, select='all', ignore.strand=TRUE, type='any'))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = paste(unique(subject.df[subjectHits,]$name2), collapse=',')) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    query.df[a$queryHits,]$inFeature <- TRUE
  }
  
  o <- suppressWarnings(GenomicRanges::distanceToNearest(query,  subject, select='all', ignore.strand=TRUE))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::top_n(-1, distance) %>%
      dplyr::ungroup() %>%
      dplyr::select(queryHits, distance) %>% 
      dplyr::distinct() %>% 
      data.frame()
    query.df[a$queryHits,]$nearestFeatureDist <- a$distance
  }
  
  query.df$nearestFeatureBoundary <- ifelse(abs(query.df$start - query.df$nearestFeatureStart) > 
                                              abs(query.df$start - query.df$nearestFeatureEnd),   
                                            query.df$nearestFeatureEnd,  
                                            query.df$nearestFeatureStart)
  
  query.df$nearestFeatureDist <- query.df$nearestFeatureDist * sign(query.df$start - query.df$nearestFeatureBoundary)
  query.df$nearestFeatureDist <- ifelse(query.df$nearestFeatureStrand=='+', query.df$nearestFeatureDist, query.df$nearestFeatureDist * -1)
  
  query.df$nearestFeatureStart    <- NULL
  query.df$nearestFeatureEnd      <- NULL
  query.df$nearestFeatureBoundary <- NULL
  
  # In exon
  o <- suppressWarnings(GenomicRanges::findOverlaps(query, subject.exons, select='all', ignore.strand=TRUE, type='any'))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = paste(unique(subject.exons.df[subjectHits,]$name2), collapse=',')) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    #query.df[a$queryHits,]$inFeatureExon  <- a$gene
    query.df[a$queryHits,]$inFeatureExon  <- TRUE
  }
  
  # In TU ort
  # There may be cases where a site overlaps two or more features which have the sampe orientation.
  # ie. +,+,+ and we want to reduce these down to a single unique sign for comparison.
  a <- query.df[is.na(query.df$inFeature),]  
  b <- query.df[! is.na(query.df$inFeature),]
  
  if(nrow(a) > 0 && nrow(b) > 0){
    b$nearestFeatureStrandCmp <- unlist(lapply(strsplit(b$nearestFeatureStrand, ','), function(x){ paste(unique(x), collapse=',')}))
    b$inFeatureSameOrt <- b$strand == b$nearestFeatureStrandCmp
    b$nearestFeatureStrandCmp <- NULL
    query.df <- dplyr::bind_rows(a, b)
  }
  
  GenomicRanges::makeGRangesFromDataFrame(query.df, keep.extra.columns = TRUE)
}
