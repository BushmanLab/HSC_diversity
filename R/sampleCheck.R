library(dplyr)
library(RMySQL)
library(gt23)
library(GenomicRanges)
library(parallel)

sampleList      <- readLines('../data/sampleList.20190718')
mergedSamples   <- read.table('../data/intSites.mergedSamples.collapsed.csv.gz', sep = ',', header = TRUE)
noMergedSamples <- read.table('../data/intSites.noMergedSamples.collapsed.csv.gz', sep = ',', header = TRUE)


n_distinct(subset(noMergedSamples, GTSP == 'GTSP0814')$posid)


prev.mergedSamples   <- read.table('../../HSC_diversity_data/intSites.mergedSamples.collapsed.csv.gz', sep = ',', header = TRUE)
prev.noMergedSamples <- read.table('../../HSC_diversity_data/intSites.noMergedSamples.collapsed.csv.gz', sep = ',', header = TRUE)

prev <- subset(prev.noMergedSamples, GTSP == 'GTSP0814')
new  <- subset(noMergedSamples, GTSP == 'GTSP0814')


GTSP0814.dbFrags <- getDBgenomicFragments(c('GTSP0814'), 
                                  sampleDB.group = 'specimen_management',
                                  intSiteDB.group = 'intsites_HSCdiversity') 



# Previous data archive method.
#----------------------------------------------------------

source('createIntSiteData.lib.R')
CPUs <- 40


# Retrieve intSite data from Bushman lab integration and sample databases.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
intSites <- getIntSiteData('specimen_management', 
                           'intsites_miseq', 
                           patients = c('pFR01', 'pFR03', 'pFR04', 'pFR05', 'pLT', 'pTST'))




length(GTSP0814.dbFrags)
length(subset(intSites, GTSP == 'GTSP0814'))


# Standardize cell type naming and subset the data to only incluse select cell types.
intSites$cellType <- toupper(as.character(intSites$cellType))
intSites$cellType <- gsub('\\s', '', intSites$cellType)
intSites$cellType <- gsub('\\-', '', intSites$cellType)
intSites$cellType <- gsub('NEUTROPHILS', 'GRANULOCYTES', intSites$cellType)
cellTypes <- c('GRANULOCYTES', 'TCELLS', 'MONOCYTES', 'BCELLS', 'NKCELLS', 'BM_BCELLS', 'BM_MONOCYTES', 'BM_GRANULOCYTES', 'BM_NKCELLS', 'BM_TCELLS')
intSites <- subset(intSites, cellType %in% cellTypes)



# Change patient codes.
convertName <- function(x){
  nameConversion <- data.frame(internal=c('pFR01',  'pFR03', 'pFR04',  'pFR05', 'pLT',    'pTST'),
                               external=c('WAS2',   'WAS4',  'WAS5',   'WAS7',  'bS/bS',  'b0/bE'))
  
  nameConversion[match(tolower(x), tolower(nameConversion$internal)),]$external
}

intSites$patient <- convertName(intSites$patient)

standardizeIntSites <- function(intSites){
  
  # Standardize break points by sample.
  cluster <- makeCluster(CPUs)
  intSites <- unlist(
    GRangesList(
      parLapply(cluster, split(intSites, intSites$sampleName), 
                function(p){
                  library(GenomicRanges)
                  library(gintools)
                  refine_breakpoints(p, counts.col = 'reads')
                })))
  
  
  # Standardize intSites by subject.
  intSites <- unlist(
    GRangesList(
      parLapply(cluster, split(intSites, intSites$patient), 
                function(p){
                  library(GenomicRanges)
                  library(gintools)
                  standardize_sites(p, counts.col = 'reads')
                })))
  
  stopCluster(cluster)
  
  intSites$posid <- paste0(seqnames(intSites), strand(intSites), start(flank(intSites, -1, start = T)))
  intSites
}


collapseReplicatesCalcAbunds <- function(intSites){
  cluster <- makeCluster(CPUs)
  o <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$GTSP), function(x){
    library(dplyr)
    library(GenomicRanges)
    x$sampleFragWidths <- paste(x$sampleName, width(x))
    x <- flank(x, -1, start = T)
    group_by(data.frame(x), posid) %>%
      mutate(reads = sum(reads)) %>%
      mutate(estAbund = n_distinct(sampleFragWidths)) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      select(seqnames, start, end, strand, posid, reads, patient, sampleName, GTSP, 
             cellType, timePoint, timePointDays, timePointMonths, estAbund) %>%
      mutate(relAbund = (estAbund / sum(estAbund))*100) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  })))
  
  stopCluster(cluster)
  o$sampleName <- NULL
  o
}



# Create an intsite object (framgent level) where no samples have been merged. 
intSites.noMergedSamples <- standardizeIntSites(intSites)
intSites.noMergedSamples.collapsed <- collapseReplicatesCalcAbunds(intSites.noMergedSamples)





tibble(prevData = n_distinct(prev$posid), 
       newData = n_distinct(new$posid), 
       reCalc = n_distinct(GTSP0814$posid),
       prevCalc = n_distinct(subset(intSites.noMergedSamples.collapsed, GTSP == 'GTSP0814')$posid))









b0/bE 	m11 	GRANULOCYTES 	GTSP0814 	32561 	34730 	2169

> all(sampleList %in% noMergedSamples$GTSP)
[1] TRUE
> sampleList %in% mergedSamples$GTSP
[1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[21]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[41]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[61]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE
[81]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE