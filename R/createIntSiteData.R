library(RMySQL)
library(GenomicRanges)
library(parallel)
library(dplyr)
source('createIntSiteData.lib.R')
CPUs <- 40

# Retrieve intSite data from Bushman lab integration and sample databases.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
intSites <- getIntSiteData('specimen_management', 
                           'intsites_HSCdiversity', 
                            patients = c('pFR01', 'pFR03', 'pFR04', 'pFR05', 'pLT', 'pTST'))



# Standardize cell type naming and subset the data to only incluse select cell types.
intSites$cellType <- toupper(as.character(intSites$cellType))
intSites$cellType <- gsub('\\s', '', intSites$cellType)
intSites$cellType <- gsub('\\-', '', intSites$cellType)
intSites$cellType <- gsub('NEUTROPHILS', 'GRANULOCYTES', intSites$cellType)
cellTypes <- c('GRANULOCYTES', 'TCELLS', 'MONOCYTES', 'BCELLS', 'NKCELLS', 'BM_BCELLS', 'BM_MONOCYTES', 'BM_GRANULOCYTES', 'BM_NKCELLS', 'BM_TCELLS')
intSites <- subset(intSites, cellType %in% cellTypes)
table(unique(p$sampleName) %in% unique(intSites$GTSP))



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


nearestFeatures <- function(intSites){
  
  # Create a splitting vector for parallelization.
  intSites$s <- ntile(1:length(intSites), CPUs)
  
  cluster <- makeCluster(CPUs)
  
  # Add nearest gene and nearest oncogene annotations.
  names(intSites) <- NULL
  intSites  <- unlist(
    GRangesList(
      parLapply(cluster, split(intSites, intSites$s), 
                function(p){
                  source('../lib/Everett.R')
                  library(dplyr)
                  
                  # Nearest gene boundary
                  p$order <- 1:length(p)
                  p <- nearestGenomicFeature(p, subject = readRDS('../data/hg38.refSeqGenesGRanges.rds'), subject.exons = readRDS('../data/hg38.refSeqGenesGRanges.exons.rds'))
                  p <- p[order(p$order)]
                  
                  # Nearest oncogene
                  # Provide a list of concogenes to nearestGenomicFeature() and then extract 
                  # the metadata from the returned object, order it, and update the metadata 
                  # for the from the nearest feature object.
                  onogoGeneList <- readRDS('../data/humanOncoGenes.rds')
                  o <- nearestGenomicFeature(p, subject = readRDS('../data/hg38.refSeqGenesGRanges.rds'), subject.exons = readRDS('../data/hg38.refSeqGenesGRanges.exons.rds'), geneList = onogoGeneList)
                  om <- data.frame(mcols(o))
                  om <- om[order(om$order),]
                  
                  pm <- data.frame(mcols(p))
                  pm <- pm[order(pm$order),]
                  
                  pm$nearestOncoFeature       <- om$nearestFeature
                  pm$nearestOncoFeatureDist   <- om$nearestFeatureDist
                  pm$nearestOncoFeatureStrand <- om$nearestFeatureStrand
                  mcols(p) <- pm
                  
                  # Nearest lymphoma associated genes.
                  onogoGeneList <- readRDS('../data/humanLymphomaGenes.rds')
                  o <- nearestGenomicFeature(p, subject = readRDS('../data/hg38.refSeqGenesGRanges.rds'), subject.exons = readRDS('../data/hg38.refSeqGenesGRanges.exons.rds'), geneList = onogoGeneList)
                  om <- data.frame(mcols(o))
                  om <- om[order(om$order),]
                  
                  pm <- data.frame(mcols(p))
                  pm <- pm[order(pm$order),]
                  
                  pm$nearestlymphomaFeature       <- om$nearestFeature
                  pm$nearestlymphomaFeatureDist   <- om$nearestFeatureDist
                  pm$nearestlymphomaFeatureStrand <- om$nearestFeatureStrand
                  
                  mcols(p) <- pm
                  
                  p
                })))
  
  stopCluster(cluster)
  
  intSites$s <- NULL
  intSites
}

# Create an intsite object (framgent level) where no samples have been merged. 
intSites.noMergedSamples <- standardizeIntSites(intSites)
intSites.noMergedSamples.collapsed <- collapseReplicatesCalcAbunds(intSites.noMergedSamples)
intSites.noMergedSamples.collapsed <- nearestFeatures(intSites.noMergedSamples.collapsed)



# Now we create a similar object where biological replicates are merged.
intSites.mergedSamples <- intSites

# Neutrophils
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2265-1')]$sampleName <- 'GTSP2264-5'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2264-5')]$GTSP <- 'GTSP2264'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2265-2')]$sampleName <- 'GTSP2264-6'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2264-6')]$GTSP <- 'GTSP2264'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2265-3')]$sampleName <- 'GTSP2264-7'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2264-7')]$GTSP <- 'GTSP2264'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2265-4')]$sampleName <- 'GTSP2264-8'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2264-8')]$GTSP <- 'GTSP2264'

# Monocytes 
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2267-1')]$sampleName <- 'GTSP2266-5'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2266-5')]$GTSP <- 'GTSP2266'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2267-2')]$sampleName <- 'GTSP2266-6'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2266-6')]$GTSP <- 'GTSP2266'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2267-3')]$sampleName <- 'GTSP2266-7'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2266-7')]$GTSP <- 'GTSP2266'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2267-4')]$sampleName <- 'GTSP2266-8'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP2266-8')]$GTSP <- 'GTSP2266'

intSites.mergedSamples[which(intSites.mergedSamples$timePoint == 'M55.1')]$timePoint <- 'M55'
intSites.mergedSamples[which(intSites.mergedSamples$timePoint == 'M55')]$timePointDays <- 1672.918

# Granulocytes
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1294-1')]$sampleName <- 'GTSP1293-5'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1293-5')]$GTSP <- 'GTSP1293'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1294-2')]$sampleName <- 'GTSP1293-6'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1293-6')]$GTSP <- 'GTSP1293'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1294-3')]$sampleName <- 'GTSP1293-7'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1293-7')]$GTSP <- 'GTSP1293'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1294-4')]$sampleName <- 'GTSP1293-8'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1293-8')]$GTSP <- 'GTSP1293'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1295-1')]$sampleName <- 'GTSP1293-9'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1293-9')]$GTSP <- 'GTSP1293'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1295-2')]$sampleName <- 'GTSP1293-10'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1293-10')]$GTSP <- 'GTSP1293'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1295-3')]$sampleName <- 'GTSP1293-11'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1293-11')]$GTSP <- 'GTSP1293'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1295-4')]$sampleName <- 'GTSP1293-12'
intSites.mergedSamples[which(intSites.mergedSamples$sampleName == 'GTSP1293-12')]$GTSP <- 'GTSP1293'

intSites.mergedSamples <- standardizeIntSites(intSites.mergedSamples)
intSites.mergedSamples.collapsed <- collapseReplicatesCalcAbunds(intSites.mergedSamples)
intSites.mergedSamples.collapsed <- nearestFeatures(intSites.mergedSamples.collapsed)


# Write out and compress intSite data.
write.table(data.frame(intSites.noMergedSamples), sep = ',', row.names = FALSE, col.names = TRUE, file = '../data/intSites.noMergedSamples.csv')
write.table(data.frame(intSites.noMergedSamples.collapsed), sep = ',', row.names = FALSE, col.names = TRUE, file = '../data/intSites.noMergedSamples.collapsed.csv')
write.table(data.frame(intSites.mergedSamples), sep = ',', row.names = FALSE, col.names = TRUE, file = '../data/intSites.mergedSamples.csv')
write.table(data.frame(intSites.mergedSamples.collapsed), sep = ',', row.names = FALSE, col.names = TRUE, file = '../data/intSites.mergedSamples.collapsed.csv')
system('gzip -9 ../data/intSites.noMergedSamples.csv')
system('gzip -9 ../data/intSites.noMergedSamples.collapsed.csv')
system('gzip -9 ../data/intSites.mergedSamples.csv')
system('gzip -9 ../data/intSites.mergedSamples.collapsed.csv')

