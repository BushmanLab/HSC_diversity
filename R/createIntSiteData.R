library(GenomicRanges)
library(gt23)
library(dplyr)
options(stringsAsFactors = FALSE)

sampleList <- readLines('../data/sampleList.20190718')

# Define samples and retrieve genomic fragments.
dbConn  <- DBI::dbConnect(RMySQL::MySQL(), group = 'specimen_management')
samples <- DBI::dbGetQuery(dbConn, "select SpecimenAccNum from gtsp where Patient in ('pFR01', 'pFR03', 'pFR04', 'pFR05', 'pLT', 'pTST')")
intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq')
            


# Standardize cell type naming and subset the data to only incluse select cell types.
intSites$cellType <- gsub('\\s|\\-', '', toupper(intSites$cellType))
intSites$cellType <- gsub('NEUTROPHILS', 'GRANULOCYTES', intSites$cellType)
intSites <- subset(intSites, cellType %in% c('GRANULOCYTES', 'TCELLS', 'MONOCYTES', 'BCELLS', 'NKCELLS'))

stopifnot(all(sampleList %in% intSites$GTSP))

intSites <- subset(intSites, GTSP %in% sampleList)

# Change patient codes.
convertName <- function(x){
  nameConversion <- data.frame(internal=c('pFR01',  'pFR03', 'pFR04',  'pFR05', 'pLT',    'pTST'),
                               external=c('WAS2',   'WAS4',  'WAS5',   'WAS7',  'bS/bS',  'b0/bE'))
  
  nameConversion[match(tolower(x), tolower(nameConversion$internal)),]$external
}

intSites$patient <- convertName(intSites$patient)


# Standardize genomic fragment boundaries.
intSites <- stdIntSiteFragments(intSites) 


# Create an intsite object (framgent level) where no samples have been merged. 
intSites.noMergedSamples <- intSites
intSites.noMergedSamples.collapsed <- collapseReplicatesCalcAbunds(intSites.noMergedSamples) %>% annotateIntSites()


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

intSites.mergedSamples.collapsed <- collapseReplicatesCalcAbunds(intSites.mergedSamples) %>% annotateIntSites()
  
# Write out and compress intSite data.
write.table(data.frame(intSites.noMergedSamples), sep = ',', row.names = FALSE, col.names = TRUE, file = '../data/intSites.noMergedSamples.csv')
write.table(data.frame(intSites.noMergedSamples.collapsed), sep = ',', row.names = FALSE, col.names = TRUE, file = '../data/intSites.noMergedSamples.collapsed.csv')
write.table(data.frame(intSites.mergedSamples), sep = ',', row.names = FALSE, col.names = TRUE, file = '../data/intSites.mergedSamples.csv')
write.table(data.frame(intSites.mergedSamples.collapsed), sep = ',', row.names = FALSE, col.names = TRUE, file = '../data/intSites.mergedSamples.collapsed.csv')
system('gzip -9 ../data/intSites.noMergedSamples.csv')
system('gzip -9 ../data/intSites.noMergedSamples.collapsed.csv')
system('gzip -9 ../data/intSites.mergedSamples.csv')
system('gzip -9 ../data/intSites.mergedSamples.collapsed.csv')

# Write out a list of samples for SRA deposition.
write(unique(intSites$GTSP), file = '../data/sampleList')
