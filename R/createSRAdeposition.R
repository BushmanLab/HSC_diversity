library(RMySQL)
library(dplyr)

options(stringsAsFactors = FALSE)
sampleList <- readLines('../data/sampleList.20190718')
sampleSearchPath <- '/data/internal/geneTherapy/processedRuns/HPC_lineage'
subjectDetails <- read.table('../data/subjectDetails.tsv', sep = '\t', header = TRUE)

d <- system(paste0('find ', sampleSearchPath, ' -name GTSP* '), intern = TRUE)
d <- d[grepl('fastq', d)]

seqSamples <- unique(sapply(d, function(x) stringr::str_extract(x, 'GTSP\\d+')))
stopifnot(all(sampleList %in% seqSamples))


dbConn <- DBI::dbConnect(RMySQL::MySQL(), group = 'specimen_management')
sampleDB <- DBI::dbGetQuery(dbConn, 'select * from gtsp')

d <- bind_rows(lapply(sampleList, function(x){
  files <- list.files(path = sampleSearchPath, pattern = x, recursive = TRUE)
  files <- files[grepl('_R1', files)]
  d <- tibble(sample = x, replicate = stringr::str_match(files, '\\-(\\d+)_')[,2])
  d$Patient <- sampleDB[match(x, sampleDB$SpecimenAccNum),]$Patient
  d$cellType <- sampleDB[match(x, sampleDB$SpecimenAccNum),]$CellType
  d$timePoint <- sampleDB[match(x, sampleDB$SpecimenAccNum),]$Timepoint
  d$age <- subjectDetails[match(d$Patient, subjectDetails$Patient),]$age
  d
}))

convertName <- function(x){
  nameConversion <- data.frame(internal=c('pFR01',  'pFR03', 'pFR04',  'pFR05', 'pLT',    'pTST'),
                               external=c('WAS2',   'WAS4',  'WAS5',   'WAS7',  'bS/bS',  'b0/bE'))
  
  nameConversion[match(tolower(x), tolower(nameConversion$internal)),]$external
}


# Clean up and fill out SRA data frame.
d$Patient <- convertName(d$Patient)
d$timePoint <- sub('_', '.', d$timePoint)
d$Sample_Name <- paste0(d$sample, '-', d$replicate)
d$isolate <- gsub('\\s', '_', paste0(d$Patient, '-', d$cellType, '-', d$timePoint, '-', d$replicate))
d$R1 <- paste0(d$sample, '-', d$replicate, '_R1.fastq.gz')
d$R2 <- paste0(d$sample, '-', d$replicate, '_R2.fastq.gz')


# Write out data which will be used to create SRA upload spread sheets.
write.table(d, file = '../data/SRA.tsv', sep = '\t', col.names = TRUE, row.names = FALSE)


# Find the paths to the sequencing files and copy they files to a common directory for SRA upload.
R1paths <- unname(sapply(d$R1, function(x) system(paste0('find ', sampleSearchPath, ' -name ', x), intern = TRUE)))
R2paths <- sub('_R1\\.', '_R2.', R1paths)
invisible(sapply(c(R1paths, R2paths), function(x) system(paste0('cp -r ', x, ' ../data/SRA/'))))

