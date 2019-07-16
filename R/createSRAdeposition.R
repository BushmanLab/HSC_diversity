library(RMySQL)
library(dplyr)

options(stringsAsFactors = FALSE)
sampleList <- readLines('../data/sampleList')
sampleSearchPath <- '/data/internal/geneTherapy/processedRuns/HPC_lineage'
subjectDetails <- read.table('../data/subjectDetails.tsv', sep = '\t', header = TRUE)

d <- system(paste0('find ', sampleSearchPath, ' -name GTSP* '), intern = TRUE)
d <- d[grepl('fastq', d)]
invisible(sapply(d, function(x) system(paste0('cp -r ', x, ' ../data/SRA/'))))

dbConn <- DBI::dbConnect(RMySQL::MySQL(), group = 'specimen_management')
sampleDB <- DBI::dbGetQuery(dbConn, 'select * from gtsp')

seqSamples <- unique(sapply(list.files('../data/SRA'), function(x) stringr::str_extract(x, 'GTSP\\d+')))

stopifnot(all(sampleList %in% seqSamples))



d <- bind_rows(lapply(sampleList, function(x){
  files <- list.files(path = '../data/SRA', pattern = x)
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

d$Patient <- convertName(d$Patient)
d$timePoint <- sub('_', '.', d$timePoint)
d$Sample_Name <- paste0(d$sample, '-', d$replicate)
d$isolate <- gsub('\\s', '_', paste0(d$Patient, '-', d$cellType, '-', d$timePoint, '-', d$replicate))
d$R1 <- paste0(d$sample, '-', d$replicate, '_R1.fastq.gz')
d$R2 <- paste0(d$sample, '-', d$replicate, '_R2.fastq.gz')

write.table(d, file = '../data/SRA.tsv', sep = '\t', col.names = TRUE, row.names = FALSE)
