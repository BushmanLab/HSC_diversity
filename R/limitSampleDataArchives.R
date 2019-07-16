options(stringsAsFactors = FALSE)
sampleList <- readLines('../data/sampleList')
invisible(sapply(list.files('../data', pattern = 'intSites', full.names = TRUE), function(x){
  o <- read.table(x, sep = ',', header = TRUE)
  write.table(subset(o, GTSP %in% sampleList), sep = ',', col.names = TRUE, row.names = FALSE, file = sub('\\.gz$', '', x))
}))