
########################################################################
######################## DATA PREPARATION ##############################
########################################################################

###Loading the packages
library(plyr)
library(tidyverse)
library(reshape2)
library(pander)
library(sqldf)
library(lazyeval)
library(lme4)
library(psych)
library(MASS)
library(scales)
library(gridExtra)
library(ggtern)
library(grid)
library(knitr)
library(kableExtra)
library(gtools)
library(ggalluvial)


#####################################
### Data loading and formatting  ###
####################################


intSites <- read.table(gzfile('../data/intSites.mergedSamples.collapsed.csv.gz'), sep = ',', header = TRUE)


for(i in c('start', 'end', 'width', 'estAbund')) intSites[[i]] <- as.numeric(intSites[[i]]) 


#
# Read in the cell sorting cross over reports which are stored in a single file with records separated with '#%'.
# The counts table is identifiable by the key word COUNTS. 
#

options(stringsAsFactors = FALSE)
crossOverReports <- readChar('../data/crossOverReports.tsv', file.info('../data/crossOverReports.tsv')$size)
crossOverReports <- unlist(strsplit(crossOverReports, '#%'))
crossOverReports <- unlist(lapply(crossOverReports, function(x){
  source     <- ifelse(grepl('#source:', x), sub('\\t+$', '', str_match(x, '#source:\\s+(.+)')[2]), '')
  cellCounts <- as.numeric(str_match_all(str_match(x, 'initialCellCounts,(.+)')[2], '([\\d\\.]+)')[[1]][,2])
  psort <- as.numeric(str_match_all(str_match(x, 'pB_postsort,(.+)')[2], '([\\d\\.]+)')[[1]][,2])
  count_psort <-as.numeric(str_match_all(str_match(x, 'postsortCount,(.+)')[2], '([\\d\\.]+)')[[1]][,2])
  VCN <-as.numeric(str_match_all(str_match(x, 'VCN,(.+)')[2], '([\\d\\.]+)')[[1]][,2])
  DNA <-as.numeric(str_match_all(str_match(x, 'DNA,(.+)')[2], '([\\d\\.]+)')[[1]][,2])
  patient    <- str_match(x, '#(\\S+)')[2]
  timePoint  <- str_match(x, '#\\S+\\s+(\\S+)')[2]
  t          <- gsub('\\t+\\n', '\n', substr(x, regexpr('COUNTS', x)+7, nchar(x)))
  m          <- read.table(tc <- textConnection(t), header = TRUE, fill = TRUE, check.names = FALSE, sep='\t'); close(tc); 
  if(nrow(m) < 5 | length(m) < 5) return(NA)
  m <- m[1:5, 1:5]
  r <- list()
  r[[paste(patient, timePoint, source, sep='|')]] <- list(
    patient   = str_match(x, '#(\\S+)\\s+(\\S+)')[2],
    timePoint = str_match(x, '#(\\S+)\\s+(\\S+)')[3],
    source    = source,
    t         = t,
    table     = m,
    cellCount = cellCounts,
    psort     = psort,
    count_psort = count_psort,
    VCN = VCN,
    DNA = DNA)
  
  message('Read in ', names(r))
  r
}), recursive = FALSE)

crossOverReports <- crossOverReports[sapply(crossOverReports, is.list)]

#
# Subset the resulting list to include only those reports which will be used in the study.
#


all_samples <- c('WAS4|m48|Blood',
                 'WAS5|m55|Blood',
                 'bS/bS|m24|Blood',
                 'b0/bE|m48|Blood')


if(! all(all_samples %in% names(crossOverReports))) stop('All the requested crossover tables were not found.')
crossOverReports <- crossOverReports[all_samples]  


### proportion for each cell type
prop_HD=c(0.2,0.2,0.2,0.2,0.2)


# Create a key to be used to update corrected abundances.
intSites$timePoint=tolower(intSites$timePoint)
intSites$key <- paste(intSites$posid, intSites$cellType, intSites$patient, intSites$timePoint)

# Put the data in a wide format count table.
assign(paste('liste_intSites_new_sup',0,sep=''),lapply(crossOverReports,function(x){  
  
  # Subset the intSite data to include only sites from a specific cell type and time point.
  i <- which(intSites$patient==x$patient &
               intSites$timePoint==x$timePoint &
               toupper(intSites$cellType) %in% toupper(names(x$table)))
  d <- intSites[i,]
  
  if(nrow(d) == 0) stop(paste0(x$patient, ' / ', x$timePoint, ' could not be found in the intSite data.'))
  
  message('Cell types in retrieved data subset(', x$patient, ' - ', x$timePoint, '): ',
          paste0(unique(d$celltype), collapse=', '))
  
  # Replace NA gene names with 'NONE' in case a function is sensitive to NA.
  if(length(which(is.na(d$nearestFeature))) > 0) d[which(is.na(d$nearestFeature)),]$nearestFeature <- 'NONE'
  
  # Reorganize the data to create an intSite / cell count table.
  d2 <- reshape2::dcast(d, posid ~ cellType, value.var='estAbund', fun.aggregate=function(x){x[1]}, fill=0)  
  
  # Add missing cell types.
  d2[names(x$table)[! toupper(names(x$table)) %in% toupper(names(d2))]] <- 0
  
  # Add nearest gene column.
  d2$gene <- intSites[match(d2$posid, intSites$posid),]$nearestFeature  
  
  # Add inFeature column.
  d2$inFeature <- intSites[match(d2$posid, intSites$posid),]$inFeature  
  
  # Reorganize the column headers to match Correction_CutData_new() input structure.
  d2 <- d2[,c(1,grep('gene',   names(d2), ignore.case = TRUE),
              grep('BCELL',  names(d2), ignore.case = TRUE),
              grep('MONO',   names(d2), ignore.case = TRUE),
              grep('GRANULO', names(d2), ignore.case = TRUE),
              grep('NKCELL', names(d2), ignore.case = TRUE),
              grep('TCELL',  names(d2), ignore.case = TRUE),
              grep('inFeature',  names(d2), ignore.case = TRUE)
  )]
  
  # Create a cell count column.
  d2$cellCount <-  apply(d2, 1, function(x){ sum(as.integer(x[3:7])) })
  
  # Rename the input columns.
  names(d2)=c("integrationSite", "gene", "Bcells", "Monocytes", "Granulocytes", "NKcells", "Tcells","inFeature", "cellCount")
  
  d2$integrationSite=as.character(d2$integrationSite)
  d2$gene=as.character(d2$gene)
  d2$key <- paste(c(x$patient, x$timePoint,x$source),collapse = "|") 
  intSites_new=d2
  
  return(intSites_new)
} ))


# Gather all the data sup0 (without threshold) in a unique dataframe.
intSites_new_concat_sup0=do.call(rbind,liste_intSites_new_sup0)


write.table(intSites_new_concat_sup0,'intSites_new_concat_sup0.csv',row.names = F,quote = F,sep=';')
