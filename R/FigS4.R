library(dplyr)
library(vegan)
library(ggplot2)
options(stringsAsFactors = FALSE)
Rscript_path <- '/home/opt/R-3.4.0/bin/Rscript'

#
# (!) genomic heat maps required ~ 100 GB of memory.
#

intSites <- read.table('../data/intSites.mergedSamples.collapsed.csv.gz', sep = ',', header = TRUE)

# Restrict intSites to samples of interest.
intSites <- bind_rows(subset(intSites, patient == 'WAS5'  & timePoint %in% c('M12.6','M36', 'M43','M55')),
subset(intSites, patient == 'WAS4'  & timePoint %in% c('M12',  'M36', 'M48','M60')),
subset(intSites, patient == 'WAS2'  & timePoint %in% c('M22',  'M48', 'M78')),
subset(intSites, patient == 'WAS7'  & timePoint %in% c('M12',  'M30', 'M48')),
subset(intSites, patient == 'bS/bS' & timePoint %in% c('M12',  'M24')),
subset(intSites, patient == 'b0/bE' & timePoint %in% c('M11',  'M36', 'M48')))

intSites <- subset(intSites, cellType %in% c('GRANULOCYTES', 'MONOCYTES', 'BCELLS', 'NKCELLS', 'TCELLS'))

# Data report
# sampleSearchPath <- '/data/internal/geneTherapy/processedRuns/HPC_lineage'
# f <- list.files(path = sampleSearchPath, recursive = TRUE)
# f <- f[grep('GTSP', f)]
# r <- group_by(intSites, patient, timePoint, cellType) %>%
#      summarise(GTSP = GTSP[1], uniqueSites = n_distinct(posid), data = any(grepl(GTSP[1], f))) %>%
#      ungroup()



       

# CellType short hands 
cellTypeShortHand <- list('GRANULOCYTES' = 'G',  'BM_GRANULOCYTES' = 'G',
                          'TCELLS'       = 'T',  'BM_TCELLS'       = 'T',
                          'MONOCYTES'    = 'M',  'BM_MONOCYTES'    = 'M',
                          'BCELLS'       = 'B',  'BM_BCELLS'       = 'B',
                          'NKCELLS'      = 'K',  'BM_NKCELLS'      = 'K')

intSites$cellType2 <- sapply(as.character(intSites$cellType), function(x) cellTypeShortHand[[x]])


# Exclude bone marrow samples from the data set.
intSites <- intSites[! grepl('^BM', intSites$cellType),]


#--------------------------------------------------------------------------------------------------


intSites$refGenome <- 'hg38'
write(c('sampleName,GTSP,patient', 'G,G,xxx', 'M,M,xxx', 'B,B,xxx', 'K,K,xxx', 'T,T,xxx'), file = 'samples')
write('seqnames,strand,position,sampleName,refGenome', file = 'sites')

set.seed(42)
# /home/everett/manuscript.geneTherapy.HSCdiversity.new/data/tmp
s <- dplyr::sample_n(subset(intSites, seqnames %in% c(paste0('chr', 1:22), 'chrX', 'chrY')), 20000)
write.table(dplyr::select(s, seqnames, strand, start, cellType2,refGenome), 
            file = 'sites', col.names = FALSE, row.names = FALSE, sep = ',', append = TRUE, quote = FALSE)


comm <- paste(Rscript_path, '../software/genomicHeatmapMaker-from_input/genomic_heatmap_from_file.R samples ', 
              '-c ../software/genomicHeatmapMaker-from_input/INSPIIRED.yml ',
              '-o genomicHeatMapOutput ',
              '-f sites ',
              '-r hg38')

if(! dir.exists('genomicHeatMapOutput')) system(comm)


comm <- paste(Rscript_path, ' ../software/EpigeneticHeatmapMaker-from_input/epi_heatmap_from_file.R samples ', 
              '-c  ../software/EpigeneticHeatmapMaker-from_input/INSPIIRED.yml',
              '-o epiGeneticHeatMapOutput',
              '-t ../data/epiCellTypes',
              '-f sites ',
              '-r hg38')

if(! dir.exists('epiGeneticHeatMapOutput')) system(comm)



genomicHeatmap <- within(
  list(), {
    heatmap_sample_info <- read.csv('samples')
    gen_heatmap <- readRDS('genomicHeatMapOutput/roc.res.rds')
    
    heatmap_scale <- seq(0.2, 0.8, 0.1)
    gen_heatmap_colors <- colorspace::diverge_hsv(
      length(heatmap_scale), h = c(240, 0), v = 1, power = 1)
    
    select_gen_features <- row.names(gen_heatmap$ROC)
    select_gen_features <- c(
      "boundary.dist", "start.dist", "general.width", "gene.width",
      "within_refSeq_gene", "refSeq_counts.10k", "refSeq_counts.100k",
      "refSeq_counts.1M", "GC.100", "GC.1k", "GC.10k", "GC.100k", "GC.1M",
      "CpG_counts.1k", "CpG_counts.10k", "CpG_density.10k", "CpG_density.100k",
      "CpG_density.1M", "DNaseI_count.1k", "DNaseI_count.10k",
      "DNaseI_count.100k", "DNaseI_count.1M")
    
    gen_heatmap$ROC <- gen_heatmap$ROC[select_gen_features,]
    
    heatmap_sample_levels <- c("G", "M", "B", "K", "T")
    
    heatmap_figure_labels <- heatmap_sample_levels
    
    stat_cuts <- c(0, 0.001, 0.01, 0.05, 1)
    gen_comp_stats <- structure(cut(
      gen_heatmap$pvalues$op[select_gen_features, 1],
      stat_cuts,
      labels = c("***", " **", " * ", "   "),
      include.lowest = TRUE),
      names = select_gen_features)
    gen_row_names <- paste0(names(gen_comp_stats), " - ", gen_comp_stats)
    
    plot_data <- gen_heatmap$ROC %>%
      reshape2::melt() %>%
      mutate(
        feat = Var1,
        comp.sym = gen_comp_stats[Var1],
        Var1 = paste0(Var1, " - ", comp.sym),
        Var1 = factor(Var1, levels = gen_row_names),
        Var2 = factor(Var2, levels = heatmap_sample_levels),
        grp = " ",
        sig = as.vector(gen_heatmap$pvalues$np[select_gen_features,]),
        sym = cut(
          sig, stat_cuts, labels = c("***", " **", " * ", "   "),
          include.lowest = TRUE))
    
    levels(plot_data$Var2) <- heatmap_figure_labels
    
    plot_data$Var1 <- gsub('\\*', '', as.character(plot_data$Var1))
    plot_data$Var1 <- gsub('\\-', '', as.character(plot_data$Var1))
    
    
    levels_var1= gsub('\\*', '',gen_row_names)
    levels_var1= gsub('\\-', '',levels_var1)
    plot_data$Var1=factor(plot_data$Var1,levels=levels_var1)
    
    plot_data$sym  <- ''
    
    gen_plot <- ggplot(plot_data, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile(color = 'black') +
      geom_text(aes(label = sym), color = "black", size = 3, nudge_y = -0.15) +
      scale_x_discrete(position = "top") +
      scale_fill_gradient2(
        breaks = c(0.2, 0.4, 0.6, 0.8),
        low = gen_heatmap_colors[1],
        mid = gen_heatmap_colors[round(length(heatmap_scale)/2)],
        high = gen_heatmap_colors[length(heatmap_scale)],
        midpoint = 0.5) +
      guides(fill = guide_colorbar(
        title.position = "left", title.hjust = 0.5,
        direction = "horizontal")) +
      labs(x = NULL, y = NULL, fill = "ROC\nScore") +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text( angle = 0, hjust = 0, vjust = 0.5, size = 12),
        axis.text.x.top = element_text(
          angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        strip.placement = "outside",
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = nrow(gen_heatmap$ROC)/ncol(gen_heatmap$ROC))
  })

ggsave(genomicHeatmap$gen_plot, device = 'png', file = 'genomicPlot.png')


epiGenomicHeatmap <- within(
  list(), {
    heatmap_sample_info <- read.csv('samples')
    gen_heatmap <- readRDS('epiGeneticHeatMapOutput/roc.res.rds')
    
    heatmap_scale <- seq(0.2, 0.8, 0.1)
    gen_heatmap_colors <- colorspace::diverge_hsv(
      length(heatmap_scale), h = c(240, 115), v = 1, power = 1)
  
    select_gen_features <- row.names(gen_heatmap$ROC)
    
    gen_heatmap$ROC <- gen_heatmap$ROC[select_gen_features,]
    
    heatmap_sample_levels <- c("G", "M", "B", "K", "T")
    
    heatmap_figure_labels <- heatmap_sample_levels
    
    stat_cuts <- c(0, 0.001, 0.01, 0.05, 1)
    gen_comp_stats <- structure(cut(
      gen_heatmap$pvalues$op[select_gen_features, 1],
      stat_cuts,
      labels = c("***", " **", " * ", "   "),
      include.lowest = TRUE),
      names = select_gen_features)
    gen_row_names <- paste0(names(gen_comp_stats), " - ", gen_comp_stats)
    
    plot_data <- gen_heatmap$ROC %>%
      reshape2::melt() %>%
      mutate(
        feat = Var1,
        comp.sym = gen_comp_stats[Var1],
        Var1 = paste0(Var1, " - ", comp.sym),
        Var1 = factor(Var1, levels = gen_row_names),
        Var2 = factor(Var2, levels = heatmap_sample_levels),
        grp = " ",
        sig = as.vector(gen_heatmap$pvalues$np[select_gen_features,]),
        sym = cut(
          sig, stat_cuts, labels = c("***", " **", " * ", "   "),
          include.lowest = TRUE))
    
    levels(plot_data$Var2) <- heatmap_figure_labels
    
    plot_data$Var1 <- gsub('*', '', plot_data$Var1)
    
    
    levels_var1= gsub('\\*', '',gen_row_names)
    levels_var1= gsub('\\-', '',levels_var1)
    plot_data$Var1=factor(plot_data$Var1,levels=levels_var1)
    
    plot_data$sym  <- ''
    
    gen_plot <- ggplot(plot_data, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile(color = 'black') +
      geom_text(aes(label = sym), color = "black", size = 3, nudge_y = -0.15) +
      scale_x_discrete(position = "top") +
      scale_fill_gradient2(
        breaks = c(0.2, 0.4, 0.6, 0.8),
        low = 'blue3',
        mid = 'white',
        high = 'yellow3',
        midpoint = 0.5) +
      guides(fill = guide_colorbar(
        title.position = "left", title.hjust = 0.5,
        direction = "horizontal")) +
      labs(x = NULL, y = NULL, fill = "ROC\nScore") +
      # custom_theme +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text( angle = 0, hjust = 0, vjust = 0.5, size = 12),
        axis.text.x.top = element_text(
          angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        strip.placement = "outside",
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = nrow(gen_heatmap$ROC)/ncol(gen_heatmap$ROC))
  })

ggsave(epiGenomicHeatmap$gen_plot, device = 'png', file = 'epiGenPlot.png')



# Oncogene enrichment
#--------------------------------------------------------------------------------------------------

# Restrict intSites to samples of interest.
early <- bind_rows(subset(intSites, patient == 'WAS5' & timePoint == 'Y1.6'),
                   subset(intSites, patient == 'WAS4'  & timePoint == 'Y1'),
                   subset(intSites, patient == 'WAS2'  & timePoint == 'M22'),
                   subset(intSites, patient == 'WAS7'  & timePoint == 'M30'),
                   subset(intSites, patient == 'bS/bS' & timePoint == 'Y1'),
                   subset(intSites, patient == 'b0/bE' & timePoint == 'M11'))

late <- bind_rows(subset(intSites,  patient == 'WAS5' & timePoint == 'M55'),
                   subset(intSites, patient == 'WAS4'  & timePoint == 'M48'),
                   subset(intSites, patient == 'WAS2'  & timePoint == 'M78'),
                   subset(intSites, patient == 'WAS7'  & timePoint == 'M48'),
                   subset(intSites, patient == 'bS/bS' & timePoint == 'M24'),
                   subset(intSites, patient == 'b0/bE' & timePoint == 'M48'))

m <- matrix(c(n_distinct(subset(early, abs(nearestOncoFeatureDist) <= 50000)$posid),
              n_distinct(subset(early, abs(nearestOncoFeatureDist) >  50000)$posid),
              n_distinct(subset(late,  abs(nearestOncoFeatureDist) <= 50000)$posid),
              n_distinct(subset(late,  abs(nearestOncoFeatureDist) >  50000)$posid)),
              byrow = TRUE, ncol = 2, dimnames = list(c('Early', 'Late'), c('< 50KB', '> 50KB')))
m
fisher.test(m)


onocGenes <- readRDS('../data/humanOncoGenes.rds')

sharedNearestGene <- base::intersect(early$nearestFeature, late$nearestFeature)

maxAbund <- 
  dplyr::group_by(intSites, nearestFeature) %>%
  dplyr::group_by(GTSP) %>%
  dplyr::mutate(sampleFrags = sum(estAbund)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(sampleFrags >= 100) %>%
  dplyr::group_by(nearestFeature) %>%
  dplyr::summarise(maxRelAbund = max(relAbund), maxAbund = max(estAbund)) %>%
  dplyr::ungroup()

early <- 
  dplyr::filter(early, nearestFeature %in% sharedNearestGene) %>%
  dplyr::mutate(earlySites = n_distinct(posid)) %>%
  dplyr::group_by(nearestFeature) %>%
  dplyr::summarise(earlySites = earlySites[1], earlyFreq = n_distinct(posid) / earlySites[1]) %>%
  dplyr::ungroup() %>% 
  dplyr::select(earlySites, nearestFeature, earlyFreq)


late <- 
  dplyr::filter(late, nearestFeature %in% sharedNearestGene) %>%
  dplyr::mutate(lateSites = n_distinct(posid)) %>%
  dplyr::group_by(nearestFeature) %>%
  dplyr::summarise(lateSites = lateSites[1], lateFreq = n_distinct(posid) / lateSites[1]) %>%
  dplyr::ungroup() %>% 
  dplyr::select(lateSites, nearestFeature, lateFreq)

genesToLabel <- c('PBX3', 'NEAT1', 'SETD2', 'RERE', 'GPATCH8')

o <- dplyr::left_join(early, late, by = "nearestFeature") %>%
     dplyr::mutate(geneLabel = ifelse(earlyFreq > 0.003 | nearestFeature %in% genesToLabel, nearestFeature, ''))

o$onco <- toupper(o$nearestFeature) %in% toupper(onocGenes)
o$maxAbund <- maxAbund[match(o$nearestFeature, maxAbund$nearestFeature),]$maxRelAbund
o$geneLabel <- sub('LOC101929163,C6orf10', 'C6orf10', o$geneLabel)
o$onco <- ifelse(o$onco, 'Oncogene', 'Non-oncogene')

library(ggrepel)
ppNum <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


bivariatePlot <-
  ggplot(o, aes(earlyFreq, lateFreq, color = onco, size = maxAbund)) +
  theme_bw() +
  geom_point() + 
  scale_color_manual(name = 'Gene category', values = c('black', 'red')) +
  scale_size(name = 'Maximum relative abundance',  range = c(0.25, 4)) + 
  geom_abline(slope=1, intercept=0, color='black', size=0.5) +
  geom_text_repel(aes(label=geneLabel), color='black', size=3,  direction='y', box.padding=0.5, point.padding=1.5, ylim = c(0.001, 0.006)) +
  xlim(c(0, 0.007)) +
  ylim(c(0, 0.007)) +
  labs(x = paste0('Earliest time point integration frequency (', ppNum(o$earlySites[1]), ' sites)'),
       y = paste0('Latest time point integration frequency (', ppNum(o$lateSites[1]), ' sites)')) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))

ggsave(bivariatePlot, device = 'png', file = 'bivariatePlot.png')



# Pooled Chao1 estimates.
#--------------------------------------------------------------------------------------------------

# For each patient / timePoint combination, calculate Chao1 of pooled samples.
d1 <- dplyr::group_by(intSites, patient, timePoint) %>%
      dplyr::summarise(Chao1 = round(estimateR(estAbund, index='chao')[2], 0)) %>%
      dplyr::ungroup()

# For each patient / timePoint combination, calculate Chao1 of pooled GRANULOCYTES and MONOCYTES samples.
d2 <- dplyr::filter(intSites, cellType %in% c('GRANULOCYTES', 'MONOCYTES')) %>%
      dplyr::group_by(patient, timePoint) %>%
      dplyr::summarise(Chao1 = round(estimateR(estAbund, index='chao')[2], 0)) %>%
      dplyr::ungroup()

# Foreach patient, calculate Chao1 for all samples >= 24 months.
d3 <- dplyr::filter(intSites, timePointMonths >=  24) %>%
      dplyr::group_by(patient) %>%
      dplyr::summarise(Chao1 = round(estimateR(estAbund, index='chao')[2], 0)) %>%
      dplyr::ungroup()
  
# Foreach patient, calculate Chao1 for all samples >= 24 months from pooled GRANULOCYTES and MONOCYTES samples.
d4 <- dplyr::filter(intSites, timePointMonths >=  24 & cellType %in% c('GRANULOCYTES', 'MONOCYTES')) %>%
      dplyr::group_by(patient) %>%
      dplyr::summarise(Chao1 = round(estimateR(estAbund, index='chao')[2], 0)) %>%
      dplyr::ungroup()
