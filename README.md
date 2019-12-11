# HSC_diversity

The present repository contains integration sites data as well as codes to generate figures 
and tables found in the article “Clonal tracking in gene therapy patients reveals a diversity 
of human hematopoietic differentiation programs”, by Six E et al.


## IS dataframes ##
Final integration sites dataframe can be found in data/intSites.mergedSamples.collapsed.csv.gz

The IS dataframe containing technical replicates data can be found in data/intSites.noMergedSamples.csv.gz

The IS dataframe containing biological replicates data can be found in data/intSites.nomergedSamples.collapsed.csv.gz

All sequence data used in this study is available at the NCBI SRA under SRP139090.

## Samples Metadata ##
For each patient and time point the cell contamination matrix, VCN input and DNA input can be found in data/crossOverReports.tsv

## Clonal Tracking pipeline ##
The IS data are corrected and normalized using the function R/AbundanceCorrection_2steps_Normalization_Steps.R

To analyze HSPC lineage output, the corrected data are clustered using the function R/Kmeans_clustering.R

## Clonal Tracking modelling algorithm ##
The function allowing to run the simulation algorithm can be found in R/Modelling/Data_generator.R

It has been used to model two distinct types of HSPC population
- a homogeneous HSPC population R/Modelling/FigS16_homogeneous.Rmd
- a heterogeneous HSPC population R/Modelling/FigS16_heterogeneous.Rmd

## Code and libraries ##
Execute the following command to create all figures and tables found in the article.  
Before running the script, please update the top of the script to include the path 
to your local R installation which contains the required libraries (listed below).

```
Rscript R/createAll.R
```

 "arrangements"
 "circlize" 
 "ComplexHeatmap" 
 "DataCombine" 
 "ggalluvial"
 "ggrepel" 
 "ggtern" 
 "gplots" 
 "grid" 
 "gridExtra"
 "gtable" 
 "gtools" 
 "kableExtra" 
 "knitr" 
 "lazyeval"
 "lme4" 
 "MASS"
 "pander" 
 "parallel" 
 "plyr"
"psych" 
"RColorBrewer" 
"reshape2" 
"rmarkdown" 
"scales"
 "sqldf" 
 "tidyverse" 
 "vegan" 