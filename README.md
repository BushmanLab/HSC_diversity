# HSC_diversity

The present repository contains integration sites data as well as codes to generate figures 
and tables found in the article “Clonal tracking in gene therapy patients reveals a diversity 
of human hematopoietic differentiation programs”, by Six E et al.

Final integration sites dataframe can be found in data/intSites.mergedSamples.collapsed.csv.gz

All sequence data used in this study is available at the NCBI SRA under SRP139090.

Execute the following command to create all figures and tables found in the article.  
Before running the script, please update the top of the script to include the path 
to your local R installation which contains the required libraries (listed below).

```
Rscript R/createAll.R
```

* rmarkdown
* parallel
* tidyverse
* reshape2
* psych
* MASS
* scales
* gridExtra
* grid
* knitr
* vegan
* kableExtra
* DataCombine
* pander
* RColorBrewer
* arrangements
* gtable
* sqldf
* lazyeval
* lme4
* ggtern
* gtools
* ggalluvial
* gplots
* ggrepel
