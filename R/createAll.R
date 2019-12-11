library(rmarkdown)
Rscript <- '/home/opt/R-3.4.0/bin/Rscript'

invisible(suppressWarnings(unlink(list.files(pattern = '\\.md$|\\.html$|_files$|^sites|^sample|\\.png$'), recursive = TRUE)))

invisible(sapply(list.files(pattern = '\\.Rmd$'), function(f){
  if(! file.exists(sub('\\.Rmd$', '.html', f)))
   system(paste0(Rscript, " -e 'library(rmarkdown); rmarkdown::render(\"", f, "\", \"html_document\")'"))
}))

system(paste0(Rscript, ' FigS3.R'))