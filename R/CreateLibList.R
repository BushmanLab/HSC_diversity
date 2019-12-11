write(file = 'requiredLibs', unique(unlist(sapply(list.files(pattern = '*.Rmd$|*.R$'), function(x){
 unlist(lapply(stringr::str_match_all(readLines(x), 'library\\(([^\\)]+)'), '[', 2))
}))))