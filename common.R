#!/usr/bin/env Rscript
options(echo=TRUE)

loadpkg <- function(x){
    if(!suppressMessages(require(x, character.only=T, quietly=T, warn.conflicts=F))) {
        stop(paste0('Loading package ', x, ' failed!'))
    }
}

## load R packages, install if needed 
ris <- function(xs){
    for(x in xs){
        if(!suppressMessages(require(x, character.only=T, quietly=T, warn.conflicts=F))){
            install.packages(x, repos='http://cran.rstudio.com/') 
            loadpkg(x)
        }
    }
}
## load Bioconductor packages, install if needed 
bis <- function(xs){
    for(x in xs){
        if(!suppressMessages(require(x, character.only=T, quietly=T, warn.conflicts=F))){
            source('http://www.bioconductor.org/biocLite.R')
            biocLite(x)
            loadpkg(x)
        }
    }
}
