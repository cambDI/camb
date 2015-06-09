<<include=FALSE>>=
  opts_chunk$set(concordance=TRUE)
opts_chunk$set(dev='pdf')
Sys.setenv(TEXINPUTS=getwd(),
           BIBINPUTS=getwd(),
           BSTINPUTS=getwd())
@

<<echo=FALSE,results='hide'>>=
  options(width=60)
suppressWarnings(library(png,quietly=TRUE,warn.conflicts=FALSE))
suppressWarnings(library(caret,quietly=TRUE,warn.conflicts=FALSE))
suppressWarnings(library(doMC,quietly=TRUE,warn.conflicts=FALSE))
library(ggplot2,quietly=TRUE)
@