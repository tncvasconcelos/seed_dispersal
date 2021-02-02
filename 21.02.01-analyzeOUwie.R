# this script will load, summarize, and analyze the results of 1000 OUwie runs on our simmaps 

## functions

## imports
require(corHMM)
require(OUwie)

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
Rsaves <- paste0(wd, "/res_ouwie/", dir("res_ouwie/"))
labels <- unlist(lapply(strsplit(dir("res_ouwie/"), "-"), function(x) paste0(x[1], x[2])))


