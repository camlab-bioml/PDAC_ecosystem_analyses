suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(SingleCellExperiment)
  library(tidyr)
  library(stringr)
})

sces <- readRDS(snakemake@input[['scelist_dis']])

# get gene universe
geneuniverse <- lapply(sces, rownames)
rm(sces)
geneuniverse <- Reduce(intersect, geneuniverse)
geneuniverse <- str_split(geneuniverse, "_", simplify = T)[,2]

saveRDS(geneuniverse, file = snakemake@output[['geneuniverse']])
