suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(magrittr)
  library(stats)
  library(gdata)
  library(posterior)
  library(bayesplot)
  library(jtools)
  library(cowplot)
  library(cmdstanr)
  library(ComplexHeatmap)
  library(PerformanceAnalytics)
  library(patchwork)
  library(ggpubr)
  library(circlize)
})

# load stan data and other parameters
celltypes <- snakemake@params[["celltypes"]]
#number.of.niches <- snakemake@params[["number_of_niches"]]
#nIter <- snakemake@params[["nIter"]]

celltypes <- gsub(" ", "_", celltypes)
celltypes <- gsub("-", "_", celltypes)
celltypes <- gsub(",", "", celltypes)

# LOAD niche factors in discovery and validation
niches.collapsed <- readRDS(snakemake@input[["microenvironment_niche_factors_collapsed"]])
niches.collapsed.validation <- readRDS(snakemake@input[["microenvironment_niche_factors_scored_val"]])

## add some dis/val information
rownames(niches.collapsed) <- paste("Discovery", rownames(niches.collapsed))
rownames(niches.collapsed.validation) <- paste("Validation", rownames(niches.collapsed.validation))

## scale data
niches.collapsed <- scale(niches.collapsed, center = TRUE, scale = TRUE)
niches.collapsed.validation <- scale(niches.collapsed.validation, center = TRUE, scale = TRUE)

## combine data
niches <- cbind(t(niches.collapsed), t(niches.collapsed.validation))
head(niches)

## scale data
#niches.to.plot <- t(scale(niches, center = TRUE, scale = TRUE))
niches.to.plot <- t(niches)
niches.to.plot <- niches.to.plot |> as.data.frame() |> rownames_to_column("niche")
head(niches.to.plot)

## save data
write.table(niches.to.plot, snakemake@output[["microenvironment_niche_factors_combined"]], sep = "\t", quote = FALSE, row.names = FALSE)

# load niche factor loadings in discovery and validation
niches.loadings.collapsed <- readRDS(snakemake@input[["niche_factor_loadings_collapsed"]]) |> as.data.frame() |> rownames_to_column("sample")
niches.loadings.collapsed.validation <- readRDS(snakemake@input[["niche_factor_loadings_scored_val"]]) |> as.data.frame() |> rownames_to_column("sample")

niches.loadings.collapsed <- niches.loadings.collapsed |> mutate(group = "Discovery")
niches.loadings.collapsed.validation <- niches.loadings.collapsed.validation |> mutate(group = "Validation")

niches.loadings.dis.val <- rbind(niches.loadings.collapsed, niches.loadings.collapsed.validation)

## save data
write.table(niches.loadings.dis.val, snakemake@output[["microenvironment_niche_factor_loadings_to_compare"]], sep = "\t", quote = FALSE, row.names = FALSE)

# LOAD intrinsice covariances in patients in discovery and validation 
cov.i.collapsed <- readRDS(snakemake@input[["intrinsic_covariance_matrices_collapsed"]])
cov.i.collapsed.validation <- readRDS(snakemake@input[["intrinsic_covariance_matrices_scored_val"]])

cov.i.cor.test.list <- lapply(celltypes, function(ct) {
  print(ct)
  print(dim(cov.i.collapsed[[paste0("cov_i_", ct)]]))
  print(dim(cov.i.collapsed.validation[[paste0("cov_i_", ct)]]))

  mat.dis <- cov.i.collapsed[[paste0("cov_i_", ct)]]
  mat.val <- cov.i.collapsed.validation[[paste0("cov_i_", ct)]]
  upper.tri.dis <- mat.dis[upper.tri(mat.dis, diag = FALSE)] |> as.vector()
  upper.tri.val <- mat.val[upper.tri(mat.val, diag = FALSE)] |> as.vector()
  
  corr.test <- cor.test(upper.tri.dis, upper.tri.val, method = "pearson", conf.level = 0.95)

  data.frame(
    cell_type = ct,
    corr = corr.test$estimate,
    p_value = corr.test$p.value,
    num_sig = nrow(mat.dis),
    conf.low = corr.test$conf.int[1],
    conf.high = corr.test$conf.int[2]
  )
})
cov.i.cor.test <- do.call(rbind, cov.i.cor.test.list)

## save data
write.table(cov.i.cor.test, snakemake@output[["intrinsic_covariance_matrices_correlation_data"]], sep = "\t", quote = FALSE, row.names = FALSE)

cov.i.ct.dis.val <- list(
  Discovery = cov.i.collapsed[[paste0("cov_i_", snakemake@params[["cov_celltype_to_plot"]])]],
  Validation = cov.i.collapsed.validation[[paste0("cov_i_", snakemake@params[["cov_celltype_to_plot"]])]]
)

saveRDS(cov.i.ct.dis.val, snakemake@output[["intrinsic_covariance_matrices_list_to_plot"]])