suppressPackageStartupMessages({
        library(magrittr)
        library(tidyverse)
        library(sjstats)
        library(ggpubr)
        library(ggsci)
        library(stringr)
        library(ComplexHeatmap)
        library(Matrix)
        library(here)
        library(scales)
        library(gridExtra)
        library(bluster)
	library(scran)
})

set.seed(123L)

# load tumour sce
sce.tumor <- readRDS(snakemake@input[["sce_tumor"]])

print(paste0("Tumor + Normal sce successfully loaded for ", snakemake@wildcards[["cohort"]], " cohort: "))
table(sce.tumor$cell_type)

sce.tumor <- sce.tumor[, sce.tumor$cell_type == paste0("Tumour_", snakemake@wildcards[["celltype"]])]
print("Number of cells in tumor sce: ")
ncol(sce.tumor)

# load infercnv results
samples <- list.files(snakemake@params[["infercnv_output_dir"]])
samples <- samples[!grepl("sample-list", samples)]
infercnv_obj_list <- lapply(samples, function(sample) {
  readRDS(paste0(snakemake@params[["infercnv_output_dir"]], sample, "/run.final.infercnv_obj"))
})
names(infercnv_obj_list) <- samples

# check what's in the infercnv object
print("The infercnv_object has the following attributes, we are insterested in the processed expression data in 'expr.data'")
names(attributes(infercnv_obj_list[[1]]))

# get tumour sample cell indices
print("Finding the indices of all non-reference cells")
all_tumour_indices_list <- lapply(infercnv_obj_list, function(infercnv_obj) {
  all_indices <- infercnv_obj@observation_grouped_cell_indices
  all_indices[[paste0("Normal_", snakemake@wildcards[["celltype"]])]] <- NULL
  all_tumour_indices <- unlist(unname(all_indices))
  all_tumour_indices
})
Reduce(sum, lapply(all_tumour_indices_list, length))

# get variance of smoothed expression of genes in tumour sample cells
print("Creating dataframe with expression variance for each cell")
infer_cnv_frame_list <- lapply(samples, function(sample) {
  infercnv_obj <- infercnv_obj_list[[sample]]
  all_tumour_indices <- all_tumour_indices_list[[sample]]
  infer_cnv_frame <- infercnv_obj@expr.data[,all_tumour_indices]
  infer_cnv_frame
})
names(infer_cnv_frame_list) <- samples

vars_df_list <- lapply(infer_cnv_frame_list, function(infer_cnv_frame) {
  vars <-colVars(infer_cnv_frame)
  vars_df <- data.frame(colnames(infer_cnv_frame), vars)
  rownames(vars_df) <- vars_df$colnames.infer_cnv_frame.
  vars_df$colnames.infer_cnv_frame. <- NULL
  vars_df
})

sce.tumor.infercnv.list <- lapply(infer_cnv_frame_list, function(infer_cnv_frame) {
  sce <- SingleCellExperiment(assays = SimpleList(logcounts = infer_cnv_frame))
  sce$sample <- str_split(colnames(infer_cnv_frame), "-", simplify = T)[,1]
  sce
})

genes.common <- Reduce(intersect, lapply(sce.tumor.infercnv.list, rownames))
sce.tumor.infercnv.list <- lapply(sce.tumor.infercnv.list, function(sce.tumor.infercnv) {
  sce.tumor.infercnv[genes.common,]
})

sce.tumor.infercnv <- Reduce(cbind, sce.tumor.infercnv.list)

# put expression variance into original tumour SCE
print("Taking those cells found in our orignial dataset so we can compare")
sce.tumor.results <- sce.tumor[,Reduce(c, lapply(vars_df_list, rownames))]
colData(sce.tumor.results)$infer_cnv_var <- Reduce(c, lapply(vars_df_list, function(vars_df) {vars_df$vars}))

sce.tumor.results.list <- lapply(vars_df_list, function(vars_df) {
  sce.tumor.results <- sce.tumor[,rownames(vars_df)]
  colData(sce.tumor.results)$infer_cnv_var <- vars_df$vars
  sce.tumor.results
})
names(sce.tumor.results.list) <- samples

# run dimensionality reduction
print("Running dimensionality reduction on original tumor sample expression")
sce.tumor.results <- scater::runPCA(sce.tumor.results)
sce.tumor.results <- scater::runUMAP(sce.tumor.results, n_neighbors = min(30, min(table(sce.tumor.results$sample)) - 1))

print("Running dimensionality reduction on infercnv smoothed tumor sample expression")
sce.tumor.infercnv <- scater::runPCA(sce.tumor.infercnv, subset_row = getTopHVGs(modelGeneVar(sce.tumor.infercnv), prop = 0.1))
sce.tumor.infercnv <- scater::runUMAP(sce.tumor.infercnv, n_neighbors = min(30, min(table(sce.tumor.infercnv$sample)) - 1))

print("Running dimensionality reduction on infercnv smoothed tumor sample expression for each sample")
sce.tumor.infercnv.list <- lapply(sce.tumor.infercnv.list, function (sce.tumor.infercnv) {
  sce.tumor.infercnv <- scater::runPCA(sce.tumor.infercnv, subset_row = getTopHVGs(modelGeneVar(sce.tumor.infercnv), prop = 0.1))
  sce.tumor.infercnv <- scater::runUMAP(sce.tumor.infercnv, n_neighbors = min(30, ncol(sce.tumor.infercnv) - 1))
  sce.tumor.infercnv
})

print("Running dimensionality reduction on original tumor sample expression for each sample")
sce.tumor.results.list <- lapply(sce.tumor.results.list, function (sce.tumor.results) {
  sce.tumor.results <- scater::runPCA(sce.tumor.results)
  sce.tumor.results <- scater::runUMAP(sce.tumor.results, n_neighbors = min(30, ncol(sce.tumor.results) - 1))
  sce.tumor.results
})

# cluster cells
print("Clustering cells using original tumor sample expression")
sce.tumor.results.list <- lapply(sce.tumor.results.list, function (sce.tumor.results) {
  nn.clusters <- clusterCells(sce.tumor.results, use.dimred="PCA", BLUSPARAM=NNGraphParam(k=40))
  colLabels(sce.tumor.results) <- nn.clusters
  sce.tumor.results
})

# compare expression variance in different clusters
print("Computing mean inferred CNV variance in different clusters")
cluster_scores_list <- lapply(sce.tumor.results.list, function (sce.tumor.results) {
  cluster_scores <- aggregate(colData(sce.tumor.results)$infer_cnv_var, list(colData(sce.tumor.results)$label), FUN=mean)
  cluster_scores[order(cluster_scores$x),][1:10,"x", drop=FALSE]
})

# load other infercnv results
cnv.regions.list <- lapply(samples, function(sample) {
  cnv.regions <- read.table(paste0(snakemake@params[['infercnv_output_dir']], sample, '/', 
                                   snakemake@params[['regional_cnv_calls']]),
                            header = T, sep = "\t")
})
names(cnv.regions.list) <- samples
print("Loaded regional CNV calls")

cnv.genes.list <- lapply(samples, function(sample) {
  cnv.genes <- read.table(paste0(snakemake@params[['infercnv_output_dir']], sample, '/',
                                 snakemake@params[['gene_cnv_calls']]),
                            header = T, sep = "\t")
})
names(cnv.genes.list) <- samples
print("Loaded gene level CNV calls")

cell.groupings.list <- lapply(samples, function(sample) {
  cell.groupings <- read.table(paste0(snakemake@params[['infercnv_output_dir']], sample, '/',
                                      snakemake@params[['cell_groupings']]),
                            header = T, sep = "\t")
})
names(cell.groupings.list) <- samples
print("Loaded cell groupings")

sce.tumor.results.list <- lapply(sce.tumor.results.list, function(sce.tumor.results) {
        cell.groupings <- cell.groupings.list[[(sce.tumor.results$sample)[1]]]
        sce.tumor.results$cell_groupings <- plyr::mapvalues(colnames(sce.tumor.results),
                from = cell.groupings$cell,
                to = cell.groupings$cell_group_name,
                warn_missing = F
        )
        sce.tumor.results
})

# cluster cells using infercnv smoothed expression
print("Clustering cells using infercnv smoothed tumor sample expression")
sce.tumor.infercnv.list <- lapply(samples, function(sample) {
  sce.tumor.infercnv <- sce.tumor.infercnv.list[[sample]]
  sce.tumor.results <- sce.tumor.results.list[[sample]]
  
  colData(sce.tumor.infercnv) <- colData(sce.tumor.results)
  
  nn.clusters <- clusterCells(sce.tumor.infercnv, use.dimred="PCA", BLUSPARAM=NNGraphParam(k=40))
  colLabels(sce.tumor.infercnv) <- nn.clusters
  sce.tumor.infercnv
})
names(sce.tumor.infercnv.list) <- samples

# get CNV state variance in each cell based on regional CNV states
print("Computing CNV state variance in each cell based on regional CNV states")
cnv.regions.var.list <- lapply(cnv.regions.list, function(cnv.regions) {
  cnv.regions.var <- cnv.regions %>%
    group_by(cell_group_name) %>%
    summarise(state_variance = var(state, na.rm = T)) %>%
    ungroup()
  cnv.regions.var
})

sce.tumor.infercnv.list <- lapply(sce.tumor.infercnv.list, function(sce.tumor.infercnv) {
  sce.tumor.infercnv$cnv_state_var <- plyr::mapvalues(sce.tumor.infercnv$cell_groupings,
                                                      from = cnv.regions.var.list[[(sce.tumor.infercnv$sample)[1]]]$cell_group_name,
                                                      to = cnv.regions.var.list[[(sce.tumor.infercnv$sample)[1]]]$state_variance) %>% as.numeric()
  sce.tumor.infercnv
})

# label cells as malignant/normal based on CNV state variance
print("Labeling cells as malignant/normal based on CNV state variance")
cnv.state.var.thres <- snakemake@params[['cnv_state_var_thres']]

sce.tumor.infercnv.list <- lapply(sce.tumor.infercnv.list, function(sce.tumor.infercnv) {
  cluster.mean.cnv.state.var <- data.frame(label = sce.tumor.infercnv$label,
                                           cnv.state.var = sce.tumor.infercnv$cnv_state_var)
  cluster.mean.cnv.state.var <- cluster.mean.cnv.state.var %>%
    group_by(label) %>%
    summarise(mean.cnv.state.var = mean(cnv.state.var, na.rm = T)) %>%
    ungroup()
  
  print(cluster.mean.cnv.state.var)
  
  cluster.mean.cnv.state.var$cell_type <- ifelse(cluster.mean.cnv.state.var$mean.cnv.state.var > cnv.state.var.thres, 
                                                 "epithelial malignant", "epithelial normal")
  
  sce.tumor.infercnv$cell_type <- plyr::mapvalues(sce.tumor.infercnv$label,
                                                  from = cluster.mean.cnv.state.var$label,
                                                  to = cluster.mean.cnv.state.var$cell_type)
  
  sce.tumor.infercnv$cell_type_cell_by_cell <- ifelse(sce.tumor.infercnv$cnv_state_var > cnv.state.var.thres, 
                                                      "epithelial malignant", "epithelial normal")
  
  sce.tumor.infercnv
})

# save cell labels
print("Saving cell labels")
labels <- lapply(sce.tumor.infercnv.list, function(sce.tumor.infercnv) {
        data.frame(
                cell_id = colnames(sce.tumor.infercnv),
                label = sce.tumor.infercnv$cell_type
        )
})
labels <- Reduce(rbind, labels)
write_csv(labels, snakemake@output[["cell_labels"]])

# save sces
print("Saving sces")
saveRDS(sce.tumor.results.list, snakemake@output[["sce_tumor_results_list"]])
saveRDS(sce.tumor.infercnv.list, snakemake@output[["sce_tumor_infercnv_list"]])