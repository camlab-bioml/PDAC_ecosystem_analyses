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

cell_type_rename <- read_csv(snakemake@input[["cell_type_rename"]])

#cell_type_rename$old_name <- gsub(" ", "_", cell_type_rename$old_name)
#cell_type_rename$old_name <- gsub("-", "_", cell_type_rename$old_name)
#cell_type_rename$old_name <- gsub(",", "", cell_type_rename$old_name)

print(cell_type_rename)

# load signature interpretation
sig.interpt <- readxl::read_xlsx(snakemake@input[["sig_interpretation"]])

#sig.interpt$signature <- gsub(" ", "_", sig.interpt$signature)
#sig.interpt$signature <- gsub("-", "_", sig.interpt$signature)
#sig.interpt$signature <- gsub(",", "", sig.interpt$signature)

# signatures to remove
ambient.sigs <- read_csv(snakemake@input[["ambient_sigs"]])$signature

#ambient.sigs <- gsub(" ", "_", ambient.sigs)
#ambient.sigs <- gsub("-", "_", ambient.sigs)
#ambient.sigs <- gsub(",", "", ambient.sigs)

# load cohort color palette
cohort_pal <- readRDS(snakemake@input[["cohort_pal"]])

# load celltype color palette
celltype_pal <- readRDS(snakemake@input[["celltype_pal"]])

celltype_pal$Cell_type <- str_split(celltype_pal$Cell_type_dis, " \\(", simplify = TRUE)[, 1]

celltype_pal_to_use <- celltype_pal$color_dis
names(celltype_pal_to_use) <- celltype_pal$Cell_type

# load stan data and other parameters
celltypes <- snakemake@params[["celltypes"]]
number.of.niches <- snakemake@params[["number_of_niches"]]
nIter <- snakemake@params[["nIter"]]

celltypes <- gsub(" ", "_", celltypes)
celltypes <- gsub("-", "_", celltypes)
celltypes <- gsub(",", "", celltypes)

# load data
niches.init.value.collapsed <- readRDS(snakemake@input[["microenvironment_niche_factors_init_value_collapsed"]])
niches.init.value.collapsed.validation <- readRDS(snakemake@input[["microenvironment_niche_factors_init_value_scored_val"]])
niche.loadings.init.value.collapsed <- readRDS(snakemake@input[["niche_factor_loadings_init_value_collapsed"]])
niche.loadings.init.value.collapsed.validation <- readRDS(snakemake@input[["niche_factor_loadings_init_value_scored_val"]])

niches.collapsed <- readRDS(snakemake@input[["microenvironment_niche_factors_collapsed"]])
niches.collapsed.validation <- readRDS(snakemake@input[["microenvironment_niche_factors_scored_val"]])
cov.i.collapsed <- readRDS(snakemake@input[["intrinsic_covariance_matrices_collapsed"]])
cov.i.collapsed.validation <- readRDS(snakemake@input[["intrinsic_covariance_matrices_scored_val"]])

# add some dis/val information
rownames(niches.init.value.collapsed) <- paste("Discovery", rownames(niches.init.value.collapsed))
rownames(niches.init.value.collapsed.validation) <- paste("Validation", rownames(niches.init.value.collapsed.validation))

rownames(niches.collapsed) <- paste("Discovery", rownames(niches.collapsed))
rownames(niches.collapsed.validation) <- paste("Validation", rownames(niches.collapsed.validation))

# combine data
niches <- cbind(t(niches.collapsed), t(niches.collapsed.validation))
niches.init.value <- rbind(niches.init.value.collapsed, niches.init.value.collapsed.validation) |> t()

cov.i.list <- lapply(celltypes, function(ct) {
  print(ct)
  print(dim(cov.i.collapsed[[paste0("cov_i_", ct)]]))
  print(dim(cov.i.collapsed.validation[[paste0("cov_i_", ct)]]))
  cov.i.collapsed[[paste0("cov_i_", ct)]] - cov.i.collapsed.validation[[paste0("cov_i_", ct)]]
})
names(cov.i.list) <- celltypes

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



# draw niches
head(niches)
niches.to.plot <- t(scale(niches, center = TRUE, scale = TRUE))

row.split <- str_split(rownames(niches.to.plot), " ", simplify = TRUE)[,1]
rownames(niches.to.plot) <- gsub("Discovery |Validation ", "", rownames(niches.to.plot))

column.split <- plyr::mapvalues(gsub(" [0-9]$| [0-9][0-9]$", "", colnames(niches.to.plot)),
  from = cell_type_rename$old_name,
  to = cell_type_rename$new_name
)
colnames(niches.to.plot) <- plyr::mapvalues(colnames(niches.to.plot),
  from = sig.interpt$signature,
  to = sig.interpt$`short interpretation`
)

col_fun <- colorRamp2(seq(min(niches.to.plot), max(niches.to.plot), length.out = 100), viridisLite::viridis(100, option = "C"))

png(snakemake@output[["microenvironment_niche_factors_plot"]], width = 17, height = 14, units = "in", res = 321)
Heatmap(niches.to.plot,
  height = nrow(niches.to.plot) * unit(0.4, "in"),
  name = "Signature\nloading",
  top_annotation = columnAnnotation(
    Celltype = column.split,
    col = list(Celltype = celltype_pal_to_use),
    show_annotation_name = TRUE
  ),
  left_annotation = rowAnnotation(
    Group = row.split,
    show_annotation_name = TRUE
  ),
  column_names_rot = 90,
  col = col_fun
) |>
  draw(column_title = paste0("Microenvironment Niche Factors", " - ", number.of.niches, " niches, ", nIter, " iterations"))
dev.off()

# draw niches init value
head(niches.init.value)
niches.init.value.to.plot <- t(niches.init.value)

row.split <- str_split(rownames(niches.init.value.to.plot), " ", simplify = TRUE)[,1]
rownames(niches.init.value.to.plot) <- gsub("Discovery |Validation ", "", rownames(niches.init.value.to.plot))

column.split <- plyr::mapvalues(gsub(" [0-9]$| [0-9][0-9]$", "", colnames(niches.init.value.to.plot)),
  from = cell_type_rename$old_name,
  to = cell_type_rename$new_name
)
colnames(niches.init.value.to.plot) <- plyr::mapvalues(colnames(niches.init.value.to.plot),
  from = sig.interpt$signature,
  to = sig.interpt$`short interpretation`
)

png(snakemake@output[["microenvironment_niche_factors_init_value_plot"]], width = 17, height = 14, units = "in", res = 321)
Heatmap(niches.init.value.to.plot,
  height = nrow(niches.init.value.to.plot) * unit(0.4, "in"),
  name = "Signature\nloading",
  top_annotation = columnAnnotation(
    celltype = column.split,
    col = list(celltype = celltype_pal_to_use),
    show_annotation_name = FALSE
  ),
  row_split = row.split,
  column_names_rot = 90,
  col = viridisLite::viridis(100, option = "C")
) |>
  draw(column_title = paste0("Initialized Microenvironment Niche Factors", " - ", number.of.niches, " niches, ", nIter, " iterations"))

# draw correlation plots
png(snakemake@output[["microenvironment_niche_factors_correlation_plot"]], width = 12, height = 10, units = "in", res = 321)
chart.Correlation(niches, histogram=TRUE, pch=19)
dev.off()

png(snakemake@output[["microenvironment_niche_factors_init_value_correlation_plot"]], width = 12, height = 10, units = "in", res = 321)
chart.Correlation(niches.init.value, histogram=TRUE, pch=19)
dev.off()

cov.i.ht.list <- lapply(names(cov.i.list), function(ct) {
  cov.i <- cov.i.list[[ct]]
  cov.i.ht <- Heatmap(cov.i, name = "Collapsed - Scored validation\nintrinsic covariance", 
                      show_row_names = FALSE, 
                      #show_column_names = FALSE, 
                      #cluster_rows = FALSE, cluster_columns = FALSE, 
                      #column_title = "cov(i)", column_title_gp = gpar(fontsize = 20), 
                      #row_title = "cov(i)", row_title_gp = gpar(fontsize = 20), 
                      #heatmap_legend_param = list(title = "cov(i)", title_gp = gpar(fontsize = 20)), 
                      col = circlize::colorRamp2(c(min(cov.i), 0, max(cov.i)), c("blue", "white", "red"))
                      )
  cov.i.ht
})

cov.i.ht.list <- lapply(cov.i.ht.list, function(cov.i.ht) draw(cov.i.ht) |> grid.grabExpr())

png(snakemake@output[["intrinsic_covariance_matrices_correlation_plot"]], width = 20, height = 10, units = "in", res = 321)
#ComplexHeatmap::draw(ht.list, merge_legends = TRUE)
#plot_grid(cov.i.ht.list, nrow = 2)
wrap_plots(cov.i.ht.list, nrow = 2)
dev.off()


# draw correlation bar plot
cell_type_rename$old_name <- gsub(" ", "_", cell_type_rename$old_name)
cell_type_rename$old_name <- gsub("-", "_", cell_type_rename$old_name)
cell_type_rename$old_name <- gsub(",", "", cell_type_rename$old_name)

cov.i.cor.test$cell_type <- plyr::mapvalues(cov.i.cor.test$cell_type,
  from = cell_type_rename$old_name,
  to = cell_type_rename$new_name
)

p.cov.corr.bar <- ggplot(cov.i.cor.test, aes(x = cell_type, y = corr, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_pal_to_use) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #geom_errorbar(aes(ymin = corr - 1.96 * sqrt((1 - corr^2) / num_sig), ymax = corr + 1.96 * sqrt((1 - corr^2) / num_sig)), width = 0.2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  labs(x = "Cell type", y = "Pearson correlation\ndiscovery intrinsic cov. and validation intrinsic cov.") +
  geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave(snakemake@output[["intrinsic_covariance_matrices_correlation_bar_plot"]], plot = p.cov.corr.bar, width = 7, height = 9, units = "in", dpi = 321)