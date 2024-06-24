suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(stringr)
  library(scales)
  library(BiocParallel)
  library(sjstats)
  library(dittoSeq)
  library(ggplotify)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(cowplot)
  library(circlize)
  library(ComplexHeatmap)
})

cohort.pal <- readRDS(snakemake@input[["cohort_pal"]])

cell_type_rename <- read_csv(snakemake@input[['cell_type_rename']])

# load signature interpretation
sig.interpt <- readxl::read_xlsx(snakemake@input[['sig_interpretation']])

# signatures to remove
ambient.sigs <- read_csv(snakemake@input[['ambient_sigs']])$signature

# load signature loading variance dataframes
sig.loading.var.dis.val <- read_tsv(snakemake@input[['sig_loading_var_dis_val']]) |> column_to_rownames("signature")

sig.loading.var.to.plot <- sig.loading.var.dis.val[!(rownames(sig.loading.var.dis.val) %in% ambient.sigs),]


names(sig.loading.var.to.plot) <- c(
        "Intra-patient heterogeneity\nDiscovery",
        "Inter-patient heterogeneity (mean)\nDiscovery",
        "Inter-patient heterogeneity (median)\nDiscovery",
        "Intra-patient heterogeneity\nValidation",
        "Inter-patient heterogeneity (mean)\nValidation",
        "Inter-patient heterogeneity (median)\nValidation")

print(head(sig.loading.var.to.plot))

sig.loading.var.to.plot <- sig.loading.var.to.plot |>
  rownames_to_column("signature") |>
  pivot_longer(cols = contains("heterogeneity"), values_to = "heterogeneity", names_to = "measure") |>
  separate_wider_delim(measure, "\n", names = c("measure", "group")) |>
  dplyr::filter(str_detect(measure, "mean", negate = TRUE)) |>
  pivot_wider(names_from = "group", values_from = "heterogeneity") |>
  mutate(celltype = str_split(signature, " Sig ", simplify = TRUE)[,1])

sig.loading.var.to.plot$celltype <-
  plyr::mapvalues(sig.loading.var.to.plot$celltype,
                  from = cell_type_rename$old_name,
                  to = cell_type_rename$new_name)


# panel A
p.sig.hetero.agree <- ggscatter(sig.loading.var.to.plot, 
          x = "Discovery", y = "Validation",
          xlab = "Gene Program loading variance (Discovery)", ylab = "Gene Program\nloading variance (Validation)",
          color = "measure", palette = "jco",
          add = "reg.line"
) +
  facet_wrap(~ measure, scales = "free") +
  stat_cor()
ggsave(snakemake@output[['figure3_a']], width = snakemake@params[['figure3_a_width']], height = snakemake@params[['figure3_a_height']], units = "in", dpi = 360)

print("Signature loading heterogeneity plot (combined) successfully created")

# panel B
p.sig.hetero.agree.ct <- ggplot(sig.loading.var.to.plot, aes(x = Discovery, y = Validation, color = measure)) +
  geom_point() + 
  stat_cor(method = "spearman") +
  geom_smooth(method = "lm") +
  facet_wrap(~ celltype, scales = "free", ncol = 4) +
  scale_color_jco() +
  labs(x = "Gene Program loading variance (Discovery)", y = "Gene Program loading variance (Validation)") +
  theme_pubr() + 
  guides(color = guide_legend(title = "Measure", override.aes = aes(label = "")))
ggsave(snakemake@output[['figure3_b']], width = snakemake@params[['figure3_b_width']], height = snakemake@params[['figure3_b_height']], units = "in", dpi = 360)

print("Signature loading heterogeneity plot (per cell type) successfully created")

# plot Figure 2
design <- "
AAA
BBB
BBB
"

p.sig.hetero.agree <- p.sig.hetero.agree + theme(legend.position = "none", aspect.ratio = 0.8, panel.spacing = unit(5, "lines"))

pdf(file = snakemake@output[["figure3_pdf"]], width = snakemake@params[["figure3_width"]], height = snakemake@params[["figure3_height"]])
p.sig.hetero.agree + p.sig.hetero.agree.ct + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))
dev.off()

print("Figure 3 PDF successfully created")

png(file = snakemake@output[["figure3_png"]], width = snakemake@params[["figure3_width"]], height = snakemake@params[["figure3_height"]], units = "in", res = 360)
p.sig.hetero.agree + p.sig.hetero.agree.ct + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))
dev.off()

print("Figure 3 PNG successfully created")