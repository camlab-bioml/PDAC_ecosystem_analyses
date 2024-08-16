suppressPackageStartupMessages({
        library(tidyverse)
        library(ComplexHeatmap)
	library(ggplot2)
	library(ggpubr)
	library(ggsci)
        library(ggrepel)
	library(viridisLite)
	library(circlize)
        library(patchwork)
        library(magick)
})

cell_type_rename <- read_csv(snakemake@input[["cell_type_rename"]])

sig.interpt <- readxl::read_xlsx(snakemake@input[["sig_interpretation"]])

# signatures to remove
ambient.sigs <- read_csv(snakemake@input[["ambient_sigs"]])$signature

ambient.sigs <- gsub(" Sig ", " ", ambient.sigs)

# load cohort color palette
cohort_pal <- readRDS(snakemake@input[["cohort_pal"]])

# load celltype color palette
celltype_pal <- readRDS(snakemake@input[["celltype_pal"]])

celltype_pal$Cell_type <- str_split(celltype_pal$Cell_type_dis, " \\(", simplify = TRUE)[, 1]

celltype_pal_to_use <- celltype_pal$color_dis
names(celltype_pal_to_use) <- celltype_pal$Cell_type

# load scematic
schematic <- image_read_pdf(snakemake@input[['schematic']]) |> image_ggplot()

# panel A
# print the schematic
ggsave(filename = snakemake@output[["figure4_a"]], plot = schematic, dpi = 360)
print("Schematic plot successfully created: ")
print(snakemake@output[["figure4_a"]])

# Load the data
dfc <- read_tsv(snakemake@input[["patient_profiles_correlation_data_frame_full_and_intra"]])

# panel B
# Are signature correlations consistent between cohorts?
## remove ambient signatures
dfc.for.plot <- dfc
dfc.for.plot <- dfc.for.plot |>
        filter(!(signature_1 %in% ambient.sigs)) |>
        filter(!(signature_2 %in% ambient.sigs))

## update cell type labels for plotting
dfc.for.plot$cell_type_1 <- plyr::mapvalues(dfc.for.plot$cell_type_1,
                                            from = cell_type_rename$old_name,
                                            to = cell_type_rename$new_name)
dfc.for.plot$cell_type_2 <- plyr::mapvalues(dfc.for.plot$cell_type_2,
                                            from = cell_type_rename$old_name,
                                            to = cell_type_rename$new_name)


## plot overall correlations
p.corr.overall <- ggplot(dfc.for.plot, aes(x = correlation_discovery, y = correlation_validation)) +
        geom_point(aes(color = gsub("Co-occurrence ", "Co-occurrence\n", same_cell_type_for_plot)), alpha = 0.5) +
        # facet_wrap(~ cell_type_1, scales = "free", nrow = 2) +
        geom_smooth(method = "lm", colour = "grey30") +
        stat_cor(method = "spearman", cor.coef.name = "rho") +
        labs(
                color = "Same celltype",
                x = "Co-occurrence (Discovery)",
                y = "Co-occurrence (Validation)"
        ) +
        theme_pubr() +
        theme(
                legend.title = element_blank(),
                legend.position = c(0.8, 0.2),
                legend.text = element_text(size = 14),
                axis.title.x = element_text(face = "bold", size = 16),
                axis.title.y = element_text(face = "bold", size = 16),
        ) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
        theme(plot.margin = unit(c(0,1,0,1),"cm"))
ggsave(snakemake@output[["figure4_b"]], width = snakemake@params[["figure4_b_width"]], height = snakemake@params[["figure4_b_height"]], units = "in", dpi = 360)

print("Figure 4B successfully created")

# panel C
# plot intra-cell type correlations
p.corr.intra <- filter(dfc.for.plot, same_cell_type) |>
        ggplot(aes(x = correlation_discovery, y = correlation_validation)) +
        geom_point(aes(colour = cell_type_1)) +
        scale_color_manual(values = celltype_pal_to_use) +
        facet_wrap(~cell_type_1, scales = "free", nrow = 2) +
        geom_smooth(method = "lm", level = 0.9) +
        stat_cor(method = "spearman", cor.coef.name = "rho", 
                 aes(label = paste(..r.label.., cut(..p.., 
                                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                                labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")),
                               sep = "~")),
                 label.x.npc = "left", label.y.npc = "top") +
        labs(
                colour = "Celltype",
                x = "Program co-occurrence in Discovery",
                y = "Program co-occurrence\nin Validation"
        ) +
        theme_pubr() +
        theme(axis.title = element_text(face = "bold"))
ggsave(snakemake@output[["figure4_c"]], width = snakemake@params[["figure4_c_width"]], height = snakemake@params[["figure4_c_height"]], units = "in", dpi = 360)

print("Figure 4C successfully created")

# panel D
# plot inter-cell type correlations
dfc.for.plot_list <- list()

for (ct in union(unique(dfc.for.plot$cell_type_1), unique(dfc.for.plot$cell_type_2))) {
        dfc.for.plot_list[[ct]] <- dfc.for.plot %>%
                filter(!same_cell_type) %>%
                filter(cell_type_1 == ct | cell_type_2 == ct) |>
                mutate(facet_cell_type = ct) |>
                mutate(cell_type = ifelse(cell_type_1 == ct, cell_type_2, cell_type_1))
}
rm(ct)

p.corr.inter <- Reduce(rbind, dfc.for.plot_list) |>
        ggplot(aes(x = correlation_discovery, y = correlation_validation)) +
        geom_point(aes(colour = cell_type)) +
        scale_color_manual(values = celltype_pal_to_use) +
        facet_wrap(~facet_cell_type, scales = "free", nrow = 2) +
        geom_smooth(method = "lm") +
        stat_cor(method = "spearman", cor.coef.name = "rho", 
                 aes(label = paste(..r.label.., cut(..p.., 
                                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                                labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")),
                               sep = "~")),
                 label.x.npc = "left", label.y.npc = "top") +
        labs(
                colour = "Celltype",
                x = "Program co-occurrence in Discovery",
                y = "Program co-occurrence in Validation"
        ) +
        theme_pubr() +
        theme(axis.title = element_text(face = "bold"))
ggsave(snakemake@output[["figure4_d"]], width = snakemake@params[["figure4_d_width"]], height = snakemake@params[["figure4_d_height"]], units = "in", dpi = 360)

print("Figure 4D successfully created")

# panel E
# plot general signature co-occurrence agreement between discovery an validation
dis.val.agree <- read_tsv(snakemake@input[["signature_cooccurrence_agreement_data_frame"]])

dis.val.agree <- dis.val.agree |>
  mutate(celltype = plyr::mapvalues(celltype, from = cell_type_rename$old_name, to = cell_type_rename$new_name)) |>
  mutate(signature = plyr::mapvalues(signature, from = sig.interpt$signature, to = sig.interpt$`short interpretation`))

dis.val.agree <- dis.val.agree |> 
  filter(!grepl("Ambient RNA|^MALAT1/NEAT1$", signature)) |>
  group_by(signature) |>
  mutate(unique_signature = paste0(signature, "_-_", row_number())) |>
  ungroup() |>
  arrange(desc(dis_val_corr))

margin_spacer <- function(x) {
  # where x is the column in your dataset
  #left_length <- nchar(levels(factor(x)))[1]
  left_length <- nchar(x)[1]
  if (left_length > 8) {
    return((left_length - 8) * 4)
  }
  else
    return(0)
}

p.cooccur.agree <- ggplot(dis.val.agree, aes(x = reorder(unique_signature, -dis_val_corr), y = dis_val_corr)) +
  geom_bar(stat = "identity", aes(fill = celltype)) +
  ylim(0, NA) +
  scale_x_discrete(labels = dis.val.agree$signature) + 
  scale_fill_manual(values = celltype_pal_to_use) +
  labs(x = "Cell state programs", y = "Correlation of co-occurrence\nbetween Discovery and Validation", fill = "Cell type") +
  geom_text(aes(label = ifelse(p_value < 0.01, "*", "")), 
            position = position_dodge(width = .9), vjust = .3, size = 20 / .pt) +
  theme_pubr() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.ticks.x = element_blank(),
        #plot.margin = margin(l = 0 + margin_spacer(dis.val.agree$signature))
        #plot.margin = margin(l = 0 + 32)
        )
ggsave(snakemake@output[["figure4_e"]], width = snakemake@params[["figure4_e_width"]], height = snakemake@params[["figure4_e_height"]], units = "in", dpi = 360)

print("Figure 4E successfully created")

# panel F
# plot signature co-occurrence in discovery vs. validation for each signature
dfc.list <- readRDS(snakemake@input[["signature_correlation_comparison_data_frame_list_examples"]])

## make sure some signatures with same inpterpretation don't get stacked on each other
dfc.list <- lapply(dfc.list, function(dfc.sig) {
  dfc.sig$celltype <- plyr::mapvalues(dfc.sig$celltype, from = cell_type_rename$old_name, cell_type_rename$new_name)
  dfc.sig$signature <- plyr::mapvalues(dfc.sig$signature, 
                                         from = sig.interpt$signature, 
                                         to = sig.interpt$`short interpretation`,
                                       warn_missing = FALSE)
  dfc.sig |> 
    filter(!grepl("Ambient RNA|^MALAT1/NEAT1$", signature)) |>
    group_by(signature) |>
    mutate(unique_signature = paste0(signature, "_-_", row_number())) |>
    ungroup() |>
    arrange(desc(correlation_discovery))
})

## draw scatter plots
p.cooccur.scatter.list <- lapply(names(dfc.list), function(sig) {
  if(nrow(dfc.list[[sig]]) == 0) {
    return()
  }
  dfc.sig <- dfc.list[[sig]] |>
    mutate(sig.of.interest = ifelse(((abs(correlation_discovery) > 0.4 | abs(correlation_validation) > 0.35) & 
                                       abs(correlation_discovery - correlation_validation) <= 0.26), signature, NA))
  plot.title <- sig.interpt |> filter(signature == sig) |> pull(`short interpretation`)
  plot.title <- paste0("Correlation with - ", plot.title)
#   plot.title <- paste0(plyr::mapvalues(gsub(" [0-9]$| [0-9][0-9]$", "", sig), 
#                                        from = cell_type_rename$old_name, 
#                                        to = cell_type_rename$new_name, warn_missing = FALSE), ":\n", plot.title)

  
  ggplot(dfc.sig, aes(x = correlation_discovery, y = correlation_validation, label = sig.of.interest, fill = celltype, color = celltype)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_point(aes(fill = celltype, color = celltype), size = 4) +
    geom_label_repel(show.legend = F, max.overlaps = 20, alpha = 1, colour = "white", 
                     arrow = arrow(type = "closed", angle = 20, length = unit(0.08, "inches"))) +
    scale_fill_manual(values = celltype_pal_to_use) + 
    scale_color_manual(values = celltype_pal_to_use) +
    labs(title = plot.title,
         x = "Correlation in Discovery",
         y = "Correlation in Validation",
         color = "Cell type") + 
    theme_pubr() +
    theme(axis.title = element_text(size = 12, face = "bold"),
          title = element_text(face = "bold")) + 
    guides(fill = "none") + 
    theme(legend.position = "none")
  
})
p.cooccur.scatter <- ggarrange(plotlist = p.cooccur.scatter.list, nrow = 2, ncol = 2, common.legend = TRUE, legend = "none")
ggsave(filename = snakemake@output[["figure4_f"]], plot = p.cooccur.scatter, width = snakemake@params[["figure4_f_width"]], height = snakemake@params[["figure4_f_height"]], units = "in", dpi = 360)

print("Figure 4F successfully created")

# plot Figure 4
design <- "
AAABBBCCCCC
AAABBBCCCCC
AAABBBCCCCC
EEEFFFDDDDD
EEEFFFDDDDD
EEEFFFDDDDD
GGGHHHIIIII
GGGHHHIIIII
GGGHHHIIIII
"

p.corr.intra <- p.corr.intra + theme(legend.position = "none", axis.title.x = element_blank())
p.corr.inter <- p.corr.inter + theme(legend.position = "none")
p.cooccur.agree <- p.cooccur.agree + theme(legend.position = c(0.5, 1), legend.direction = "horizontal")

pdf(file = snakemake@output[["figure4_pdf"]], width = snakemake@params[["figure4_width"]], height = snakemake@params[["figure4_height"]])
schematic + p.corr.overall + p.corr.intra + p.corr.inter + p.cooccur.scatter.list[[1]] + p.cooccur.scatter.list[[2]] + p.cooccur.scatter.list[[3]] + p.cooccur.scatter.list[[4]] + p.cooccur.agree + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 18))
dev.off()

print("Figure 4 PDF successfully created")

png(file = snakemake@output[["figure4_png"]], width = snakemake@params[["figure4_width"]], height = snakemake@params[["figure4_height"]], units = "in", res = 360)
schematic + p.corr.overall + p.corr.intra + p.corr.inter + p.cooccur.scatter.list[[1]] + p.cooccur.scatter.list[[2]] + p.cooccur.scatter.list[[3]] + p.cooccur.scatter.list[[4]] + p.cooccur.agree + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 18))
dev.off()

print("Figure 4 PNG successfully created")