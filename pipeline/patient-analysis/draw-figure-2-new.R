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
        library(gridExtra)
        library(Cairo)
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
ggsave(filename = snakemake@output[["figure2_a"]], plot = schematic, dpi = 600)
print("Schematic plot successfully created: ")
print(snakemake@output[["figure2_a"]])

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

dfc.for.plot$same_cell_type_for_plot <- gsub("between", "across", dfc.for.plot$same_cell_type_for_plot)

## plot overall correlations
p.corr.overall <- ggplot(dfc.for.plot, aes(x = correlation_discovery, y = correlation_validation)) +
        geom_point(aes(color = gsub("Co-occurrence ", "Co-occurrence\n", same_cell_type_for_plot)), alpha = 0.5) +
        # facet_wrap(~ cell_type_1, scales = "free", nrow = 2) +
        geom_smooth(method = "lm", colour = "grey30") +
        stat_cor(method = "spearman", cor.coef.name = "rho", label.x = -0.55, label.y = -0.55) +
        labs(
                color = "Same celltype",
                x = "Co-occurrence (Discovery)",
                y = "Co-occurrence (Validation)"
        ) +
        theme_pubr() +
        theme(
                legend.title = element_blank(),
                legend.position = c(0.85, 0.18),
                legend.text = element_text(size = 11),
                legend.background = element_blank(),
                axis.title.x = element_text(face = "bold", size = 14),
                axis.title.y = element_text(face = "bold", size = 14),
        ) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 5))) +
        theme(plot.margin = unit(c(0,1,0,1),"cm"))

# Function to calculate confusion matrix and metrics
calculate_metrics <- function(data) {
  TP <- sum(data$correlation_discovery > 0 & data$correlation_validation > 0)
  FP <- sum(data$correlation_discovery > 0 & data$correlation_validation <= 0)
  FN <- sum(data$correlation_discovery <= 0 & data$correlation_validation > 0)
  TN <- sum(data$correlation_discovery <= 0 & data$correlation_validation <= 0)
  
  c(
    PPV = TP / (TP + FP),       # Positive Predictive Value
    NPV = TN / (TN + FN),       # Negative Predictive Value
    #FPR = FP / (FP + TN),
    #FNR = FN / (FN + TP),
    #TNR = TN / (FP + TN),
    Recall = TP / (TP + FN),
    Precision = TP / (TP + FP)
  )
}

# Calculate metrics for each group
same_type_metrics <- calculate_metrics(subset(dfc.for.plot, same_cell_type == TRUE))
diff_type_metrics <- calculate_metrics(subset(dfc.for.plot, same_cell_type == FALSE))
#overall_metrics <- calculate_metrics(dfc.for.plot)

# Create combined metrics table
metrics <- data.frame(
  Metric = c("PPV", "NPV", "Recall", "Precision"),
  `Same Type` = sprintf("%.1f%%", same_type_metrics * 100),
  `Different Type` = sprintf("%.1f%%", diff_type_metrics * 100)#,
  #`Overall` = sprintf("%.1f%%", overall_metrics * 100)
)

# Create and add table to plot
table_theme <- ttheme_minimal(
  base_size = 8,
  padding = unit(c(1, 2), "mm")
)

metrics_table <- tableGrob(metrics, rows = NULL, theme = table_theme)

p.corr.overall <- p.corr.overall + 
  annotation_custom(
    metrics_table,
    ymin = 0.5, ymax = 1,
    xmin = -0.5, xmax = 0.15
  )

ggsave(snakemake@output[["figure2_b"]], width = snakemake@params[["figure2_b_width"]], height = snakemake@params[["figure2_b_height"]], units = "in", dpi = 600)

print("Figure 2B successfully created")

# add-on panel B (potential new panel B)
# Get all pairs containing each signature and compute FDR
library(dplyr)
library(tidyr)

# Function to compute metrics for one signature
compute_signature_metrics <- function(sig, df) {
        # Get all pairs containing this signature
        sig_pairs <- df %>%
                filter(signature_1 == sig | signature_2 == sig)

        # Count true and false positives
        true_positives <- sum(sig_pairs$correlation_discovery > 0 &
                sig_pairs$correlation_validation > 0)
        false_positives <- sum(sig_pairs$correlation_discovery > 0 &
                sig_pairs$correlation_validation < 0)
        true_negatives <- sum(sig_pairs$correlation_discovery < 0 &
                sig_pairs$correlation_validation < 0)
        false_negatives <- sum(sig_pairs$correlation_discovery < 0 &
                sig_pairs$correlation_validation > 0)

        # Calculate FDR and other metrics
        fdr <- if ((false_positives + true_positives) > 0) {
                false_positives / (false_positives + true_positives)
        } else {
                NA
        }
        frr <- if ((false_negatives + true_positives) > 0) {
                false_negatives / (false_negatives + true_positives)
        } else {
                NA
        }
        sensitivity <- if ((true_positives + false_negatives) > 0) {
                true_positives / (true_positives + false_negatives)
        } else {
                NA
        }
        specificity <- if ((true_negatives + false_positives) > 0) {
                true_negatives / (true_negatives + false_positives)
        } else {
                NA
        }
        ppv <- if ((true_positives + false_positives) > 0) {
                true_positives / (true_positives + false_positives)
        } else {
                NA
        }
        npv <- if ((true_negatives + false_negatives) > 0) {
                true_negatives / (true_negatives + false_negatives)
        } else {
                NA
        }

        # Get cell type and mean discovery correlation
        cell_type <- df %>%
                summarize(
                        type1 = first(cell_type_1[signature_1 == sig]),
                        type2 = first(cell_type_2[signature_2 == sig])
                ) %>%
                transmute(cell_type = coalesce(type1, type2)) %>%
                pull(cell_type)

        mean_discovery_corr <- mean(sig_pairs$correlation_discovery)

        # Return results
        data.frame(
                signature = sig,
                celltype = cell_type,
                true_positives = true_positives,
                false_positives = false_positives,
                true_negatives = true_negatives,
                false_negatives = false_negatives,
                positive_predictive_value = ppv,
                negative_predictive_value = npv,
                false_discovery_rate = fdr,
                false_rejection_rate = frr,
                sensitivity = sensitivity,
                specificity = specificity,
                mean_discovery_correlation = mean_discovery_corr,
                total_pairs = nrow(sig_pairs)
        )
}

# Get unique signatures
all_signatures <- unique(c(dfc.for.plot$signature_1, dfc.for.plot$signature_2))

# Apply function to all signatures
signature_metrics_summary <- do.call(rbind, lapply(all_signatures, function(sig) {
        compute_signature_metrics(sig, dfc.for.plot)
}))
signature_metrics_summary <- signature_metrics_summary |>
        mutate(signature = plyr::mapvalues(signature, from = sig.interpt$signature, to = sig.interpt$`short interpretation`))

# save signature metrics summary
write_tsv(signature_metrics_summary, snakemake@output[["signature_metrics_summary"]])

# Sort by FDR
signature_fdr_summary <- signature_metrics_summary %>%
        arrange(false_discovery_rate)

# View results
print(signature_fdr_summary)

# panel C
# plot intra-cell type correlations
p.corr.intra <- filter(dfc.for.plot, same_cell_type) |>
        ggplot(aes(x = correlation_discovery, y = correlation_validation)) +
        geom_point(aes(colour = cell_type_1), size = 1, alpha = 0.7) +
        scale_color_manual(values = celltype_pal_to_use) +
        facet_wrap(~cell_type_1, scales = "free", nrow = 2) +
        geom_smooth(method = "lm", level = 0.9, color = "grey30", color = "grey70", alpha = 0.2) +
        stat_cor(method = "spearman", cor.coef.name = "rho", 
                 aes(label = paste(..r.label.., cut(..p.., 
                                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                                labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")),
                               sep = "~")),
                 label.x.npc = "left", label.y.npc = "top") +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 4),
                           labels = scales::number_format(accuracy = 0.1)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                           labels = scales::number_format(accuracy = 0.1)) +
        labs(
                colour = "Celltype",
                x = "Program co-occurrence in Discovery",
                y = "Program co-occurrence in Validation"
        ) +
        theme_pubr() +
        theme(axis.title = element_text(face = "bold"),
              strip.text = element_text(face = "bold"),
              strip.background = element_rect(fill = "white", color = "white", size = 0.5),
              panel.grid.major = element_line(color = "grey90", size = 0.5),
              panel.grid.minor = element_blank())
ggsave(snakemake@output[["figure2_c"]], width = snakemake@params[["figure2_c_width"]], height = snakemake@params[["figure2_c_height"]], units = "in", dpi = 600)

print("Figure 2C successfully created")

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
        geom_point(aes(colour = cell_type), size = 1, alpha = 0.7) +
        scale_color_manual(values = celltype_pal_to_use) +
        facet_wrap(~facet_cell_type, scales = "free", nrow = 2) +
        geom_smooth(method = "lm", color = "grey30", fill = "grey70", alpha = 0.2) +
        stat_cor(method = "spearman", cor.coef.name = "rho", 
                 size = 3,
                 aes(label = paste(..r.label.., cut(..p.., 
                                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                                labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")),
                               sep = "~")),
                 label.x.npc = "left", label.y.npc = "top") +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 4),
                           labels = scales::number_format(accuracy = 0.1)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                           labels = scales::number_format(accuracy = 0.1)) +
        labs(
                colour = "Celltype",
                x = "Program co-occurrence in Discovery",
                y = "Program co-occurrence in Validation"
        ) +
        theme_pubr() +
        theme(axis.title = element_text(face = "bold"),
              strip.text = element_text(face = "bold"),
              strip.background = element_rect(fill = "white", color = "white", size = 0.5),
              panel.grid.major = element_line(color = "grey90", size = 0.5),
              panel.grid.minor = element_blank())
ggsave(snakemake@output[["figure2_d"]], width = snakemake@params[["figure2_d_width"]], height = snakemake@params[["figure2_d_height"]], units = "in", dpi = 600)

print("Figure 2D successfully created")

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
  arrange(desc(dis_val_corr)) # sort by correlation

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
  labs(x = "Programs", y = "Co-occurrence reproducibility\nbetween Discovery and Validation", fill = "Cell type") +
  geom_text(aes(label = ifelse(p_value < 0.01, "*", "")), 
            position = position_dodge(width = .9), vjust = 1.2, size = 16 / .pt) +
  theme_pubr() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 13, face = "bold")
        #plot.margin = margin(l = 0 + margin_spacer(dis.val.agree$signature))
        #plot.margin = margin(l = 0 + 32)
        )

# Get max y value from original plot data
y_max <- max(layer_data(p.cooccur.agree)$y)

# Add PPV track - first sort data to match plot order
dis.val.agree <- dis.val.agree |>
        mutate(ct_sig = paste0(celltype, "_-_", signature)) |>
        arrange(desc(dis_val_corr))
signature_fdr_summary <- signature_fdr_summary |>
        mutate(ct_sig = paste0(celltype, "_-_", signature))

# Create PPV data with proper matching
ppv_data <- data.frame(
    ct_sig = dis.val.agree$ct_sig,
    x = seq_len(nrow(dis.val.agree)),
    ppv = signature_fdr_summary$positive_predictive_value[
        match(dis.val.agree$ct_sig, signature_fdr_summary$ct_sig)
    ]
) |>
    # Sort to match main plot order
    arrange(match(ct_sig, dis.val.agree$ct_sig[order(-dis.val.agree$dis_val_corr)]))
print(ppv_data)

# Create final plot with PPV track
p.cooccur.agree.with.ppv <- p.cooccur.agree +
    geom_col(
        data = ppv_data,
        aes(x = x, y = -0.15 * ppv),
        fill = "darkgrey",
        width = 0.8
    ) +
    scale_y_continuous(
        limits = c(-0.15 * 1.0, y_max),
        expand = c(0, 0),
        sec.axis = sec_axis(
            ~./(-0.15),
            name = "PPV",
            labels = scales::percent_format(accuracy = 1),
            breaks = c(0, 0.5, 1.0)
        )
    )
ggsave(snakemake@output[["figure2_e"]], width = snakemake@params[["figure2_e_width"]], height = snakemake@params[["figure2_e_height"]], units = "in", dpi = 600)

print("Figure 2E successfully created")

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
  plot.title <- paste0("Co-occurrence with - ", plot.title)
  plot.title <- gsub("Germinal center B cell", "GC B cell", plot.title)
  
  ggplot(dfc.sig, aes(x = correlation_discovery, y = correlation_validation, label = sig.of.interest, fill = celltype, color = celltype)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_point(aes(fill = celltype, color = celltype), size = 4) +
    geom_label_repel(show.legend = F, max.overlaps = 18, alpha = 1, colour = "white", 
                     arrow = arrow(type = "closed", angle = 20, length = unit(0.08, "inches"))) +
    scale_fill_manual(values = celltype_pal_to_use) + 
    scale_color_manual(values = celltype_pal_to_use) +
    labs(title = plot.title,
         x = "Co-occurrence in Discovery",
         y = "Co-occurrence in Validation",
         color = "Cell type") + 
    theme_pubr() +
    theme(axis.title = element_text(size = 12, face = "bold"),
          title = element_text(face = "bold")) + 
    guides(fill = "none") + 
    theme(legend.position = "none")
  
})
p.cooccur.scatter <- ggarrange(plotlist = p.cooccur.scatter.list, nrow = 2, ncol = 2, common.legend = TRUE, legend = "none")
ggsave(filename = snakemake@output[["figure2_f"]], plot = p.cooccur.scatter, width = snakemake@params[["figure2_f_width"]], height = snakemake@params[["figure2_f_height"]], units = "in", dpi = 600)

print("Figure 2F successfully created")

# plot Figure 2
design <- "
AAAABBCCCCC
AAAABBCCCCC
AAAABBCCCCC
DDDDDEEEFFF
DDDDDEEEFFF
DDDDDEEEFFF
GGGHHHIIIII
GGGHHHIIIII
GGGHHHIIIII
"

design <- "
AAAABB
AAAABB
CCCCEE
CCCCEE
DDDDFF
DDDDFF
GGGGGG
GGGGGG
"

p.corr.intra <- p.corr.intra + theme(legend.position = "none")
p.corr.inter <- p.corr.inter + theme(legend.position = "none")
p.cooccur.agree.with.ppv <- p.cooccur.agree.with.ppv + theme(legend.position = c(0.5, 1), legend.direction = "horizontal")

pdf(file = snakemake@output[["figure2_pdf"]], width = snakemake@params[["figure2_width"]], height = snakemake@params[["figure2_height"]])
schematic + p.corr.overall + p.corr.intra + p.corr.inter + p.cooccur.scatter.list[[3]] + p.cooccur.scatter.list[[4]] + p.cooccur.agree.with.ppv + 
        plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 18))
dev.off()

print("Figure 2 PDF successfully created")

CairoPNG(filename = snakemake@output[["figure2_png"]], width = snakemake@params[["figure2_width"]], height = snakemake@params[["figure2_height"]], units = "in", res = 600)
schematic + p.corr.overall + p.corr.intra + p.corr.inter + p.cooccur.scatter.list[[3]] + p.cooccur.scatter.list[[4]] + p.cooccur.agree.with.ppv + 
        plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 18))
dev.off()

print("Figure 2 PNG successfully created")