suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(readr)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
  library(corrplot)
  library(patchwork)
})

# load cell type rename csv
cell_type_rename <- read_csv(snakemake@input[["cell_type_rename"]])

# load signature loadings in patients and in single cells
sig.profiles <- read_tsv(snakemake@input[["sig_profiles"]])
sig.loadings <- read_tsv(snakemake@input[["sig_loading_mtx"]])

names(sig.profiles) <- gsub(" Rep | RepVal ", " ", names(sig.profiles))
names(sig.loadings) <- gsub(" Rep | RepVal ", " ", names(sig.loadings))

# load df.redim
df.redim <- read_tsv(snakemake@input[["dimred_with_top_two_sig_loadings"]])

# load signature interpretation
sig.interpt <- readxl::read_xlsx(snakemake@input[["sig_interpretation"]])

# signatures to remove
ambient.sigs <- read_csv(snakemake@input[["ambient_sigs"]])$signature
ambient.sigs <- gsub(" Sig ", " ", ambient.sigs)

celltype <- snakemake@wildcards[["subtype"]]
condition <- snakemake@wildcards[["condition"]]
profile.flavour <- snakemake@wildcards[["profile"]]
dim.red.plot <- snakemake@params[["dim_red_plot"]]


# tidy up df.redim
df.redim <- df.redim |>
  mutate(
    top_sig = plyr::mapvalues(top_sig, from = sig.interpt$signature, to = sig.interpt$`short interpretation`, warn_missing = FALSE),
    second_sig = plyr::mapvalues(second_sig, from = sig.interpt$signature, to = sig.interpt$`short interpretation`, warn_missing = FALSE)
  ) |>
  mutate(
    top_sig = ifelse(grepl("Ambient RNA", top_sig), NA, top_sig),
    second_sig = ifelse(grepl("Ambient RNA", second_sig), NA, second_sig)
  )

# plot top two signatures in UMAP
p.top.sig <- ggplot(df.redim, aes(x = UMAP_1, y = UMAP_2, color = top_sig)) +
  geom_point(alpha = 0.3, size = 0.4, shape = 19) +
  # scale_color_manual(values = c(celltype.pal.dis, celltype.pal.val)) +
  theme_pubr() +
  labs(x = paste0(dim.red.plot, " 1"), y = paste0(dim.red.plot, " 2"), colour = "Top signature") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8, shape = 16))) +
  theme(
    axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold")
  )

p.second.sig <- ggplot(df.redim, aes(x = UMAP_1, y = UMAP_2, color = second_sig)) +
  geom_point(alpha = 0.3, size = 0.4, shape = 19) +
  # scale_color_manual(values = cohort.pal) +
  theme_pubr() +
  labs(x = paste0(dim.red.plot, " 1"), y = paste0(dim.red.plot, " 2"), colour = "Second signature") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8, shape = 16))) +
  theme(
    axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold")
  )

png(filename = snakemake@output[["top_two_sigs_umap"]], width = 12, height = 5, units = "in", res = 360)
p.top.sig + p.second.sig + plot_layout(guides = "auto", nrow = 1) & theme(legend.position = "right")
dev.off()

# plot signature loading patterns in UMAP
df.redim.pattern <- df.redim |>
  select(Cell_ID, Sample, Cohort, Cell_type, UMAP_1, UMAP_2)

names(sig.loadings) <- gsub("cohort", "Cohort", names(sig.loadings))
names(sig.loadings) <- gsub("sample", "Sample", names(sig.loadings))
names(sig.loadings) <- gsub("cell_id", "Cell_ID", names(sig.loadings))

df.redim.pattern <- left_join(df.redim.pattern, sig.loadings, by = c("Cell_ID", "Sample", "Cohort"))

## tidy up df.redim.pattern
df.redim.pattern <- df.redim.pattern |>
  pivot_longer(contains(celltype), names_to = "Signature", values_to = "Loading") |>
  group_by(Signature) |>
  mutate(Loading = scale(Loading)) |>
  ungroup() |>
  filter(!(Signature %in% ambient.sigs)) |>
  mutate(
    Signature = plyr::mapvalues(
      Signature, 
      from = sig.interpt$signature, 
      to = sig.interpt$`short interpretation`, 
      warn_missing = FALSE
      )
  )

## plot
p.loading.pattern <- ggplot(df.redim.pattern, aes(x = UMAP_1, y = UMAP_2, color = Loading)) +
  geom_point(alpha = 0.4, size = 0.2, shape = 19) +
  facet_wrap( ~ Signature) +
  scale_color_viridis_c(option = "cividis", direction = -1) +
  theme_pubr() +
  labs(x = "", y = "", colour = "Signature loading") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8, shape = 16)), x = "none", y = "none") +
  theme(
    axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold")
  )

ggsave(
  snakemake@output[["sig_loading_pattern_umap"]],
  plot = p.loading.pattern,
  device = "png",
  width = 12, height = 12, units = "in", dpi = 360, bg = "white"
  )


# plot signature loading correlations
png(filename = snakemake@output[["sig_loading_pattern_corrplot"]], width = 7, height = 7, units = "in", res = 360)
corrplot(cor(sig.profiles %>% select(contains(celltype)) %>% select(order(names(.))) %>% as.matrix() %>% scale()))
dev.off()

png(filename = snakemake@output[["sig_loading_pattern_single_cell_corrplot"]], width = 7, height = 7, units = "in", res = 360)
corrplot(cor(sig.loadings %>% select(contains(celltype)) %>% select(order(names(.))) %>% as.matrix() %>% scale()))
dev.off()
