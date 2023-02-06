suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
  library(corrplot)
  library(circlize)
  library(dendextend)
})

profile.flavor = snakemake@wildcards[['profile']]
condition = snakemake@wildcards[['condition']]
compartment = snakemake@wildcards[['compartment']]

# load patient profiles
sigprofile.df <- read_tsv(snakemake@input[['patient_profiles']])
sigprofile.df.long <- read_tsv(snakemake@input[['patient_profiles_long']])
sigprofile.corr <- read_tsv(snakemake@input[['patient_profiles_corr']]) %>% 
  as.matrix()
rownames(sigprofile.corr) <- colnames(sigprofile.corr)

sigprofile.mtx <- sigprofile.df %>% select(!c(cohort, sample)) %>% as.matrix()

# simple heatmap ---------------------------------------------------------------
row_ha <- rowAnnotation(Cohort = sigprofile.df$cohort)
col_ha <- columnAnnotation(Celltype = (sigprofile.df %>% select(!c(cohort, sample)) %>% names() %>% str_split(., pattern = " ", simplify = T))[,1])

png(snakemake@output[['patient_profiles_heatmap_grouped']], width = 12, height = 7, units = "in", res = 321)
Heatmap(sigprofile.mtx %>% scale(), 
        name = paste("Norm sig", profile.flavor, sep = " "),
        #left_annotation = row_ha, 
        row_split = sigprofile.df$cohort,
        #top_annotation = col_ha,
        column_split = (sigprofile.df %>% select(!c(cohort, sample)) %>% names() %>% str_split(., pattern = " ", simplify = T))[,1]
        )
dev.off()

png(snakemake@output[['patient_profiles_heatmap_clustered']], width = 12, height = 7, units = "in", res = 321)
Heatmap(sigprofile.mtx %>% scale(), 
        name = paste("Norm sig", profile.flavor, sep = " "),
        left_annotation = row_ha, 
        #row_split = sigprofile.df$cohort,
        top_annotation = col_ha,
        #column_split = (sigprofile.df %>% select(!c(cohort, sample)) %>% names() %>% str_split(., pattern = " ", simplify = T))[,1]
        )
dev.off()

# signature correlation across samples - corrplot ------------------------------
png(snakemake@output[['patient_profiles_corrplot']], width = 12, height = 12, units = "in", res = 321)
corrplot(corr = sigprofile.corr)
dev.off()

# faceted stacked bar plot for the top frequency flavor ------------------------
ggplot(sigprofile.df.long, aes(x = sample, y = get(profile.flavor), fill = signature)) + 
  geom_bar(position = "stack", stat = "identity") + 
  facet_wrap(~ celltype, ncol = 1, scales = "free") +
  theme_pubr(x.text.angle = 45) + 
  theme(axis.text.x = element_blank(),
        legend.position = "right") + 
  labs(y = profile.flavor)
ggsave(snakemake@output[['patient_profiles_stacked_bar_plot']], device = "png", height = 15, width = 10, units = "in", dpi = 321)

# plot circlized bars ----------------------------------------------------------
sigprofile.plist <- lapply(unique(sigprofile.df.long$sample), function(s) {
  # ----- This section prepare a dataframe for plotting ---- #
  sigprofile.df.long.sample <- sigprofile.df.long %>% filter(sample == s) %>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0))
  sigprofile.df.long.sample$celltype <- factor(sigprofile.df.long.sample$celltype, levels = unique(sigprofile.df.long.sample$celltype))
  
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 1
  to_add <- data.frame(matrix(NA, empty_bar*nlevels(sigprofile.df.long.sample$celltype), ncol(sigprofile.df.long.sample)))
  colnames(to_add) <- colnames(sigprofile.df.long.sample)
  to_add$celltype <- rep(levels(sigprofile.df.long.sample$celltype), each = empty_bar)
  sigprofile.df.long.sample <- rbind(sigprofile.df.long.sample, to_add)
  sigprofile.df.long.sample <- sigprofile.df.long.sample %>% arrange(celltype)
  sigprofile.df.long.sample$id <- seq(1, nrow(sigprofile.df.long.sample))
  
  # Get the name and the y position of each label
  number_of_bar <- nrow(sigprofile.df.long.sample)
  angle <- 90 - 360 * (sigprofile.df.long.sample$id-0.5) /number_of_bar     # I subtract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  sigprofile.df.long.sample$hjust <- ifelse( angle < -90, 1, 0)
  sigprofile.df.long.sample$angle <- ifelse(angle < -90, angle+180, angle)
  # ----- ------------------------------------------- ---- #
  
  # prepare a data frame for base lines
  base_data <- sigprofile.df.long.sample %>% 
    group_by(celltype) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  
  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1.5
  grid_data$start <- grid_data$end - 1
  #grid_data <- grid_data[-1,]
  
  # get position/value for the indicator lines
  line_guides <- summary(sigprofile.df.long.sample[[profile.flavor]])
  line_max <- line_guides['Max.']
  
  sigprofile.df.long.sample[[profile.flavor]][sigprofile.df.long.sample[[profile.flavor]] == 0] = line_max*0.01
  
  # make the circlized plot
  p <- ggplot(sigprofile.df.long.sample, aes(x = as.factor(id), y= get(profile.flavor), fill = celltype)) + 
    geom_bar(aes(x = as.factor(id), y= get(profile.flavor), fill = celltype), stat = "identity", alpha = 0.5) + 
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sure barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = line_max, xend = start, yend = line_max), 
                 colour = "grey", alpha=1, size=0.3, inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = line_max*0.75, xend = start, yend = line_max*0.75), 
                 colour = "grey", alpha=1, size=0.3, inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = line_max*0.5, xend = start, yend = line_max*0.5), 
                 colour = "grey", alpha=1, size=0.3, inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = line_max*0.25, xend = start, yend = line_max*0.25), 
                 colour = "grey", alpha=1, size=0.3, inherit.aes = FALSE ) +
    # Add text showing the value of each 100/75/50/25 lines
    annotate("text", 
             x = rep(max(sigprofile.df.long.sample$id),4) + 0.5, 
             y = seq(from = line_max*0.25, to = line_max, by = line_max*0.25) + line_max*0.05, 
             label = format(seq(from = line_max*0.25, to = line_max, by = line_max*0.25), digits = 2),
             color="grey", size=2, angle=0, fontface="bold", hjust=1) +
    geom_bar(aes(x = as.factor(id), y = get(profile.flavor), fill = celltype), stat = "identity", alpha = 0.5) + 
    # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
    ylim(-line_guides['Max.'], NA) + 
    # Custom the theme: no axis title and no cartesian grid
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
    ) +
    # This makes the coordinate polar instead of cartesian.
    coord_polar(start = 0) +
    # Add the labels, using the dataframe that we have created before
    geom_text(data=sigprofile.df.long.sample,
              aes(x = id, y = get(profile.flavor)+line_max*0.02, label = signature, hjust=hjust),
              color="black", fontface="bold", alpha=0.6, size=2.5, angle= sigprofile.df.long.sample$angle, inherit.aes = FALSE ) +
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -line_max*0.1, xend = end, yend = -line_max*0.1), 
                 colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
    geom_text(data=base_data, aes(x = title, y = -line_max*0.15, label = celltype), hjust = c(1,1,1,1,0,0,0,0), 
              colour = "black", alpha=0.8, size=1.5, fontface="bold", inherit.aes = FALSE)
  
  p
})
names(sigprofile.plist) <- unique(sigprofile.df.long$sample)

png(snakemake@output[['patient_profiles_circlized_bar_plot']], width = 30, height = 30, unit = "in", res = 321)
cowplot::plot_grid(plotlist = sigprofile.plist, ncol = floor(sqrt(length(sigprofile.plist))))
dev.off()




















