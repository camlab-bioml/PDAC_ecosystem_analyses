suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(stringr)
  library(corrplot)
  library(ComplexHeatmap)
  library(patchwork)
})

w.corr <- read_tsv(snakemake@input[['gene_loading_corr']])
rownames(w.corr) <- names(w.corr)

w.dist <- read_tsv(snakemake@input[['gene_loading_dist']])
rownames(w.dist) <- names(w.dist)

celltype = snakemake@wildcards[['subtype']]

w.dist.sim.list <- readRDS(snakemake@input[['simulated_dist_list']])
w.corr.sim.list <- readRDS(snakemake@input[['simulated_corr_list']])

# plot distance between signatures
p.list <- lapply(seq(length(w.dist.sim.list)), function(sig.dis.id) {
  dis.dist.sim.list <- w.dist.sim.list[[sig.dis.id]]
  the.list <- lapply(seq(length(dis.dist.sim.list)), function(sig.val.id) {
    w.dist.sim <- dis.dist.sim.list[[sig.val.id]]
    
    ggplot(data = data.frame(simulated.distance = w.dist.sim,
                             five.percent = quantile(w.dist.sim, 0.05),
                             one.percent = quantile(w.dist.sim, 0.01),
                             real.distance = as.matrix(w.dist)[paste("discovery", sig.dis.id, sep = " "), paste("validation", sig.val.id, sep = " ")]), 
           mapping = aes(x = simulated.distance)) +
      geom_histogram(aes(y = after_stat(density)), color = "darkblue", fill = "lightblue") +
      geom_vline(aes(xintercept = mean(five.percent)), 
                 color = "blue", linetype = "dashed", linewidth = 1) +
      geom_text(aes(x = mean(five.percent), label = paste0("\n", "5%"), y = 0), 
                colour = "blue", angle = 90, hjust = 0, vjust = 0.5) +
      geom_vline(aes(xintercept = mean(one.percent)), 
                 color = "orange", linetype = "dashed", linewidth = 1) +
      geom_text(aes(x = mean(one.percent), label = paste0("\n", "1%"), y = 0), 
                colour = "orange", angle = 90, hjust = 0, vjust = 0.5) +
      geom_vline(aes(xintercept = mean(real.distance)), 
                 color = "red", linetype = "solid", linewidth = 1) +
      geom_text(aes(x = mean(real.distance), label = paste0("\n", "real distance: ", round(mean(real.distance),0)), y = 0), 
                colour = "red", angle = 90, hjust = 0, vjust = 0.5) +
      geom_density(alpha = .2, fill = "#FF6666") +
      theme_pubr() + 
      labs(title = paste0(celltype, " Discovery.", sig.dis.id, " vs. ", "Validation.", sig.val.id),
           x = "Simulated distance") +
      labs(title = paste0("Validation.", sig.val.id),
           x = NULL,
           y = NULL) 
  })
  names(the.list) <- paste("validation", seq(length(dis.dist.sim.list)), sep = " ")
  the.list
})
names(p.list) <- paste("discovery", seq(length(w.dist.sim.list)), sep = " ")


pdf(file = snakemake@output[['gene_loading_sim_dist_plot']], width = 15, height = 10, onefile = T)
for (plot.sig.id in seq(length(p.list))) {
  print(wrap_plots(plotlist = p.list[[paste0("discovery ", plot.sig.id)]]) +
          plot_annotation(title = paste0(celltype, " Discovery.", plot.sig.id, " vs. "),
                          subtitle = paste0("x = Simulated distance, y = Density"),
                          theme = theme(plot.title = element_text(size = 18, face = "bold"))))
}
dev.off()


## plot correlation between signatures
p.list <- lapply(seq(length(w.corr.sim.list)), function(sig.dis.id) {
  dis.corr.sim.list <- w.corr.sim.list[[sig.dis.id]]
  the.list <- lapply(seq(length(dis.corr.sim.list)), function(sig.val.id) {
    w.corr.sim <- dis.corr.sim.list[[sig.val.id]]
    
    ggplot(data = data.frame(simulated.correlation = w.corr.sim,
                             five.percent = quantile(w.corr.sim, 0.95),
                             one.percent = quantile(w.corr.sim, 0.99),
                             real.correlation = as.matrix(w.corr)[paste("discovery", sig.dis.id, sep = " "), paste("validation", sig.val.id, sep = " ")]), 
           mapping = aes(x = simulated.correlation)) +
      geom_histogram(aes(y = after_stat(density)), color = "darkblue", fill = "lightblue") +
      geom_vline(aes(xintercept = mean(five.percent)), 
                 color = "blue", linetype = "dashed", linewidth = 1) +
      geom_text(aes(x = mean(five.percent), label = paste0("\n", "5%"), y = 0), 
                colour = "blue", angle = 90, hjust = 0, vjust = 0.5) +
      geom_vline(aes(xintercept = mean(one.percent)), 
                 color = "orange", linetype = "dashed", linewidth = 1) +
      geom_text(aes(x = mean(one.percent), label = paste0("\n", "1%"), y = 0), 
                colour = "orange", angle = 90, hjust = 0, vjust = 0.5) +
      geom_vline(aes(xintercept = mean(real.correlation)), 
                 color = "red", linetype = "solid", linewidth = 1) +
      geom_text(aes(x = mean(real.correlation), label = paste0("\n", "real correlation: ", round(mean(real.correlation),4)), y = 0), 
                colour = "red", angle = 90, hjust = 0, vjust = 0.5) +
      geom_density(alpha = .2, fill = "#FF6666") +
      theme_pubr() + 
      labs(title = paste0(celltype, " Discovery.", sig.dis.id, " vs. ", "Validation.", sig.val.id),
           x = "Simulated correlation") +
      labs(title = paste0("Validation.", sig.val.id),
           x = NULL,
           y = NULL) 
  })
  names(the.list) <- paste("validation", seq(length(dis.corr.sim.list)), sep = " ")
  the.list
})
names(p.list) <- paste("discovery", seq(length(w.corr.sim.list)), sep = " ")


pdf(file = snakemake@output[['gene_loading_sim_corr_plot']], width = 15, height = 10, onefile = T)
for (plot.sig.id in seq(length(p.list))) {
  print(wrap_plots(plotlist = p.list[[paste0("discovery ", plot.sig.id)]]) +
          plot_annotation(title = paste0(celltype, " Discovery.", plot.sig.id, " vs. "),
                          subtitle = paste0("x = Simulated correlation, y = Density"),
                          theme = theme(plot.title = element_text(size = 18, face = "bold"))))
}
dev.off()


