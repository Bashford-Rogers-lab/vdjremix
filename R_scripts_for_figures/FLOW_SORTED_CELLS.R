# Publication-ready analysis script for FLOW Figure 4
# Panel A (schematic) is not generated in R here.
# This script covers panel B (PCA-style sample projection) and panel C (genotype-by-population heatmap).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(labdsv)
  library(readr)
  library(vdjremix)
})

set.seed(1974)


files <- list(
  eigenvectors = file.path(base_dir, "Eigenvectors_BCR_FLOW_BCR.txt"),
  metadata = file.path(base_dir, "File information.txt")
)

output_files <- list(
  pca_plot = file.path(base_dir, "PCA_dimred_patients_FLOW_rWGNCA.pdf"),
  mapping_table = file.path(base_dir, "All_mapping_dimred_patients_FLOW_rWGNCA.txt"),
  heatmap = file.path(base_dir, "heatmap_new.pdf")
)

# ---- Load and align data ----
module_df <- read.delim(files$eigenvectors, check.names = FALSE)

# The original workflow dropped column 2 from the exported eigenvector table.
# Keep that behaviour here, but make it explicit.
if (ncol(module_df) < 3) {
  stop("The eigenvector table has fewer than 3 columns; cannot drop the non-module column used in the original workflow.")
}

# Use the first column as sample ID when it is not already stored as row names.
if (is.null(rownames(module_df)) || all(rownames(module_df) == as.character(seq_len(nrow(module_df))))) {
  rownames(module_df) <- module_df[[1]]
}

module_df <- module_df[, -2, drop = FALSE]
module_matrix <- module_df %>%
  mutate(across(everything(), ~ as.numeric(as.character(.x)))) %>%
  as.matrix()
rownames(module_matrix) <- rownames(module_df)
module_matrix <- module_matrix[order(rownames(module_matrix)), , drop = FALSE]

metadata <- read.delim(files$metadata, check.names = FALSE) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  mutate(
    age_numeric = as.numeric(`Age.at.visit`),
    age_group = case_when(
      between(age_numeric, 20, 40.9) ~ "20-40",
      between(age_numeric, 41, 60.9) ~ "40-60",
      between(age_numeric, 61, 80.9) ~ "60-80",
      age_numeric >= 81              ~ "80+",
      TRUE                           ~ NA_character_
    )
  )

rownames(metadata) <- metadata$Sequencing.ID
metadata <- metadata[order(rownames(metadata)), , drop = FALSE]

common_ids <- intersect(rownames(module_matrix), rownames(metadata))
module_matrix <- module_matrix[common_ids, , drop = FALSE]
metadata <- metadata[common_ids, , drop = FALSE]
stopifnot(identical(rownames(module_matrix), rownames(metadata)))

analysis_df <- cbind(
  metadata,
  as.data.frame(module_matrix, check.names = FALSE)
)
analysis_df$cell_type <- analysis_df$Population

# ---- Panel B: sample projection plot ----

Plot_dimred_for_patients<-function(p){
  library(umap)
  library(labdsv)
  
  x <- pca(t(mat),dim=3)
  dimred = x$loadings
  dimred_umap= x$loadings
  #dimred= mat
  #umap = umap(dimred , n_components = 2)
  #dimred_umap = umap$ layout
  db <- kmeans(dimred_umap, 8)
  cluster.umap = as.numeric(db$ cluster)
  
  clusters = sort(unique(cluster.umap))
  cluster_match = match(cluster.umap, clusters)
  id_raw = rownames(mat)
  cell_types = unique(df1$Population)
  genotypes = unique(df1$Genotype)
  cell_type = df1$Population
  genotype = df1$Genotype
  
  mapping = cbind(rownames(mat), cell_type, genotype, cluster.umap)
  out_file_table = concat(c(out_dir, "All_mapping_dimred_patients_",batch,".txt"))
  write.table(mapping, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  fileout1=concat(c(out_dir, "PCA_dimred_patients_",batch,".pdf"))
  w=2.9
  pdf(file=fileout1, height=w*2*1, width=w*3)
  par(mfrow= c(2,3), mar = c(5,5,3,3))
  
  rownames(dimred) = rownames(mat)
  
  library(RColorBrewer)
  pches = rep(c(21:25), 100)
  cols =  rep(add.alpha (c(brewer.pal(6, "Dark2"),brewer.pal(8, "Paired"),brewer.pal(8, "Set2")), alpha = 0.5), 100)
  cols = cols[-5]
  
  col_match = match(cell_type, cell_types)
  pch_match = match(genotype, genotypes)
  
  cols_geno = c("white", "grey")
  
  #cluster by cell type
  plot(dimred_umap[,c(1,2)], pch = pches[col_match], col = cols[col_match], bg = cols[col_match], cex =1.2,xlab = "Dim1", ylab = "Dim2", main = "coloured by cell type")
  
  plot(dimred_umap[,c(1,2)], pch = pches[col_match], col = cols[col_match], bg = cols_geno[pch_match], cex =1.2,xlab = "Dim1", ylab = "Dim2", main = "coloured by cell type and genotype")
  
  plot(c(1,2),c(1,2), col = "white", main = '', xlab = "", ylab = "", cex = 1,lwd = 2,axes = F)
  legend("bottomleft", cell_types, pch = pches,cex= 0.8, bty="n", pt.bg = cols, col = NA, pt.lwd = 1, text.font = 1)
  
  #cluster by genotype
  legend("bottomright", genotypes, pch = 21,cex= 0.8, bty="n", pt.bg = cols_geno, col = "black", pt.lwd = 1, text.font = 1)
  
  # plot(dimred_umap[,c(1,2)], pch = pches[cluster_match], col = cols[cluster_match], bg = cols[cluster_match], cex =1.2,xlab = "Dim1", ylab = "Dim2", main = "clustered features")
  # plot(c(1,2),c(1,2), col = "white", main = '', xlab = "", ylab = "", cex = 1,lwd = 2,axes = F)
  # legend("topleft", apply(cbind("cluster ",clusters), 1, paste, collapse = ""), pch = pches,cex= 0.8, bty="n", pt.bg = cols, col = cols,  pt.lwd = 1, text.font = 1)
  
  dev.off()
  
}
Plot_dimred_for_patients(p)


# ---- Panel C: genotype effect heatmap across populations ----
run_population_module_heatmap <- function(analysis_df) {
  genotype_colors <- c(
    "Minor Hom. (TT)" = "#1B9E77",
    "Major Hom. (II)" = "#D95F02",
    "Not Significant" = "#F0F0F0"
  )

  all_population_results <- analysis_df %>%
    group_by(Population) %>%
    group_map(function(pop_data, pop_info) {
      population_name <- pop_info$Population

      long_df <- pop_data %>%
        pivot_longer(
          cols = starts_with("Module_"),
          names_to = "Module",
          values_to = "Module_Score"
        ) %>%
        mutate(Module_num = as.integer(gsub("Module_", "", Module)))

      anova_results <- long_df %>%
        group_by(Module, Module_num) %>%
        do({
          fit <- aov(Module_Score ~ Genotype, data = .)
          data.frame(p_value = summary(fit)[[1]][["Pr(>F)"]][1])
        }) %>%
        ungroup() %>%
        mutate(padj = p.adjust(p_value, method = "BH"))

      long_df %>%
        group_by(Module, Module_num, Genotype) %>%
        summarise(mean_score = mean(Module_Score, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = Genotype, values_from = mean_score) %>%
        left_join(anova_results, by = c("Module", "Module_num")) %>%
        mutate(
          Significance = case_when(
            padj < 0.05 & !is.na(`Minor Hom. (TT)`) & !is.na(`Major Hom. (II)`) & `Minor Hom. (TT)` > `Major Hom. (II)` ~ "Minor Hom. (TT)",
            padj < 0.05 & !is.na(`Minor Hom. (TT)`) & !is.na(`Major Hom. (II)`) & `Major Hom. (II)` > `Minor Hom. (TT)` ~ "Major Hom. (II)",
            TRUE ~ "Not Significant"
          ),
          Population = population_name
        ) %>%
        select(Module, Module_num, Significance, padj, Population)
    }) %>%
    bind_rows()

  module_lookup <- all_population_results %>%
    distinct(Module, Module_num) %>%
    arrange(Module_num)

  heatmap_df <- tidyr::expand_grid(
    module_lookup,
    Population = unique(analysis_df$Population)
  ) %>%
    left_join(all_population_results, by = c("Module", "Module_num", "Population")) %>%
    replace_na(list(Significance = "Not Significant", padj = 1)) %>%
    mutate(Module = factor(Module, levels = module_lookup$Module))

  heatmap_plot <- ggplot(heatmap_df, aes(x = Module, y = Population, fill = Significance)) +
    geom_tile(color = "black", linewidth = 0.25) +
    scale_fill_manual(values = genotype_colors, name = "Significantly higher in (BH)") +
    labs(
      x = "Module",
      y = "Population",
      title = "Module score significance across populations"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(face = "bold")
    )

  ggsave(output_files$heatmap, plot = heatmap_plot, width = 22, height = 3)
  heatmap_plot
}

# ---- Helper used in the original plotting workflow ----
add.alpha <- function(col, alpha = 1) {
  if (missing(col)) {
    stop("Please provide a vector of colours.")
  }
  apply(sapply(col, col2rgb) / 255, 2, function(x) {
    rgb(x[1], x[2], x[3], alpha = alpha)
  })
}

# ---- Run figure panels ----
pca_plot <- run_sample_projection(module_matrix, metadata)
print(pca_plot)

heatmap_plot <- run_population_module_heatmap(analysis_df)
print(heatmap_plot)
