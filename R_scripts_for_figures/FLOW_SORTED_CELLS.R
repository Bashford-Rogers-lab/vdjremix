######1. Schematic
#####2. PCA 
Draw_box_plot<-function(box,x,width,c,lwd,line_col){
  segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
  segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
  segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
  rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
  segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}
Means_factor = function(factor, x){
  m = NULL
  for(i1 in c(1:length(levels(factor)))){
    x1 = x[which(factor==levels(factor)[i1])]
    x1 = x1[which(x1!=-1)]
    m = c(m, mean(x1))}
  return(m)}
Medians_factor = function(factor, x){
  m = NULL
  for(i1 in c(1:length(levels(factor)))){
    x1 = x[which(factor==levels(factor)[i1])]
    x1 = x1[which(x1!=-1)]
    m = c(m, median(x1))}
  return(m)}
concat = function(v) {
  res = ""
  for (i in 1:length(v)){res = paste0(res,v[i])}
  res
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha)) }

###########
out_dir = "~/Library/CloudStorage/OneDrive-Nexus365/DPhil Project/Ageing/FLOW/FLOW_final_matrices/"
batch = "FLOW_rWGNCA"
##############
file = concat(c(out_dir, "Eigenvectors_BCR_FLOW_BCR.txt"))
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
mat = p
heatmap(p)
p <- p[,-c(2)]

###load age and other metadata
file = concat(c(out_dir,"File information.txt"))
p1 <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
df1 = data.frame(p1)
library(dplyr)
df1 <- df1 %>% mutate(age_group = as.numeric(as.character("Age.at.visit")),
                      age_group = case_when(
                        between(age_group, 20, 40.9) ~ "20-40",
                        between(age_group, 41, 60.9) ~ "40-60",
                        between(age_group, 61, 80.9) ~ "60-80",
                        age_group >= 81 ~ "80+",
                        TRUE ~ NA_character_))

rownames(df1) <- df1$Sequencing.ID
p <- p[ order(row.names(p)), ]
df1 <- df1[ order(row.names(df1)), ]
## check 
rownames(df1) == rownames(p)
df1$cell_type <- df1$Population
mat=p


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

        
###heatmap
genotype_colors <- c("Minor Hom." = "#1B9E77", "Major Hom." = "#D95F02", "Not Significant" = "#F0F0F0")
genotype_colors <- c("Minor Hom. (TT)" = "#1B9E77", "Major Hom. (II)" = "#D95F02", "Not Significant" = "#F0F0F0") # Green and Orange from Set2, Off-White

library(dplyr)
library(tidyr)
library(ggplot2)

all_population_results <- df1 %>%
  group_by(Population) %>%
  group_map(function(pop_data, pop_info) {
    population_name <- pop_info$Population
    
    # Pivot longer and extract numeric part of module
    df_long <- pop_data %>%
      pivot_longer(
        cols = starts_with("Module_"),
        names_to = "Module",
        values_to = "Module_Score"
      ) %>%
      mutate(Module_num = as.integer(gsub("Module_", "", Module)))  # numeric version
    
    # ANOVA
    anova_results <- df_long %>%
      group_by(Module, Module_num) %>%
      do(model = aov(Module_Score ~ Genotype, data = .)) %>%
      mutate(p.value = summary(model)[[1]][["Pr(>F)"]][1]) %>%
      dplyr::select(Module, Module_num, p.value)
    
    # BH correction
    bh_adjusted_results <- anova_results %>%
      mutate(padj = p.adjust(p.value, method = "BH"))
    
    # Mean scores + significance calls
    significant_results <- df_long %>%
      group_by(Module, Module_num, Genotype) %>%
      summarise(mean_score = mean(Module_Score, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Genotype, values_from = mean_score) %>%
      left_join(bh_adjusted_results, by = c("Module", "Module_num")) %>%
      mutate(Significance = case_when(
        padj < 0.05 & !is.na(`Minor Hom. (TT)`) & !is.na(`Major Hom. (II)`) & `Minor Hom. (TT)` > `Major Hom. (II)` ~ "Minor Hom. (TT)",
        padj < 0.05 & !is.na(`Minor Hom. (TT)`) & !is.na(`Major Hom. (II)`) & `Major Hom. (II)` > `Minor Hom. (TT)` ~ "Major Hom. (II)",
        TRUE ~ "Not Significant"
      )) %>%
      dplyr::select(Module, Module_num, Significance, padj) %>%
      mutate(Population = population_name)
  }) %>%
  bind_rows()

# Get ordered Module names + numbers
module_lookup <- all_population_results %>%
  distinct(Module, Module_num) %>%
  arrange(Module_num)

# Ensure all Module-Population combinations are present, keep Module_num
all_combinations <- expand_grid(
  module_lookup,
  Population = unique(df1$Population)
)

# Build heatmap data with numeric ordering
heatmap_matrix_data <- all_combinations %>%
  left_join(all_population_results, by = c("Module", "Module_num", "Population")) %>%
  replace_na(list(Significance = "Not Significant", padj = 1)) %>%
  mutate(Module = reorder(Module, Module_num))   # enforce numeric order

# Create the matrix heatmap
matrix_heatmap_plot <- ggplot(heatmap_matrix_data, aes(y = Population, x = Module, fill = Significance)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = genotype_colors, name = "Significantly Higher In (BH)") +
  labs(
    y = "Population",
    x = "Module",
    title = "Module Score Significance Across Populations (BH Corrected)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.y = element_blank()
  )

# Display the matrix heatmap
print(matrix_heatmap_plot)


pdf("heatmap_new.pdf", 22, 3)
matrix_heatmap_plot
dev.off()

