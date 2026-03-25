library(vdjremix)
mega_matrix1= readRDS(file = concat(c(output_directory, "VDJ_REMIX_megamatrix.rds")))

res <- run_vdjremix(
  feature_matrix = VDJ_REMIX_megamatrix,
  na_frequencies = c(0.8),
  min_cluster_size = 10,
  n_repeats = 20,
  n_cores = parallel::detectCores() - 1
)

#key outputs
modules = res$modules            # feature modules
mat_remix = res$eigengenes         # module eigengenes 
var_explained = res$variance_explained # variance explained by each module
loadings = res$loadings           # feature contributions within modules
rownames(mat_remix) = rownames(VDJ_REMIX_megamatrix)
### which associate with PDAC group
#### tumour

ids_association = c(groups_PCA[[1]], groups_PCA[[2]])
factor = c(rep("ME", length(groups_PCA[[1]])), rep("AE", length(groups_PCA[[2]])))
factor= factor(factor)
mat_stat = mat_remix
colnames(mat_stat) = paste("PDAC_mod.", colnames(mat_stat))
fit = manova(formula = mat_stat[ids_association,] ~ factor )

p1 = summary.aov(fit)
nam = gsub(" Response ","",names(p1))
p_value = NULL
means = NULL
i1 = 0
for(i in p1){
  i1 = i1+1
  p_value = c(p_value, i$'Pr(>F)'[1]) 
  if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
  }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
}
p_value[which(is.na(p_value))] = 2
names(p_value) = nam
print(min(p_value))
#print(length(which(p_value<0.05)))
colnames(means) = paste("mean.group.", c(1:length(means[1,])))
combined_p_value = cbind(p_value ,means)
rownames(combined_p_value) = nam
p.group = rep("REMIX_PDAC", length(nam))
summary = cbind(nam, p.group,combined_p_value)

signif_modules = names(which(p_value<0.05))
modules[gsub("PDAC_mod. ","", signif_modules)]
p_value[signif_modules]
### plot modules and stats on components
Module_stats<-function(signif_modules, modules, mega_matrix1){
  features = unlist(modules[gsub("PDAC_mod. ","", signif_modules)])
  factor= factor(factor)
  mat_stat = mega_matrix1[,features]
  colnames(mat_stat) = paste("PDAC_mod.", colnames(mat_stat))
  fit = manova(formula = mat_stat[ids_association,] ~ factor )
  
  p1 = summary.aov(fit)
  nam = gsub(" Response ","",names(p1))
  p_value = NULL
  means = NULL
  i1 = 0
  for(i in p1){
    i1 = i1+1
    p_value = c(p_value, i$'Pr(>F)'[1]) 
    if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
    }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
  }
  p_value[which(is.na(p_value))] = 2
  names(p_value) = nam
  print(min(p_value))
  #print(length(which(p_value<0.05)))
  colnames(means) = paste("mean.group.", c(1:length(means[1,])))
  combined_p_value = cbind(p_value ,means)
  rownames(combined_p_value) = nam
  p.group = rep("REMIX_PDAC", length(nam))
  summary = cbind(nam, p.group,combined_p_value)
  
  summary[which(p_value<0.05),]
}


library(ggpubr)
library(ggplot2)


sample_colors <- c("Myeloid enriched PDAC" = "#E7298AFF", "Myeloid enriched PDAC" = "#17BECFFF")


for (mod_name in signif_modules) {
  

  clean_name <- gsub("PDAC_mod. ", "", mod_name)
  

  df_plot <- data.frame(
    score = c(mat_remix[groups_PCA[[1]], clean_name], 
              mat_remix[groups_PCA[[2]], clean_name]),
    group = factor(c(rep("Myeloid enriched PDAC", length(groups_PCA[[1]])), 
                     rep("Adaptive enriched PDAC", length(groups_PCA[[2]]))))
  )
  

  df_plot <- df_plot[!is.na(df_plot$score), ]
  

  p <- ggboxplot(
    df_plot,
    x = "group",
    y = "score",
    color = "group",
    fill = "group",
    alpha = 0.2,          # Soft fill for the box
    add = "jitter",       # Add individual points
    add.params = list(size = 1.2, alpha = 0.7),
    width = 0.5,          # Tighter box width
    palette = sample_colors
  ) +
    # Add significance comparison
    stat_compare_means(
      comparisons = list(c("Group 1", "Group 2")),
      method = "t.test",   # or "wilcox.test"
      label = "p.signif",
      bracket.size = 0.5
    ) +
    # Clean labels and theme
    labs(
      title = mod_name,
      subtitle = paste0("ANOVA p = ", formatC(p_value[mod_name], format = "e", digits = 2)),
      x = "Patient Group",
      y = "Module Score"
    ) +
    theme_pubr() +   
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) +
    theme(legend.position = "none") 
  # 5. Save the output
  ggsave(
    filename = paste0("Boxplot_", clean_name, "_gg.pdf"), 
    plot = p, 
    width = 4, 
    height = 5
  )
}
