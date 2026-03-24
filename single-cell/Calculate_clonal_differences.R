
Boxplot_custom<-function(groups, main, width_plot, colsx){
  factors = names(groups)
  max = max(c(unlist(groups), unlist(groups))*1.35)
  min = 0
  if(min(unlist(groups))<0){
    min = -1}
  b = (max-min)*0.034
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = width_plot
  max_scale = min(c(max,100))
  range = max-min
  if(range>50){scale = c(-100:100)*20}
  if(range>200){scale = c(-100:100)*50}
  if(range<=50){scale = c(-100:100)*10}
  if(range<=30){scale = c(-100:100)*5}
  if(range <10){scale = c(-100:100)*2.5}
  if(range <5){scale = c(-100:100)*1}
  if(range <4){scale = c(-100:100)*0.5}
  if(range <1.5){scale = c(-100:1000)*0.2}
  if(range <0.5){scale = c(-100:100)*0.1}
  if(range <0.1){scale = c(-100:100)*0.01}
  if(range <0.01){scale = c(-100:100)*0.001}
  cex = 0.9
  Fun<-function(x){x}
  
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.38
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  
  for(i in c(1:l)){
    points1=as.numeric(groups[[i]])
    box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
    Draw_box_plot(box1,i,width,colsx[i],1, colsx1[i])
    points(rep(i, length(points1)),points1, pch =21, col=colsx[i],bg = colsx1[i], cex = 0.7)
  }
}

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


Draw_box_plot<-function(box,x,width,c,lwd,line_col){
	segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
	segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
	segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
	rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
	segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}
	
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha)) }

concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}

######################
batch = "PDAC150Ka"
analysis = "VDJ_clonality"

######################
## get patient groups for plotting: 
input_directory = "/PDAC150K/DATA/OUTPUTS/"
output_directory = "/GITHUB/PDAC150K/DATA/OUTPUTS/"

input_directory_groups  = "/PDAC150K/DATA/PROPORTIONS/"
groups_PCA = readRDS(file = concat(c(input_directory_groups, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))

library(yarrr)
library(RColorBrewer)

Get_megamatrix<-function(output_directory, batch){
  file_BCR=concat(c(output_directory, "VDJ_Clonality_information_BCR_", batch,".txt"))
  p_BCR <- as.matrix(read.csv(file_BCR, head=TRUE, sep="\t"))
  
  file = concat(c(input_directory, "VDJ_repertoire_feature_information_BCR_",batch,".rds"))
  VDJ_BCR = readRDS(file= file)
  
  file_TCR=concat(c(output_directory, "VDJ_Clonality_information_TCR_CD4TCR_CD8_", batch,".txt"))
  p_TCR <- as.matrix(read.csv(file_TCR, head=TRUE, sep="\t"))
  
  cn = colnames(p_TCR)
  exclude = c(intersect(grep("TCR_CD8_", cn), grep("CD4", cn)), 
              intersect(grep("TCR_CD4_", cn), grep("CD8", cn)))
  
  p_TCR = p_TCR[,setdiff(cn, cn[exclude])]
  
  file = concat(c(input_directory, "VDJ_repertoire_feature_information_TCR_",batch,".rds"))
  VDJ_TCR = readRDS(file= file)

  #calculated previously using Cell_type_pathways_per_sample.R file
  file = "Pathways_mean_pathway_score_per_cell_type_per_sample_PDAC150Ka_raw_pathway_level.rds"
  pathways = readRDS(file= file)
  
  pw = names(pathways)
  include = c("KEGG_APOPTOSIS", "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY" ,"KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY" ,
              "HALLMARK_HYPOXIA", "GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_T_CELL_PROLIFERATION",
              "GO_B_CELL_PROLIFERATION" , "GO_RESPONSE_TO_INTERFERON_BETA",
              "GO_RESPONSE_TO_INTERFERON_ALPHA" , "GO_RESPONSE_TO_INTERFERON_GAMMA","GO_CELLULAR_RESPONSE_TO_INTERFERON_BETA"  ,
              "GO_INTERFERON_GAMMA_PRODUCTION",
              "GO_B_CELL_RECEPTOR_SIGNALING_PATHWAY"   )
  pw = pw[which(pw %in% include)]
  pathways = pathways[pw]
  ### combine all files to one mega matrix
  names(VDJ_BCR) = paste("BCR", names(VDJ_BCR))
  names(VDJ_TCR) = paste("TCR", names(VDJ_TCR))
  colnames(p_BCR) = paste("BCR", colnames(p_BCR))
  colnames(p_TCR) = paste("TCR", colnames(p_TCR))
  
  VDJ_BCR = VDJ_BCR[grep("gene usage", names(VDJ_BCR),invert = T)]
  
  list_mats1 = c(list(p_BCR), list(p_TCR))
  list_mats2 = c(VDJ_BCR)
  
  ids = rownames(p_BCR)
  mega_matrix = NULL
  for(i in c(1:length(list_mats1))){
    m = list_mats1[[i]][ids, ]
    if(length(mega_matrix)==0){
      mega_matrix = m
    }else{
      mega_matrix = cbind(mega_matrix, m)
    }
  }
  
  for(i in c(1:length(list_mats2))){
    m = list_mats2[[i]]
    #chain = "BCR"
    #if(i==2){chain = "TCR"}
    #names(m) = paste(chain, names(m))
    for(j in c(1:length(m))){
      m1 = m[[j]][gsub("_","-", ids, fixed = T),]
      colnames(m1) = apply(cbind(names(list_mats2)[i], names(m)[j], colnames(m1)), 1, paste, collapse = "|")
      if(length(grep("usage",names(list_mats2)[i] ))!=0){
        for(i1 in c(1:length(ids))){
          if(sum(m1[i1,])>=5){
            m1[i1,] = m1[i1,]*100/sum(m1[i1,])
          }else{m1[i1,] = NA}
        }
      }
      if(length(mega_matrix)==0){
        mega_matrix = m1
      }else{
        mega_matrix = cbind(mega_matrix, m1)
      }
    }
  }
  
  for(i in c(1:length(pathways))){
    m = pathways[[i]][gsub("_","-", ids, fixed = T), ]
    colnames(m) = apply(cbind(names(pathways)[i],colnames(m)), 1, paste, collapse = "|")
    mega_matrix = cbind(mega_matrix, m)
  }
  
  #### filter mega_matrix
  
  a = apply(mega_matrix, 2, function(x){length(which(is.na(x)==F))})
  mega_matrix1 = mega_matrix[,which(a>=19)]
  
  a = apply(mega_matrix1, 2, function(x){length(unique(x))})
  mega_matrix1 = mega_matrix1[,which(a>6)]
  
  exclude = colnames(mega_matrix1)[grep("NA", colnames(mega_matrix1))]
  exclude = exclude[grep("GO", exclude, invert = T)]
  exclude = exclude[grep("KEGG", exclude, invert = T)]
  exclude = exclude[grep("REACTOME", exclude, invert = T)]
  
  exclude = colnames(mega_matrix1)[grep("|-", colnames(mega_matrix1), fixed = T)]
  mega_matrix1 = mega_matrix1[, which(colnames(mega_matrix1) %in% exclude==F)]
  
  a = apply(mega_matrix1, 2, function(x){length(which(is.na(x)==F))})
  
  w1 = intersect(colnames(mega_matrix1)[grep("GO_B_CELL_RECEPTOR_SIGNALING_PATHWAY", colnames(mega_matrix1))], colnames(mega_matrix1)[grep("T cell", colnames(mega_matrix1))])
  w2 = intersect(colnames(mega_matrix1)[grep("GO_B_CELL_RECEPTOR_SIGNALING_PATHWAY", colnames(mega_matrix1))], colnames(mega_matrix1)[grep("NK", colnames(mega_matrix1))])
  w3 = intersect(colnames(mega_matrix1)[grep("GO_B_CELL_PROLIFERATION", colnames(mega_matrix1))], colnames(mega_matrix1)[grep("T cell", colnames(mega_matrix1))])
  w4 = intersect(colnames(mega_matrix1)[grep("GO_B_CELL_PROLIFERATION", colnames(mega_matrix1))], colnames(mega_matrix1)[grep("NK", colnames(mega_matrix1))])
  w5 = intersect(colnames(mega_matrix1)[grep("GO_T_CELL_PROLIFERATION", colnames(mega_matrix1))], colnames(mega_matrix1)[grep("B cell", colnames(mega_matrix1))])
  w6 = intersect(colnames(mega_matrix1)[grep("KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY", colnames(mega_matrix1))], colnames(mega_matrix1)[grep("T cell", colnames(mega_matrix1))])
  w7 = intersect(colnames(mega_matrix1)[grep("KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY", colnames(mega_matrix1))], colnames(mega_matrix1)[grep("NK", colnames(mega_matrix1))])
  w8 = intersect(colnames(mega_matrix1)[grep("KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", colnames(mega_matrix1))], colnames(mega_matrix1)[grep("B cell", colnames(mega_matrix1))])
  exclude = c(w1,w2,w3,w4,w5,w6,w7,w8)
  mega_matrix1 = mega_matrix1[, which(colnames(mega_matrix1) %in% exclude==F)]
  
  saveRDS(file = concat(c(output_directory, "VDJ_REMIX_megamatrix.rds")), mega_matrix1)
}

##### apply REMIX to mega_matrix1
library(vdjremix)
mega_matrix1= readRDS(file = concat(c(output_directory, "VDJ_REMIX_megamatrix.rds")))

res <- run_vdjremix(
  feature_matrix = mega_matrix1,
  na_frequencies = c(0.1, 0.2),
  min_cluster_size = 3,
  n_repeats = 20,
  n_cores = 8
)

#key outputs
modules = res$modules            # feature modules
mat_remix = res$eigengenes         # module eigengenes 
var_explained = res$variance_explained # variance explained by each module
loadings = res$loadings           # feature contributions within modules
rownames(mat_remix) = rownames(mega_matrix1)
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


### do same for blood

#### blood
ids_association = c(groups_PCA[[1]], groups_PCA[[2]])
ids_association = gsub("biopsy", "blood",ids_association)
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




##plot 
library(ggpubr)
library(ggplot2)

sample_colors <- c("Myeloid enriched PDAC" = "#E7298AFF", "Adaptive enriched PDAC" = "#17BECFFF")

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
    alpha = 0.2,          
    add = "jitter",       
    add.params = list(size = 1.2, alpha = 0.7),
    width = 0.5,          # Tighter box width
    palette = sample_colors
  ) +
    stat_compare_means(
      comparisons = list(c("Group 1", "Group 2")),
      method = "t.test",   
      label = "p.signif",
      bracket.size = 0.5
    ) +
    labs(
      title = mod_name,
      subtitle = paste0("ANOVA p = ", formatC(p_value[mod_name], format = "e", digits = 2)),
      x = "Patient Group",
      y = "Module Score"
    ) +
    theme_pubr() +   
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) +# Clean publication-ready theme
    theme(legend.position = "none") # Hide legend as X-axis labels suffice
  
}
