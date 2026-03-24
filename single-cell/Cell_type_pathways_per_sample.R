# srun -p short --cpus-per-task 8 --pty bash
##data availability: https://www.nature.com/articles/s41467-024-55424-2#data-availability 

##### run 
module purge
module load Anaconda3/2024.02-1
conda init
source ~/.bashrc
conda --version
conda activate seurat_env
R
#### 

library("Seurat")
library(ggplot2)
library(future) 

library(fgsea) ## not compatible when seurat is loaded!
library(ggplot2)
library(BiocParallel)
library(org.Hs.eg.db)

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


a = 1
if(a==1){
	file="Samples_PDAC150K_2.txt"
	p <- as.matrix(read.csv(file, head=T, sep="\t"))
	p=p[which(p[,"To_use_in_PDAC150K"]=="Yes"),]
	sample_id = as.character(p[,"Sample_Name"])
	sample_output_id = as.character(p[,"Sample_Name"])
	BCR.location = as.character(p[,"Location_of_BCR"])
	BCR.location = gsub("all_contig","filtered_contig",BCR.location)
	TCR.location = as.character(p[,"Location_of_TCR"])
	TCR.location = gsub("all_contig","filtered_contig",TCR.location)
	Overall_sample_group = as.character(p[,"Patient"])
	Site = as.character(p[,"Sample_type"])
	batch = "PDAC150Ka"
	PLOTS = "PLOTS/"
	out_dir = "/10X_GENOMICS/PDAK150K_WORKING_DATA/"
	out_dir_raw = "/10X_GENOMICS/PDAK150K_WORKING_DATA/"
	}
	

######################################

All_cell_apoptosis_markers_PancrImmune<-function(){
  ######## all cells: check apoptosis etc
  type = "all"
  analysis = "Pathways"
  pbmc = readRDS(file=concat(c("PDAC_RPCA_integrated_refined_annotation_scaled.rds")))
  
  file = concat(c("Overall_UMAP_annotations_PDAC150Ka_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  ids_x = rownames(p)
  inter = intersect(ids_x, rownames(pbmc@meta.data))
  length(inter)
  length(rownames(pbmc@meta.data))

  cell_type = pbmc@meta.data$cell_refined_annotation
  names(cell_type) = rownames(pbmc@meta.data)
  cell_types = sort(unique(cell_type))
  Patient = pbmc@meta.data$Patient
  Patient = gsub("_CD45p1","",Patient)
  Patient = gsub("_CD45p2","",Patient)
  Sample.Type = pbmc@meta.data$Sample.Type
  orig.ident  = apply(cbind(Patient, Sample.Type),1,paste,collapse = "-")
  names(orig.ident) = rownames(pbmc@meta.data)
  orig.idents = sort(unique( orig.ident))
  cell_ids = rownames(pbmc@meta.data)
  overall_cell_type = pbmc@meta.data$cell_refined_annotation
  names(overall_cell_type) = rownames(pbmc@meta.data)
  overall_cell_types = sort(unique( overall_cell_type))
  names(Patient) = cell_ids
  names(Sample.Type) = cell_ids

  pat_sample = apply(cbind(Patient, Sample.Type), 1, paste, collapse = "-")
  pat_samples = sort(unique(pat_sample))
  pbmc@meta.data$pat_sample = pat_sample
  
  library(fgsea) ## not compatible when seurat is loaded!
  pathways.hallmark <- gmtPathways("msigdb.v6.2.symbols.gmt.txt")
  
  ### get pathways to measure
  n = names(pathways.hallmark)
  pathways0 = c("KEGG_APOPTOSIS" , "REACTOME_APOPTOSIS","GO_EXECUTION_PHASE_OF_APOPTOSIS","GO_POSITIVE_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS",
               "REACTOME_EXTRINSIC_PATHWAY_FOR_APOPTOSIS","REACTOME_APOPTOSIS_INDUCED_DNA_FRAGMENTATION" ,"REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS" ,
               "HALLMARK_APOPTOSIS","GO_CELLULAR_COMPONENT_DISASSEMBLY_INVOLVED_IN_EXECUTION_PHASE_OF_APOPTOSIS")
  pathways1 = c("GO_T_CELL_PROLIFERATION","GO_B_CELL_PROLIFERATION" ,"GO_LEUKOCYTE_PROLIFERATION", "GO_CELL_PROLIFERATION" ,
                "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"    ,                 
                "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY" ,                    
                "GO_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY"         ,
                "GO_REGULATION_OF_B_CELL_RECEPTOR_SIGNALING_PATHWAY"         ,
                "GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY"                       ,"GO_B_CELL_RECEPTOR_SIGNALING_PATHWAY","HALLMARK_HYPOXIA")
  pathways2 = intersect(n[grep("HYPOX", n)], n[grep("GO_", n)])
  pathways3 = intersect(n[grep("INTERFERON", n)], n[grep("GO_", n)])
  
  pathways_all =c(pathways0, pathways1, pathways2, pathways3)
  
  pathways.hallmark1 = pathways.hallmark[which(names(pathways.hallmark) %in% pathways_all)]
  scale_data= pbmc@assays$ RNA@ scale.data
  genes = rownames(scale_data)
  
  for(g in c(1:length(pathways.hallmark1))){
    gene_list = pathways.hallmark1[[g]]
    gene_list = intersect(genes, gene_list)
    if(length(gene_list)>=3){
      pbmc <- AddModuleScore(object = pbmc,features = list(gene_list),ctrl = 5,name = concat(c(names(pathways.hallmark1)[g])),assay = 'RNA',slot = 'scale.data')
      print(g)
    }
  }
  names_pathways = names(pathways.hallmark1)
  list_mean_per_cell_type_pathway= NULL
  
  for(g in c(1:length(names_pathways))){
    mean_per_cell_type_pathway = matrix(data = NA, nrow = length(pat_samples), ncol = length(overall_cell_types), dimnames = c(list(pat_samples), list(overall_cell_types)))
    score = pbmc@meta.data[,concat(c(names_pathways[g],"1"))]
    for(s in c(1:length(pat_samples))){
      for(c in c(1:length(overall_cell_types))){
        w = intersect(which (pat_sample==pat_samples[s]), which(overall_cell_type==overall_cell_types[c]))
        if(length(w)>=5){
          mean_per_cell_type_pathway[pat_samples[s], overall_cell_types[c]]= median(score[w])
        }
      }
    }
    list_mean_per_cell_type_pathway = c(list_mean_per_cell_type_pathway, list(mean_per_cell_type_pathway))
  }
  names(list_mean_per_cell_type_pathway) = names_pathways

  concat(c(out_dir,PLOTS,"/",analysis,"_mean_pathway_score_per_cell_type_per_sample_", batch,"_raw_pathway_level.rds"))
  
  saveRDS(file = concat(c(out_dir,PLOTS,"/",analysis,"_mean_pathway_score_per_cell_type_per_sample_", batch,"_raw_pathway_level.rds")), list_mean_per_cell_type_pathway)
  
}
  
  
  
  
