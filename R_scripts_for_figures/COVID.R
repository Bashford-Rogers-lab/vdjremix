library(mgcv)
library(dplyr)
set.seed(1974)

#get the matrix with broad meta data 
# meta = read.delim("Clinical data_COMBAT.txt")
# df <- read.delim("Eigenvectors_COVID_BCR.txt")
# df$Generic.name <- gsub("^PROG_(.*)_productive$", "\\1", rownames(df))
# 
# df1 <- left_join(new, df, by="Generic.name")
# df1$INDIVIDUAL <- as.integer(df1$INDIVIDUAL)
# df2 <- left_join(df1, broad_meta, by="INDIVIDUAL")
# df2$AGE <- as.numeric(df2$AGE)
# 
# completed <- df2[complete.cases(df2$AGE),]
# completed <- completed[complete.cases(completed$Module_1),]
# 
# #####basic model
# pdf("Module_scores_with_age_bigger.pdf")
# par(cex.lab = 1, cex.axis = 1.5)
# for (i in colnames(completed)[grep("Module", colnames(completed))]) {
#   mod.ss <- smooth.spline(completed$AGE, completed[[i]], all.knots = T)
#   plot(completed$AGE, completed[[i]], xlab = "age", ylab = "Module score")
#   lines(mod.ss, col = "blue")
# }
# dev.off()




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
out_dir = "~/OneDrive - Nexus365/DPhil Project/Ageing/COMBAT/"
batch = "COV_rWGNCA"
##############
file = concat(c("~/OneDrive - Nexus365/DPhil Project/Ageing/COMBAT/Eigenvectors_BCR_COVID_BCR_with_uncorrelated_features.txt"))
#file = concat(c("~/OneDrive - Nexus365/DPhil Project/Ageing/COMBAT/Imputed_DATA_FINAL_COVID_BCR.txt"))
mat <- as.data.frame(read.csv(file, head=TRUE, sep="\t", row.names = 1))
mat <- mat[1:78,]
#mat <- mat[1:(nrow(mat)-5), ]
mat1 <- as.matrix(sapply(mat, function(x) as.numeric(as.character(x))))
rownames(mat1) <- rownames(mat)
mat = mat1

mat <- mat[ order(row.names(mat)), ]
p = mat
heatmap(p)

###load age and other metadata
library(dplyr)
mat <- as.data.frame(mat)
broad_meta = read.delim("~/Library/CloudStorage/OneDrive-Nexus365/DPhil Project/Ageing/COMBAT/COMBAT_BCR_IDs2.txt")
broad_meta$Sequencing.ID

mat$Sequencing.ID = rownames(mat)
full <- left_join(mat, broad_meta, by = "Sequencing.ID")
rownames(full) <- full$Sequencing.ID
full <- full[ order(row.names(full)), ]
p <- as.data.frame(p)

full <- full %>% mutate(age_group = as.numeric(as.character(Age)),
                      age_group = case_when(
                        between(age_group, 20, 40.9) ~ "20-40",
                        between(age_group, 41, 60.9) ~ "40-60",
                        between(age_group, 61, 80.9) ~ "60-80",
                        age_group >= 81 ~ "80+",
                        TRUE ~ NA_character_))


full$covid_status <- ifelse(full$Source %in% c("COVID_HCW_MILD", "COVID_MILD"), "mild",
                            ifelse(full$Source %in% c("COVID_SEV", "COVID_CRIT"), "crit",
                                   full$Source))
table(full$covid_status)

rownames(df1) <- rownames(p)
rownames(df2) <- rownames(p)


#recovery test
library(dplyr)
library(clinfun) 
library(readr)   

df1 = full
names(df1)[30:36] <- paste0("Module_", names(df1)[30:36])

ordered_groups <- c("COVID_MILD","COVID_HCW_MILD", "HV" )

df_trend <- df1 %>%
  filter(Source %in% ordered_groups) %>%
  mutate(Source = factor(Source, levels = ordered_groups, ordered = TRUE))

module_cols <- colnames(df_trend)[grepl("Module_", colnames(df_trend))]


# -----------------------------------------------------------------
trend_results <- list() 

for (mod_col in module_cols) {
  
  temp_df <- df_trend %>%
    dplyr::select(Source, Score = all_of(mod_col)) %>% 
    filter(!is.na(Score))
  
   if (n_distinct(temp_df$Source) < 2) {
    message(paste("Skipping module", mod_col, "- data present in fewer than 2 groups after removing NAs."))
    next 
  }
  
  jt_test_result <- tryCatch({
    jonckheere.test(x = temp_df$Score, g = temp_df$Source, alternative = "decreasing")
  }, error = function(e) {
    warning(paste("Jonckheere test failed for module", mod_col, "with error:", e$message))
    NULL
  })
  
  if (!is.null(jt_test_result)) {
    trend_results[[mod_col]] <- data.frame(
      module = mod_col,
      p_value = jt_test_result$p.value,
      statistic_JT = jt_test_result$statistic
    )
  }
}

# 3. Combine Results and Perform Multiple Testing Correction
# ----------------------------------------------------------
trend_results_df <- bind_rows(trend_results)

if(nrow(trend_results_df) > 0) {
  trend_results_df$FDR <- p.adjust(trend_results_df$p_value, method = "BH")
  
  # Print the final results, ordered by significance
  print("Trend Analysis Results:")
  print(trend_results_df %>% arrange(FDR))
  
  # Print just the significant ones
  print("Modules with a significant increasing trend (FDR < 0.05):")
  print(trend_results_df %>% filter(FDR < 0.05) %>% arrange(FDR))
  
} else {
  print("No results were generated. This might mean no module had sufficient data across groups to be tested.")
}


# 4. Visualization (Example for a top module)
# -------------------------------------------
# After running the above, let's assume 'Module_X' was your top hit
# Replace "Module_X" with a real module name from your significant results.
# If no modules were significant, pick one with a low p-value to visualize the trend.

if (nrow(trend_results_df) > 0) {
  library(ggplot2)
  
  # Let's take the top module from the results
  module_to_plot <- trend_results_df %>% arrange(FDR) %>% pull(module) %>% .[1]
  
  # Get the FDR value to display in the title
  fdr_value <- trend_results_df %>% filter(module == module_to_plot) %>% pull(FDR)
  
  trend_plot <- ggplot(df_trend, aes(x = Source, y = !!sym(module_to_plot), fill = Source)) +
    geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.7, shape = 21, color = "black") +
    # Add the linear trend line for visual confirmation
    geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_fill_manual(values = c("HV" = "#8dd3c7", "COVID_MILD" = "#ffffb3", "COVID_HCW_MILD" = "#bebada")) + 
    labs(
      title = paste("", gsub("_", " ", module_to_plot)),
      subtitle = paste("Jonckheere-Terpstra Test for trend, FDR =", format.pval(fdr_value, digits = 3)),
      x = "Group",
      y = "Module Score"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
  
  print(trend_plot)
}

###publication figure
pdf("trend_test.pdf", 4, 6)

if (nrow(trend_results_df) > 0) {
  library(ggplot2)
  
  # Filter for modules with a significant FDR value
  significant_modules <- trend_results_df %>% 
    filter(FDR < 0.05) %>% 
    pull(module)
  
  # Loop through each significant module to create a plot
  for (module_to_plot in significant_modules) {
    
    # Get the FDR value for the current module
    fdr_value <- trend_results_df %>% 
      filter(module == module_to_plot) %>% 
      pull(FDR)
    
    trend_plot <- ggplot(df_trend, aes(x = Source, y = !!sym(module_to_plot), fill = Source)) +
      geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.15, height = 0, alpha = 0.7, shape = 21, color = "black") +
      # Add the linear trend line for visual confirmation
      geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black", linetype = "dashed", linewidth = 0.8) +
      scale_fill_manual(values = c("HV" = "#8dd3c7", "COVID_MILD" = "#ffffb3", "COVID_HCW_MILD" = "#bebada")) +
      labs(
        title = paste("", gsub("_", " ", module_to_plot)),
        subtitle = paste("Jonckheere-Terpstra Test for trend, FDR =", format.pval(fdr_value, digits = 3)),
        x = "Group",
        y = "Module Score"
      ) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 10, color = "black", hjust = 1, angle = 45, vjust = 0.5)
      )
    
    print(trend_plot)
  }
}

dev.off()
print(trend_results_df %>% filter(FDR < 0.05) %>% arrange(FDR))

                         
##### PLS-DA
## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', out.width = '90%')


# --- 1. Setup: Load Libraries and Prepare Data ---

library(dplyr)
library(tidyr)
library(mixOmics)
library(ggplot2)

# --- NEW: Define additional variables and create a clean mapping for labels ---
additional_vars <- c(
  "Percentage_max_cluster_size..IGHG1",
  "V_gene_replacement_clonal_expansion..d5_norm",
  "J_gene_freq_by_uniq_VDJ_IGHD.IGHM_unmutated..IGHJ1",
  "J_gene_freq_by_uniq_VDJ_IGHE..IGHJ4",
  "J_gene_freq_by_uniq_VDJ_IGHE..IGHJ6",
  "J_gene_freq_by_uniq_VDJ_IGHG4..IGHJ2",
  "J_gene_freq_by_uniq_VDJ_IGHG4..IGHJ4"
)

# Create cleaner labels for these variables for plotting later
additional_var_labels <- c(
  "Max_Clust_Size_IGHG1",
  "V_Gene_Repl_d5_norm",
  "J_Freq_IGHD_M_unmut_IGHJ1",
  "J_Freq_IGHE_IGHJ4",
  "J_Freq_IGHE_IGHJ6",
  "J_Freq_IGHG4_IGHJ2",
  "J_Freq_IGHG4_IGHJ4"
)

# Create a named vector for easy renaming later
label_map <- setNames(additional_var_labels, additional_vars)

# --- Combine original modules with the new variables ---
module_cols_original <- paste0("Module_", 1:29)
all_predictors <- c(module_cols_original, additional_vars)

###filter
#remove - HCW_MILD, group critical and sev together, keep mild, hv, sepsis? (try both), remove batch control


set.seed(5249) # for reproducibility, remove for normal use
full$covid_status <- ifelse(full$Source %in% c("COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"), "Covid",
                            ifelse(full$Source %in% c("HV"), "health",
                                   full$Source))
##test1 - sepsis vs health
df1 <- full %>%
  dplyr::filter(Source != "COVID_HCW_MILD", Source != "Batch control") %>%
  dplyr::filter(covid_status %in% c("Sepsis", "health"))


#test 2 - covid vs health
df1 <- full %>%
  dplyr::filter(Source != "COVID_HCW_MILD", Source != "Batch control") %>%
  dplyr::filter(covid_status %in% c("Covid", "health"))


full$covid_status <- ifelse(full$Source %in% c("COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"), "cov",
                            ifelse(full$Source %in% c("HV"), "health",
                                   full$Source))


df1 <- full %>%
  dplyr::filter(Source != "COVID_HCW_MILD", Source != "Batch control") %>%
  dplyr::filter(covid_status %in% c("sev", "COVID_MILD", "HV"))

df1 <- full %>%
  dplyr::filter(Source != "COVID_HCW_MILD", Source != "Batch control") %>%
  dplyr::filter(covid_status %in% c("sev", "HV"))

full$covid_status <- ifelse(full$Source %in% c("COVID_SEV", "COVID_CRIT"), "sev",
                                   full$Source)

df1 <- full %>%
  dplyr::filter(Source != "COVID_HCW_MILD", Source != "Batch control") %>%
  dplyr::filter(covid_status %in% c("cov", "health"))

outcome_var <- "covid_status"

# --- Create a clean data matrix with the EXPANDED predictor set ---
plsda_data <- df1 %>%
  dplyr::select(all_of(outcome_var), all_of(all_predictors)) %>%
  tidyr::drop_na()

# Create the final X and Y for mixOmics
X <- plsda_data %>%
  dplyr::select(all_of(all_predictors)) %>%
  as.matrix()

Y <- plsda_data %>%
  pull(.data[[outcome_var]]) %>%
  as.factor() %>%
  droplevels()

dim(X) # check the dimensions of the X dataframe
summary(Y) # check the distribution of class labels


## ---- fig.cap = "FIGURE 1: Barplot of the variance each principal component explains of the SRBCT gene expression data."----
pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE) # run pca method on data
plot(pca.srbct)  # barplot of the eigenvalues (explained variance per component)


## ---- fig.cap = "FIGURE 2: Preliminary (unsupervised) analysis with PCA on the SRBCT gene expression data"----------
plotIndiv(pca.srbct, group = plsda_data$covid_status, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') # onto the PCA subspace


## -------------------------------------------------------------------------------------------------------------------
srbct.splsda <- mixOmics::splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later


## ---- fig.show = "hold", out.width = "49%", fig.cap = "FIGURE 3: Sample plots of the SRBCT gene expression data after a basic PLS-DA model was operated on this data. (a) depicts the samples with the confidence ellipses of different class labels while (b) depicts the prediction background generated by these samples. Both plots use the first two components as axes."----
# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2,
          group = plsda_data$covid_status, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = plsda_data$covid_status, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")


## ---- fig.cap = "FIGURE 4: Tuning the number of components in PLS-DA on the SRBCT gene expression data. For each component, repeated cross-validation (10 × 3−fold CV) is used to evaluate the PLS-DA classification performance (OER and BER), for each type of prediction distance; `max.dist`, `centroids.dist` and `mahalanobis.dist`."----
# undergo performance evaluation in order to tune the number of components to use
perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold",
                          folds = 5, nrepeat = 100, # use repeated cross-validation
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")


## -------------------------------------------------------------------------------------------------------------------
perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()


## ---- fig.cap = "FIGURE 5:  Tuning keepX for the sPLS-DA performed on the SRBCT gene expression data. Each coloured line represents the balanced error rate (y-axis) per component across all tested keepX values (x-axis) with the standard deviation based on the repeated cross-validation folds. As sPLS-DA is an iterative algorithm, values represented for a given component (e.g. comp 1 to 2) include the optimal keepX value chosen for the previous component (comp 1)."----
# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 100, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 2) # allow for paralleliation to decrease runtime


plot(tune.splsda.srbct, col = color.jet(4)) # plot output of variable number tuning


## -------------------------------------------------------------------------------------------------------------------
tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()



## -------------------------------------------------------------------------------------------------------------------
tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()


## -------------------------------------------------------------------------------------------------------------------
optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]


## -------------------------------------------------------------------------------------------------------------------
# form final model with optimised values for component and variable count
final.splsda <- mixOmics::splsda(X, Y,
                       ncomp = optimal.ncomp,
                       keepX = optimal.keepX)


## ---- fig.show = "hold", out.width = "49%", fig.cap = "FIGURE 6:  Sample plots from sPLS-DA performed on the SRBCT gene expression data including 95% confidence ellipses. Samples are projected into the space spanned by the first three components. (a) Components 1 and 2 and (b) Components 1 and 3. Samples are coloured by their tumour subtypes."----
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = plsda_data$covid_status, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = plsda_data$covid_status, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## # set the styling of the legend to be homogeneous with previous plots
## legend=list(legend = levels(Y), # set of classes
##             col = unique(color.mixo(Y)), # set of colours
##             title = "Tumour Type", # legend title
##             cex = 0.7) # legend size
##
## # generate the CIM, using the legend and colouring rows by each sample's class
## cim <- cim(final.splsda, row.sideColors = color.mixo(Y),
##            legend = legend)


## ---- fig.cap = "FIGURE 8:  Stability of variable selection from the sPLS-DA on the SRBCT gene expression data. The barplot represents the frequency of selection across repeated CV folds for each selected gene for component 1 (a), 2 (b) and 3 (c)."----
# form new perf() object which utilises the final model
perf.splsda.srbct <- perf(final.splsda,
                          folds = 10, nrepeat = 100, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)


# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.srbct$features$stable[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2)
plot(perf.splsda.srbct$features$stable[[2]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2)
plot(perf.splsda.srbct$features$stable[[3]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(c) Comp 3', las =2)


## ---- fig.cap = "FIGURE 9: Correlation circle plot representing the genes selected by sPLS-DA performed on the SRBCT gene expression data. Gene names are truncated to the first 10 characters. Only the genes selected by sPLS-DA are shown in components 1 and 2."----
var.name.short <- substr(X[, 2], 1, 10) # form simplified gene names

plotVar(final.splsda, comp = c(1,2), var.names = list(var.name.short), cex = 3) # generate correlation circle plot


## -------------------------------------------------------------------------------------------------------------------
train <- sample(1:nrow(X), 35) # randomly select 50 samples in training
test <- setdiff(1:nrow(X), train) # rest is part of the test set

# store matrices into training and test set:
X.train <- X[train, ]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]



## -------------------------------------------------------------------------------------------------------------------
# train the model
train.splsda.srbct <- mixOmics::splsda(X.train, Y.train, ncomp = optimal.ncomp, keepX = optimal.keepX)


## -------------------------------------------------------------------------------------------------------------------
# use the model on the Xtest set
predict.splsda.srbct <- predict(train.splsda.srbct, X.test, dist = "mahalanobis.dist")



## -------------------------------------------------------------------------------------------------------------------
# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda.srbct$covid_status$mahalanobis.dist[,2]

table(factor(predict.comp2, levels = levels(Y)), Y.test)


## -------------------------------------------------------------------------------------------------------------------
# evaluate the prediction accuracy for the first three components
predict.comp3 <- predict.splsda.srbct$covid_status$mahalanobis.dist[,3]
table(factor(predict.comp3, levels = levels(Y)), Y.test)



## ---- fig.show = "hold", out.width = "49%", fig.cap = "FIGURE 10: ROC curve and AUC from sPLS-DA on the SRBCT gene expression data on component 1 (a) and all three components (b) averaged across one-vs.-all comparisons."----
auc.splsda = auroc(final.splsda, roc.comp = 1, print = T) # AUROC for the first component



auc.splsda = auroc(final.splsda, roc.comp = 2, print = T)# AUROC for all three components



##plot for publication
plot_object <- plotIndiv(final.splsda, comp = c(1,2),
                         group = plsda_data$covid_status,
                         ind.names = FALSE,
                         ellipse = TRUE,
                         legend = TRUE,
                         title = '')

custom_colors <- c("health" = "#FFA500", "Sepsis" = "#C9F")
custom_colors <- c("health" = "#FFA500", "Covid" = "#56ae57")
# 2. Build the final plot with corrected aesthetics
publication_plot <- ggplot(data = plot_object$df,
                           aes(x = x, y = y, fill = group, color = group)) +

  # Add the semi-transparent filled ellipses
  stat_ellipse(geom = "polygon", alpha = 0.4) +

  # Add the filled circles (points) with a darker border
  geom_point(shape = 21, size = 3, stroke = 1) +

  # Manually set the colors for both fill and color aesthetics
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +

  # Use a classic theme for a clean look
  theme_classic() +

  # Add labels and customize font sizes
  labs(x = "PLS1", y = "PLS2") +
  theme(
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_blank(), # Removes the legend title
    legend.text = element_text(size = 12)
  )

# 3. Print the final plot
print(publication_plot)

pdf("health_vs_covid.pdf", 5, 4)
print(publication_plot)
dev.off()


pdf("health_vs_sepsis.pdf", 5, 4)
print(publication_plot)
dev.off()


                         
