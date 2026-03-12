setwd("~/OneDrive - Nexus365/DPhil Project/Ageing/AUTOIMMUNE")

##data cleaning
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
out_dir = "~/Library/CloudStorage/OneDrive-Nexus365/DPhil Project/Ageing/AUTOIMMUNE/"
batch = "Autoimm_rWGNCA"
##############
file = concat(c(out_dir, "Eigenvectors_AI_BCR.txt"))
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
mat = p
heatmap(p)

###load age and other metadata
file = concat(c(out_dir,"Overall_autoimmune paper summary locations_share(Table S3).csv"))
p1 <- as.matrix(read.csv(file, sep=","))
df1 = data.frame(p1)
library(dplyr)
df1 <- df1 %>%
  mutate(
    age_numeric = as.numeric(`Patient.age.at.time.of.sampling`), 
    age_group = case_when(
      between(age_numeric, 20, 40) ~ "20-40",
      between(age_numeric, 41, 60) ~ "41-60", 
      between(age_numeric, 61, 80) ~ "61-80", 
      age_numeric >= 81 ~ "81+",              
      TRUE ~ NA_character_ 
    )
  )

table(df1$age_group)

df1 <- df1[1:238,]
rownames(df1) <- df1$Sample
p <- p[ order(row.names(p)), ]
df1 <- df1[ order(row.names(df1)), ]

# setdiff(rownames(df1), rownames(p))  
# [1] "ISO5_9_F4"
df1 <- df1[rownames(df1) %in% rownames(p), , drop = FALSE]
## check 
rownames(df1) == rownames(p)

# --- Efficient Data Loading and Merging ---
setwd("~/OneDrive - Nexus365/DPhil Project/Ageing/AUTOIMMUNE")

##data cleaning
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
out_dir = "~/Library/CloudStorage/OneDrive-Nexus365/DPhil Project/Ageing/AUTOIMMUNE/"
batch = "Autoimm_rWGNCA"

# 1. Load and prepare the primary metadata (df1) and eigenvector (p) data
df1 <- read.csv(concat(c(out_dir,"Overall_autoimmune paper summary locations_share(Table S3).csv")),
                stringsAsFactors = FALSE) # Better to load as strings initially
p <- as.matrix(read.csv(concat(c(out_dir, "Eigenvectors_AI_BCR.txt")), head=TRUE, sep="\t"))

# Clean up sample names and ensure order matches before merging
df1 <- df1[1:238,]
rownames(df1) <- df1$Sample
p <- p[ order(row.names(p)), ]
df1 <- df1[ order(row.names(df1)), ]
df1 <- df1[rownames(df1) %in% rownames(p), , drop = FALSE]

# Convert 'p' to a data frame and join with 'df1' using a common ID
p_df <- as.data.frame(p) %>%
  tibble::rownames_to_column(var = "row_id")

df1 <- df1 %>%
  tibble::rownames_to_column(var = "row_id") %>%
  left_join(p_df, by = "row_id")


# 2. Load and prepare the MMF therapy data into its OWN data frame
tp_mmf <- read.delim("~/Downloads/Time_since_MMF.txt")
tp_mmf <- tp_mmf %>%
  rename(row_id = samples) # Rename 'samples' to 'row_id' for a clean join

# 3. Load and prepare the RTX therapy data into its OWN data frame
tp_rtx <- read.delim("~/Downloads/Time_since_RTX.txt")
tp_rtx <- tp_rtx %>%
  rename(row_id = samples) # Rename 'samples' to 'row_id' for a clean join


# 4. Sequentially join the therapy data to the main df1 data frame
# The first join adds the MMF columns.
# The second join adds the RTX columns to the SAME data frame.
df1 <- df1 %>%
  left_join(tp_mmf, by = "row_id") %>%
  left_join(tp_rtx, by = "row_id")

# 5. Verify the result (optional but recommended)
# You should now see columns like 'time_since_MMFs' AND 'time_since_RTXs' in the same data frame
glimpse(df1)
# Make sure all required libraries are loaded
library(dplyr)
library(ggplot2)
library(ggpubr)
# library(grid) # No longer needed as grid annotations are removed

# --- 1. Data Preparation ---

df1 <- df1 %>%
  mutate(disease_broad = case_when(
    grepl("AAV|EGPA", Disease.group) ~ "AAV",
    grepl("SLE", Disease.group) ~ "SLE",
    grepl("Healthy", Disease.group) ~ "Healthy",
    TRUE ~ "Other"
  ))

# --- Define Constants ---
TARGET_MODULE <- "Module_9"
MODULE_LABEL <- "Module 9"

disease_colors <- c("AAV" = "#228B22", "SLE" = "#FFB347")
healthy_color <- "#F49595"

# --- 2. Create the Five Data Subsets for Plotting ---
df_aav_diagnosis <- df1 %>%
  filter(disease_broad == "AAV" & grepl("At diagnosis|At flare", Disease.and.flare.status))
df_sle_diagnosis <- df1 %>%
  filter(disease_broad == "SLE" & grepl("At diagnosis|At flare", Disease.and.flare.status))
df_mmf <- df1 %>%
  filter(therapies.x == "MMF", disease_broad %in% c("AAV", "SLE"),
         !is.na(time_since_MMFs), time_since_MMFs <= 500, time_since_MMFs >= 0)
df_rtx <- df1 %>%
  filter(therapies.y == "RTX", disease_broad %in% c("AAV", "SLE"),
         !is.na(time_since_RTXs), time_since_RTXs <= 500, time_since_RTXs >= 0)
df_healthy <- df1 %>%
  filter(disease_broad == "Healthy")

# --- 3. Determine Common Y-axis Range ---
y_values <- c(
  df_aav_diagnosis[[TARGET_MODULE]], df_sle_diagnosis[[TARGET_MODULE]],
  df_mmf[[TARGET_MODULE]], df_rtx[[TARGET_MODULE]], df_healthy[[TARGET_MODULE]]
)
# <<< CHANGE: Removed the extra padding on the y-axis
y_range <- range(y_values, na.rm = TRUE)

# --- 4. Create the Five Individual Plots ---
# (Plotting code is unchanged, but now uses the tighter y_range)
p_aav_diag <- ggplot(df_aav_diagnosis, aes(x = "AAV Pre-Treatment\n(day 0)", y = .data[[TARGET_MODULE]])) +
  geom_boxplot(fill = disease_colors["AAV"], width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(shape = 21, fill = "white", width = 0.15, size = 2) + coord_cartesian(ylim = y_range) +
  labs(x = NULL, y = paste(MODULE_LABEL, "Score")) + theme_classic2() + theme(axis.text.x = element_text(size = 9))
p_sle_diag <- ggplot(df_sle_diagnosis, aes(x = "SLE Pre-Treatment\n(day 0)", y = .data[[TARGET_MODULE]])) +
  geom_boxplot(fill = disease_colors["SLE"], width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(shape = 21, fill = "white", width = 0.15, size = 2) + coord_cartesian(ylim = y_range) +
  labs(x = NULL, y = "") + theme_classic2() + theme(axis.text.x = element_text(size = 9), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

# <<< CHANGE: Adjusted cor.coef.coord to work with the un-padded y-axis
p_mmf <- ggscatter(df_mmf, x = "time_since_MMFs", y = TARGET_MODULE, cor.coef = TRUE, cor.coef.coord = c(50, y_range[2]*0.9), cor.method = "pearson", xlab = "Time since MMF therapy (d)", ylab = "", color = "disease_broad", palette = disease_colors, add = "reg.line", conf.int = FALSE, add.params = list(color = "grey30"), shape = "disease_broad", ggtheme = theme_classic2()) +
  coord_cartesian(ylim = y_range) + theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
p_rtx <- ggscatter(df_rtx, x = "time_since_RTXs", y = TARGET_MODULE, cor.coef = TRUE, cor.coef.coord = c(50, y_range[2]*0.9), cor.method = "pearson", xlab = "Time since RTX therapy (d)", ylab = "", color = "disease_broad", palette = disease_colors, add = "reg.line", conf.int = FALSE, add.params = list(color = "grey30"), shape = "disease_broad", ggtheme = theme_classic2()) +
  coord_cartesian(ylim = y_range) + theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
p_healthy <- ggplot(df_healthy, aes(x = "Healthy", y = .data[[TARGET_MODULE]])) +
  geom_boxplot(fill = healthy_color, width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(shape = 21, fill = "white", width = 0.15, size = 2) + coord_cartesian(ylim = y_range) +
  labs(x = NULL, y = "") + theme_classic2() + theme(axis.text.x = element_text(size = 9), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

# --- 5. P-value Calculation (Results printed to Console) ---
diag_group <- df1 %>%
  filter(disease_broad %in% c("AAV", "SLE"), grepl("At diagnosis|At flare", Disease.and.flare.status))
early_mmf_group <- df1 %>%
  filter(disease_broad %in% c("AAV", "SLE"), !is.na(time_since_MMFs) & time_since_MMFs > 0 & time_since_MMFs <= 90)
early_rtx_group <- df1 %>%
  filter(disease_broad %in% c("AAV", "SLE"), !is.na(time_since_RTXs) & time_since_RTXs > 0 & time_since_RTXs <= 90)

# Perform and Print Test 1: Diagnosis vs Early MMF
test_mmf <- wilcox.test(diag_group[[TARGET_MODULE]], early_mmf_group[[TARGET_MODULE]])
message("\n--- Statistical Comparison: Diagnosis vs. Early MMF ---")
message(paste("N (Diagnosis/Flare):", nrow(diag_group)))
message(paste("N (1-3 Months Post-MMF):", nrow(early_mmf_group)))
message(paste("Wilcoxon p-value:", format.pval(test_mmf$p.value, digits = 3)))
message("------------------------------------------------------")

# Perform and Print Test 2: Diagnosis vs Early RTX
test_rtx <- wilcox.test(diag_group[[TARGET_MODULE]], early_rtx_group[[TARGET_MODULE]])
message("\n--- Statistical Comparison: Diagnosis vs. Early RTX ---")
message(paste("N (Diagnosis/Flare):", nrow(diag_group)))
message(paste("N (1-3 Months Post-RTX):", nrow(early_rtx_group)))
message(paste("Wilcoxon p-value:", format.pval(test_rtx$p.value, digits = 3)))
message("------------------------------------------------------\n")


# --- 6. Arrange Plots ---
# This is now the final plot object, without the extra annotations.
final_plot <- ggarrange(
  p_aav_diag, p_sle_diag, p_mmf, p_rtx, p_healthy,
  ncol = 5, nrow = 1, align = "h",
  widths = c(1.5, 1.5, 2.5, 2.5, 1.5),
  common.legend = TRUE, legend = "bottom"
)

# Print the final plot to the RStudio plot pane
print(final_plot)


pdf("Module_9_dynamics.pdf", 10, 6)
final_plot
dev.off()

# Heatmap
library(dplyr)
library(tidyr)
library(ggplot2)

module_cols_original <- paste0("Module_", 1:20)
module_labels <- paste0("Module ", 1:20)

if (is.factor(df1$Disease.and.flare.status)) {
  df1$Disease.and.flare.status <- as.character(df1$Disease.and.flare.status)
}
df1$Disease.and.flare.status[df1$Disease.and.flare.status == "NA"] <- NA

# Subset all data to at diagnosis
df1 <- df1 %>% filter(grepl("At diagnosis|Healthy", df1$Disease.and.flare.status, ignore.case = TRUE))

healthy_data <- df1 %>% filter(Disease.and.flare.status == "Healthy" & !is.na(Disease.and.flare.status))
all_unique_diseases <- unique(df1$Disease.and.flare.status)
diseases_to_compare <- setdiff(all_unique_diseases, "Healthy")
diseases_to_compare <- diseases_to_compare[!is.na(diseases_to_compare)]

# --- STATS FOR COLORS: HEALTH vs DISEASE (Your existing code) ---
comparison_results <- list()
for (idx in 1:length(module_cols_original)) {
  mod_col_original <- module_cols_original[idx]
  mod_label <- module_labels[idx]
  healthy_values <- healthy_data[[mod_col_original]]
  healthy_values <- healthy_values[!is.na(healthy_values)]
  for (disease_group in diseases_to_compare) {
    if (length(healthy_values) < 2) {
      p_val <- 1; direction <- "neutral"; mean_disease <- NA
    } else {
      disease_data <- df1 %>% filter(Disease.and.flare.status == disease_group & !is.na(Disease.and.flare.status))
      disease_values <- disease_data[[mod_col_original]]
      disease_values <- disease_values[!is.na(disease_values)]
      mean_disease <- mean(disease_values, na.rm = TRUE)
      if (length(disease_values) >= 2) {
        test_result <- wilcox.test(disease_values, healthy_values, exact = FALSE, conf.int = FALSE)
        p_val <- test_result$p.value
        direction <- if (mean(disease_values) > mean(healthy_values)) "higher" else "lower"
      } else {
        p_val <- 1; direction <- "neutral"
      }
    }
    comparison_results[[length(comparison_results) + 1]] <- data.frame(
      module = mod_label, disease_group = disease_group, p_value = p_val,
      mean_disease = mean_disease, mean_healthy = mean(healthy_values),
      direction = direction, stringsAsFactors = FALSE
    )
  }
}
results_df <- bind_rows(comparison_results)
results_df$FDR <- p.adjust(results_df$p_value, method = "BH")
results_df$category <- "Not significant"
results_df$category[results_df$FDR < 0.05 & results_df$direction == "higher"] <- "Significantly higher\nFDR < 0.05"
results_df$category[results_df$FDR < 0.005 & results_df$direction == "higher"] <- "Significantly higher\nFDR < 0.005"
results_df$category[results_df$FDR < 0.05 & results_df$direction == "lower"] <- "Significantly lower\nFDR < 0.05"
results_df$category[results_df$FDR < 0.005 & results_df$direction == "lower"] <- "Significantly lower\nFDR < 0.005"


# --- NEW LOGIC: STATS FOR ASTERISK - KRUSKAL-WALLIS AMONG DISEASES ---

# 1. Isolate disease-only data and prepare it for testing
disease_only_df <- df1 %>%
  filter(Disease.and.flare.status != "Healthy" & !is.na(Disease.and.flare.status)) %>%
  # Pivot to long format for easy formula-based testing (value ~ group)
  pivot_longer(
    cols = all_of(module_cols_original),
    names_to = "module_original",
    values_to = "score"
  )

# 2. Loop through each module, perform Kruskal-Wallis test, and collect p-values
kw_results_list <- list()
for (mod_col in module_cols_original) {
  
  # Filter data for the current module
  module_data_for_test <- disease_only_df %>%
    filter(module_original == mod_col) %>%
    # Drop levels to ensure kruskal.test doesn't see empty disease groups
    mutate(Disease.and.flare.status = factor(Disease.and.flare.status))
  
  # Perform Kruskal-Wallis test if there are at least 2 disease groups with data
  if (length(unique(module_data_for_test$Disease.and.flare.status)) >= 2) {
    kw_test_result <- tryCatch({
      kruskal.test(score ~ Disease.and.flare.status, data = module_data_for_test)
    }, error = function(e) { NULL })
    
    p_value <- if (!is.null(kw_test_result)) kw_test_result$p.value else 1
  } else {
    p_value <- 1 # Not enough groups to compare
  }
  
  kw_results_list[[mod_col]] <- data.frame(module_original = mod_col, p_value = p_value)
}

# 3. Combine results and apply FDR correction
kw_results_df <- bind_rows(kw_results_list)
kw_results_df$FDR_across_diseases <- p.adjust(kw_results_df$p_value, method = "BH")

# 4. Generate the final asterisk data frame for plotting
asterisk_col_name <- " " # Dummy column name for plotting
asterisk_data <- kw_results_df %>%
  # Map original module names to the plot labels
  mutate(module = gsub("_", " ", module_original)) %>%
  # Create the label: "*" if significant, "" otherwise
  mutate(label = if_else(FDR_across_diseases < 0.05, "*", "")) %>%
  # Assign it to our new dummy column name
  mutate(disease_group = asterisk_col_name) %>%
  # Make sure the 'module' factor levels match the main plot for correct y-axis alignment
  mutate(module = factor(module, levels = rev(module_labels))) %>%
  select(module, disease_group, label)


# --- PLOTTING PREP & EXECUTION (Modified to use the new asterisk logic) ---
# ... (rest of your plotting prep code) ...
category_levels_ordered <- c("Significantly higher\nFDR < 0.005", "Significantly higher\nFDR < 0.05", "Significantly lower\nFDR < 0.05", "Significantly lower\nFDR < 0.005", "Not significant")
category_colors <- c("Significantly higher\nFDR < 0.005" = "#8B0000", "Significantly higher\nFDR < 0.05" = "#F8BB96", "Significantly lower\nFDR < 0.005" = "#00008B", "Significantly lower\nFDR < 0.05" = "#8EC4E5", "Not significant" = "#FFFFFF")

results_df$category <- factor(results_df$category, levels = category_levels_ordered)
results_df$module <- factor(results_df$module, levels = rev(module_labels))

desired_disease_order_from_image <- c("Behcets At diagnosis", "SLE At diagnosis", "AAV MPO+ At diagnosis", "AAV PR3+ At diagnosis", "EGPA ANCA- At diagnosis", "EGPA MPO+ At diagnosis", "CD At diagnosis", "IgAV At diagnosis")
final_disease_groups_in_plot <- intersect(desired_disease_order_from_image, unique(results_df$disease_group))
if(length(final_disease_groups_in_plot) == 0) { final_disease_groups_in_plot <- sort(unique(results_df$disease_group)) }

all_x_levels <- c(final_disease_groups_in_plot, asterisk_col_name)
results_df$disease_group <- factor(results_df$disease_group, levels = all_x_levels)

heatmap_plot <- ggplot(results_df, aes(x = disease_group, y = module)) +
  geom_tile(aes(fill = category), color = "grey50", linewidth = 0.38) +
  # This geom_text now uses the asterisk_data from the Kruskal-Wallis tests
  geom_text(data = asterisk_data, aes(label = label), size = 8, color = "black", fontface = "bold", vjust = 0.75) +
  scale_fill_manual(values = category_colors, name = "Significance compared to health", labels = function(breaks) gsub("\n", " ", breaks), breaks = category_levels_ordered, drop = FALSE) +
  scale_x_discrete(drop = FALSE, labels = c(final_disease_groups_in_plot, "")) +
  labs(x = "Diseases and Flare Status", y = "Module Scores", title = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10, colour = "black"),
    axis.text.y = element_text(size=10, colour = "black"),
    axis.title.x = element_text(face="bold", size=10, margin = margin(t = 10)),
    axis.title.y = element_text(face="bold", size=10, margin = margin(r = 10)),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "grey80"),
    legend.position = "right",
    legend.title = element_text(face="bold", size=11),
    legend.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

pdf("heatmap_significance_with_inter-disease_asterisks.pdf", 6.5, 5)
print(heatmap_plot)
dev.off()

print(heatmap_plot)
