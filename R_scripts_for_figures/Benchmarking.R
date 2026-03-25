# df = feature matrix (samples × features; imputed)
# 
# p = VDJ-REMIX eigengenes (samples × K)
# 
# df1 = metadata with age_numeric, age_group
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

rownames(df1) <- df1$Sample
p <- p[ order(row.names(p)), ]
df1 <- df1[ order(row.names(df1)), ]

df1 <- df1[rownames(df1) %in% rownames(p), , drop = FALSE]
## check 
rownames(df1) == rownames(p)

df <- read.delim("Imputed_DATA_FINAL_AI_BCR.txt") 

p <- p[ order(row.names(p)), ]
df <- df[ order(row.names(df)), ]

df <- df[rownames(df) %in% rownames(p), , drop = FALSE]
stopifnot(all(rownames(df) == rownames(p)))
stopifnot(all(rownames(df1) == rownames(p)))

set.seed(1234)

K <- ncol(p)

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()  # optional

datExpr <- as.data.frame(df)  # WGCNA likes data.frame

# Basic QC
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Choose soft-threshold power
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)

# Pick power (you can automate, but start with the recommended)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6  # fallback

# Build network & modules
wgcna <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  corType = "pearson",
  minModuleSize = 10,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

wgcna_colors <- labels2colors(wgcna$colors)

# Module eigengenes
MEs <- moduleEigengenes(datExpr, colors = wgcna_colors)$eigengenes
MEs <- orderMEs(MEs)

# Align to your sample order
MEs <- MEs[rownames(p), , drop = FALSE]



####MOFA
#BiocManager::install("MOFA2")
library(MOFA2)

X <- t(as.matrix(df))  # features × samples
data_list <- list(AIRR = X)

MOFAobject <- create_mofa(data_list)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- K
train_opts <- get_default_training_options(MOFAobject)
train_opts$maxiter <- 1000
train_opts$convergence_mode <- "fast"
train_opts$seed <- 1234
train_opts$verbose <- TRUE
train_opts$stochastic <- FALSE
train_opts$gpu_mode <- FALSE

# Option 1: set likelihoods directly during prepare_mofa (works in many versions)
MOFAobject <- prepare_mofa(
  MOFAobject,
  model_options = model_opts,
  training_options = train_opts
)


MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)
Z <- get_factors(MOFAobject)
str(Z)
names(Z)
F_mofa <- Z$group1
dim(F_mofa)
head(F_mofa[,1:3])





E_remix <- p

# WGCNA might have >K modules; use top-K by variance
var_me <- apply(MEs, 2, var)
E_wgcna <- MEs[, names(sort(var_me, decreasing = TRUE))[1:min(K, ncol(MEs))], drop = FALSE]

# MOFA already trained to K
E_mofa <- F_mofa


###bench
W_mofa <- get_weights(MOFAobject)$AIRR  # features × factors

loading_entropy <- function(w) {
  p <- abs(w) / sum(abs(w))
  -sum(p * log(p + 1e-12))
}

ent_mofa <- apply(W_mofa, 2, loading_entropy)

df_ent <- data.frame(Factor = seq_along(ent_mofa), Entropy = ent_mofa)
library(ggplot2)
ggplot(df_ent, aes(x = Factor, y = Entropy)) +
  geom_point() +
  labs(y = "Loading entropy (lower = more concentrated)", x = "MOFA factor")

kME <- signedKME(datExpr, MEs)


##therapy
y_raw <- df1$Disease.and.flare.status

y_bin <- ifelse(
  grepl("^SLE", y_raw), "SLE",
  ifelse(grepl("^AAV", y_raw), "AAV", NA)
)

y_bin <- factor(y_bin)
table(y_bin, useNA = "ifany")

keep <- !is.na(y_bin)

y <- y_bin[keep]

E_remix <- E_remix[keep, , drop = FALSE]
E_wgcna <- E_wgcna[keep, , drop = FALSE]
E_mofa  <- E_mofa[keep,  , drop = FALSE]

# Alignment check
stopifnot(
  all(rownames(E_remix) == rownames(df1)[keep]),
  all(rownames(E_wgcna) == rownames(df1)[keep]),
  all(rownames(E_mofa)  == rownames(df1)[keep])
)
cv_binary_accuracy <- function(E, y, folds = 5, repeats = 10, seed = 123) {
  set.seed(seed)
  y <- factor(y)
  n <- nrow(E)
  
  accs <- numeric(folds * repeats)
  idx <- 1
  
  for (r in seq_len(repeats)) {
    fold_id <- sample(rep(1:folds, length.out = n))
    
    for (k in seq_len(folds)) {
      tr <- which(fold_id != k)
      te <- which(fold_id == k)
      
      x_tr <- as.matrix(E[tr, , drop = FALSE])
      x_te <- as.matrix(E[te, , drop = FALSE])
      y_tr <- y[tr]
      y_te <- y[te]
      
      # Skip folds without both classes
      if (length(unique(y_tr)) < 2) {
        accs[idx] <- NA
        idx <- idx + 1
        next
      }
      
      fit <- glmnet::cv.glmnet(
        x_tr, y_tr,
        family = "binomial",
        type.measure = "class"
      )
      
      # ---- EXPLICIT glmnet predict (avoid MOFA masking) ----
      pred_mat <- glmnet:::predict.glmnet(
        fit$glmnet.fit,
        x_te
      )  # samples × lambdas
      
      # Find lambda.min column
      lambda_idx <- which.min(abs(fit$glmnet.fit$lambda - fit$lambda.min))
      
      prob <- as.numeric(pred_mat[, lambda_idx])
      
      pred <- ifelse(prob > 0.5, levels(y_tr)[2], levels(y_tr)[1])
      
      accs[idx] <- mean(pred == y_te, na.rm = TRUE)
      idx <- idx + 1
    }
  }
  
  c(
    mean = mean(accs, na.rm = TRUE),
    sd   = sd(accs,   na.rm = TRUE),
    n    = sum(!is.na(accs))
  )
}
res_remix <- cv_binary_accuracy(E_remix, y)
res_wgcna <- cv_binary_accuracy(E_wgcna, y)
res_mofa  <- cv_binary_accuracy(E_mofa,  y)

bench <- data.frame(
  Method  = c("VDJ-REMIX", "WGCNA", "MOFA2"),
  MeanAcc = c(res_remix["mean"], res_wgcna["mean"], res_mofa["mean"]),
  SD      = c(res_remix["sd"],   res_wgcna["sd"],   res_mofa["sd"]),
  Nfolds  = c(res_remix["n"],    res_wgcna["n"],    res_mofa["n"])
)

bench


###entropy
entropy_topN_norm <- function(w, N = 50, eps = 1e-12) {
  w <- abs(w)
  w <- w[is.finite(w)]
  if (length(w) == 0) return(NA_real_)
  
  # take top-N (or fewer if not enough)
  w <- sort(w, decreasing = TRUE)
  N_use <- min(N, length(w))
  w <- w[seq_len(N_use)]
  
  s <- sum(w)
  if (s == 0) return(NA_real_)
  
  p <- w / s
  
  H <- -sum(p * log(p + eps))
  H_norm <- H / log(N_use)  # in [0,1]
  return(H_norm)
}


vdj_loadings_list <- readRDS("~/Downloads/Module_FeaturePCALoadings_AI_184_BCR.rds")
entropy_vdj <- sapply(vdj_loadings_list, loading_entropy)

df_entropy_vdj <- data.frame(
  Method = "VDJ-REMIX",
  Dimension = names(entropy_vdj),
  Entropy = entropy_vdj
)

# datExpr: samples × features
# MEs: samples × modules
kME_mat <- WGCNA::signedKME(datExpr, MEs)  # features × modules (usually)

entropy_wgcna <- apply(kME_mat, 2, entropy_topN_norm, N = 50)

df_entropy_wgcna <- data.frame(
  Method = "WGCNA",
  Dimension = colnames(kME_mat),
  Entropy = entropy_wgcna
)


W_mofa <- get_weights(MOFAobject)$AIRR  

entropy_mofa <- apply(W_mofa, 2, loading_entropy)

df_entropy_mofa <- data.frame(
  Method = "MOFA2",
  Dimension = colnames(W_mofa),
  Entropy = entropy_mofa
)
df_entropy <- rbind(
  df_entropy_vdj,
  df_entropy_wgcna,
  df_entropy_mofa
)

library(ggplot2)

ggplot(df_entropy, aes(x = Method, y = Entropy, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  theme_classic() +
  ylab("Loading entropy (lower = more interpretable)") +
  xlab(NULL)




##variance explained
# --- sanity checks ---
stopifnot(is.matrix(X) || is.data.frame(X))
X <- as.matrix(X)

stopifnot(is.matrix(E_remix), is.matrix(E_wgcna), is.matrix(E_mofa))
stopifnot(identical(rownames(df), rownames(E_remix)))
stopifnot(identical(rownames(df), rownames(E_wgcna)))
stopifnot(identical(rownames(df), rownames(E_mofa)))

# --- helper: weighted mean R2 across features for ONE latent dimension ---
weighted_mean_r2 <- function(z, X, eps = 1e-12) {
  z <- as.numeric(z)
  ok <- is.finite(z)
  if (sum(ok) < 3) return(NA_real_)
  z <- z[ok]
  Xok <- X[ok, , drop = FALSE]
  
  vj <- apply(Xok, 2, var, na.rm = TRUE)
  vj[!is.finite(vj)] <- 0
  keep_feat <- vj > eps
  if (sum(keep_feat) < 2) return(NA_real_)
  vj <- vj[keep_feat]
  Xok <- Xok[, keep_feat, drop = FALSE]
  
  zc <- z - mean(z, na.rm = TRUE)
  denom_z <- sum(zc^2)
  if (!is.finite(denom_z) || denom_z < eps) return(NA_real_)
  
  Xc <- scale(Xok, center = TRUE, scale = FALSE)
  denom_x <- sqrt(colSums(Xc^2))
  denom_x[denom_x < eps] <- NA_real_
  
  cor_xz <- as.numeric(crossprod(Xc, zc)) / (denom_x * sqrt(denom_z))
  r2 <- cor_xz^2
  r2[!is.finite(r2)] <- NA_real_
  
  keep <- is.finite(r2) & is.finite(vj) & (vj > 0)
  if (sum(keep) < 2) return(NA_real_)
  sum(vj[keep] * r2[keep]) / sum(vj[keep])
}

# --- SAFE per-dimension variance explained (no apply()) ---
variance_profile <- function(E, X) {
  stopifnot(is.matrix(E))
  out <- vapply(seq_len(ncol(E)),
                function(j) weighted_mean_r2(E[, j], X = X),
                numeric(1))
  names(out) <- colnames(E)
  out
}

K <- ncol(E_remix)

# WGCNA: keep top-K by eigengene variance if more than K
if (ncol(E_wgcna) > K) {
  var_me <- apply(E_wgcna, 2, var)
  E_wgcna_use <- E_wgcna[, names(sort(var_me, decreasing = TRUE))[1:K], drop = FALSE]
} else E_wgcna_use <- E_wgcna

# MOFA: keep first K if more
if (ncol(E_mofa) > K) {
  E_mofa_use <- E_mofa[, 1:K, drop = FALSE]
} else E_mofa_use <- E_mofa

# --- compute ---
r2_remix <- variance_profile(E_remix, X)
r2_wgcna <- variance_profile(E_wgcna_use, X)
r2_mofa  <- variance_profile(E_mofa_use, X)

df_var <- rbind(
  data.frame(Method="VDJ-REMIX", Dimension=names(r2_remix), R2=as.numeric(r2_remix)),
  data.frame(Method="WGCNA",     Dimension=names(r2_wgcna), R2=as.numeric(r2_wgcna)),
  data.frame(Method="MOFA2",     Dimension=names(r2_mofa),  R2=as.numeric(r2_mofa))
)

# --- plot ---
library(ggplot2)

p_var <- ggplot(df_var, aes(x=Method, y=R2, fill=Method)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) +
  geom_jitter(width=0.15, size=1.2, alpha=0.6) +
  theme_classic() +
  ylab(expression(paste("Global variance explained (variance-weighted mean ", R^2, " across features)"))) +
  xlab(NULL)

print(p_var)

# --- stats: Kruskal-Wallis + Dunn (Holm) ---
library(FSA)
kw <- kruskal.test(R2 ~ Method, data=df_var)
kw

dt <- FSA::dunnTest(R2 ~ Method, data=df_var, method="holm")
dt$res











