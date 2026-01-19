# VDJ-REMIX

**VDJ-REMIX** (Repertoire Module Identification and eXploration) is an R package for the **unsupervised, interpretable modular analysis of immune receptor repertoire (AIRR-seq) feature matrices**.

It is designed to address the unique statistical challenges of repertoire summary data, including **high dimensionality, heterogeneous feature types, collinearity, and missingness**, which are not well handled by standard dimensionality reduction approaches.

VDJ-REMIX implements a **network-based modularisation framework**, inspired by Weighted Gene Correlation Network Analysis (WGCNA), but **refactored and optimised specifically for BCR and TCR repertoire feature matrices**.

---

## Key features

- **Tailored to AIRR-seq summary features**  
  Supports repertoire metrics such as V/J gene usage, isotype composition, somatic hypermutation (SHM), clonal diversity, and expansion indices.

- **Robust correlation estimation under missingness**  
  Pairwise correlations are estimated using repeated subsampling of shared non-missing observations, avoiding bias introduced by fully imputed matrices.

- **Automated, interpretable module discovery**  
  Features are grouped into biologically coherent modules.

- **Interpretable module summarisation**  
  Modules are designed to be sparse, internally coherent, and directly linkable to immunological mechanisms.

---

## Installation

VDJ-REMIX is under active development and can be installed from GitHub:

```r
# install.packages("devtools")
devtools::install_github("sakinaamin/vdjremix")
library(vdjremix)
data(example_airr_features)
dim(example_airr_features)


#use case - The entire VDJ-REMIX workflow can be executed using a single wrapper function, while retaining access to all intermediate results
res <- run_vdjremix(
  feature_matrix = example_airr_features,
  na_frequencies = c(0.1, 0.2),
  min_cluster_size = 5,
  n_repeats = 500,
  n_cores = 2
)


#key outputs
res$modules            # feature modules
res$eigengenes         # module eigengenes 
res$variance_explained # variance explained by each module
res$loadings           # feature contributions within modules


#For full worked example:
browseVignettes("vdjremix")
```

Contact

For questions or feedback, please contact:

Sakina (sakina.amin@sjc.ox.ac.uk), 
Lauren (lauren_overend@live.co.uk) or
Rachael (rachael.bashford-rogers@bioch.ox.ac.uk)
