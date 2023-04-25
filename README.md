# scMEB
R package for the scMEB methods.

This package provides a novel and fast method for detecting single-cell DEGs 
without prior cell clustering results. The proposed method utilizes a small 
part of known non-DEGs (stably expressed genes) to build a minimum enclosing 
ball and defines the DEGs based on the distance of a mapped gene to the center 
of the hypersphere in a feature space. The method on this package is described 
in the article 'scMEB: A fast and clustering-independent method for detecting 
differentially expressed genes in single-cell RNA-seq data' [1]. 

# Installation
This package can be installed through `devtools` in R:
```{r}
library("devtools")
devtools::install_github("FocusPaka/scMEB")
```

# User's Guide
Please refer to the 
[vignetee](https://github.com/FocusPaka/scMEB/blob/main/vignettes/scMEB.Rmd) 
for detailed function instructions.

# Reference
1. Zhu, J.D, Yang, Y.L. scMEB: A fast and clustering-independent method for 
detecting differentially expressed genes in single-cell RNA-seq data. 
(2023, pending publication)



