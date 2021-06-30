# UCell: Robust and scalable single-cell gene signature scoring


`UCell` is an R package for calculating gene signatures in single-cell datasets. UCell scores, based on the Mann-Whitney U statistic, are robust to dataset size and heterogeneity, and their calculation demands relatively less computing time and memory than other available methods, enabling the processing of large datasets (>10^5 cells). UCell can be applied to any cell vs. gene data matrix, and includes functions to directly interact with Seurat objects. 


![UCell_figure](https://github.com/carmonalab/UCell/blob/master/docs/Figure1.png?raw=true)


Find the installation instructions for the package and usage vignettes below.

### Package Installation

To install `UCell` directly from its GitHub repository, run the following code from within R or RStudio:
```
library(remotes)
remotes::install_github("carmonalab/UCell")
```

### Test the package

Load sample data and test your installation:
```
library(UCell)

my.matrix <- UCell::sample.matrix
gene.sets <- list(Tcell_signature = c("CD2","CD3E","CD3D"),
			Myeloid_signature = c("SPI1","FCER1G","CSF1R"))

scores <- ScoreSignatures_UCell(my.matrix, features=gene.sets)
head(scores)
```

### Examples and tutorials

Run UCell demos to learn about the functionalities of the package:

* [Single-cell gene signature scoring with UCell](https://carmonalab.github.io/UCell/UCell_matrix_vignette.html)

* [Using UCell with Seurat objects](https://carmonalab.github.io/UCell/UCell_Seurat_vignette.html)

* [Using UCell and Seurat to identify different T cell subtypes/states in human tumors](https://carmonalab.github.io/UCell/UCell_vignette_TILstates.html)

### New in version 1.1.0

You can now specify positive and negative (up- or down-regulated) genes in signatures. For example, build signatures as:

```
markers <- list()
markers$Tcell_gd <- c("TRDC+", "TRGC1+", "TRGC2+", "TRDV1+","TRAC-","TRBC1-","TRBC2-")
markers$Tcell_NK <- c("FGFBP2+", "SPON2+", "KLRF1+", "FCGR3A+", "CD3E-","CD3G-")
markers$Tcell_CD4 <- c("CD4","CD40LG")
markers$Tcell_CD8 <- c("CD8A","CD8B")
markers$Tcell_Treg <- c("FOXP3","IL2RA")
```

Genes without a + or - sign are assumed to be positive genes for that signature.

The UCell score is calculated as:
$$
U = max(0, U^+ - w_neg * U^-)
$$ 

where U^+ and U^- are respectively the U scores for the positive and negative set, and 'w_neg' is a weight on the negative set.

When no negative set of genes is present, U = U^+, therefore the behavior is identical to previous UCell versions.   

### Documentation

See a description of the functions implemented in UCell at: [UCell functions](docs/functions.md)


### Citation

UCell: robust and scalable single-cell gene signature scoring. Massimo Andreatta & Santiago J Carmona
https://doi.org/10.1101/2021.04.13.439670
