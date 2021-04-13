# Robust and scalable single-cell gene signature enrichment analysis


`UCell` is an R package for calculating gene signatures in single-cell datasets. UCell enrichment scores, based on the Mann-Whitney U statistic, are robust to dataset size and heterogeneity, and their calculation demands relatively less computing time and memory than other available methods, enabling the processing of large datasets (>10^5 cells). UCell can be applied to any cell vs. gene data matrix, and includes functions to directly interact with Seurat objects. 

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

* [Signature enrichment analysis with UCell](https://carmonalab.github.io/UCell/UCell_matrix_vignette.html)

* [Using UCell with Seurat objects](https://carmonalab.github.io/UCell/UCell_Seurat_vignette.html)


![UCell_figure](https://github.com/carmonalab/UCell/blob/master/docs/Figure1.png?raw=true)

### Documentation

See a description of the functions implemented in UCell at: [UCell functions](docs/functions.md)


### Reference

Coming soon.
