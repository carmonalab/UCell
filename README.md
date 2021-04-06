# Robust and scalable single-cell gene signature enrichment analysis


`UCell` is an R package for calculating gene signatures in single-cell datasets. UCell enrichment scores, based on the Mann-Whitney U statistic, are robust to dataset size and heterogeneity, and their calculation demands relatively less computing time and memory than other available methods, enabling the processing of large datasets (>10^5 cells). UCell can be applied to any cell vs. gene data matrix, and includes functions to directly interact with Seurat objects. 


Find the installation instructions for the package below, and a vignette detailing its functions at [Tutorial (html)](https://carmonalab.github.io/UCell/tutorial.html)

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


```
