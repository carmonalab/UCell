# UCell - Functions

* `ScoreSignatures_UCell`    Calculate module enrichment scores from single-cell data. Given a gene vs. cell matrix, calculates module/signature enrichment scores on single-cell level using Mann-Whitney U statistic. Returned scores are normalized U statistic (equivalent to AUC - Area Under the Curve). This score depends only on the gene expression ranks of individual cell, and therefore is robust across datasets.

* `ScoreSignatures_UCell_sce`    A wrapper for UCell to interact directly with SingleCellExperiment objects. Give a `sce` object, UCell scores for the specified gene signatures are returned in altExp: `altExp(sce, 'UCell')`

* `AddModuleScore_UCell`   A wrapper for UCell to interact directly with Seurat objects. Given a Seurat object, calculates module/signature enrichment scores on single-cell level using Mann-Whitney U statistic. Returned scores are normalized U statistic (equivalent to AUC - Area Under the Curve). In contrast to Seurat's `AddModuleScore` (based on population average gene expression binning) this score depend only on the gene expression ranks of individual cell, and therefore is robust across datasets.

* `StoreRankings_UCell`   Calculate and store gene rankings for a single-cell dataset. Given a gene vs. cell matrix, calculates the rankings of expression for all genes in each cell. It can then be applied to the function `ScoreSignatures_UCell` to evaluate gene signatures on the gene expression ranks of individual cells.   

Find more information, syntax and examples using the R help function e.g. `?ScoreSignatures_UCell`

