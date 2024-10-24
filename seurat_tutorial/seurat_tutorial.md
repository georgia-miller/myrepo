seurat_tutorial
================
Georgia Miller
2024-10-24

This is following the Seurat tutorial from
<https://satijalab.org/seurat/articles/pbmc3k_tutorial>

This is to analysing a dataset of Peripheral Blood Mononuclear Cells
(PBMCs) available from 10X genomics. 2,700 single cells were sequenced
on Illumina NextSeq500.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 'SeuratObject' was built under R 4.4.0 but the current version is
    ## 4.4.1; it is recomended that you reinstall 'SeuratObject' as the ABI
    ## for R may have changed

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, t

``` r
library(patchwork)
```

Read10X reads in output of the cellranger pipeline form 10X, returns
Unique Molecular Identified (UMI) count matrix. Values in matrix
represents number of molecules for each feature (gene, row) detected in
each cell (column). More recent versions use h5 file format, read in
using Read10X_h5() function.

``` r
# load the PBMC dataset
pbmc_data <- Read10X(data.dir = "/Users/k2477939/Documents/myrepo/seurat_tutorial/filtered_gene_bc_matrices/hg19/")
```
