Seurat tutorial
================
Georgia Miller
2024-10-24

## This is following the Seurat tutorial from <https://satijalab.org/seurat/articles/pbmc3k_tutorial>

This is to analysing a dataset of Peripheral Blood Mononuclear Cells
(PBMCs) available from 10X genomics. 2,700 single cells were sequenced
on Illumina NextSeq500.

# Read in data

Read10X reads in output of the cellranger pipeline form 10X, returns
Unique Molecular Identified (UMI) count matrix. Values in matrix
represents number of molecules for each feature (gene, row) detected in
each cell (column). More recent versions use h5 file format, read in
using Read10X_h5() function.

``` r
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/k2477939/Documents/myrepo/seurat_tutorial/filtered_gene_bc_matrices/hg19/")
```

## Set up Seurat object

Use count matrix to create a Seurat object which is a container for both
data (e.g. count matrix) and analysis (PCA, clustering results) for
sc-dataset.

``` r
# Initialise Seurat object with raw non-normalised data. 
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

``` r
pbmc
```

    ## An object of class Seurat 
    ## 13714 features across 2700 samples within 1 assay 
    ## Active assay: RNA (13714 features, 0 variable features)
    ##  1 layer present: counts

## Examine data

Examine count matrix data. . values are 0s, this sparse-matrix
representation saves memory and speed

``` r
# Examine some genes in the first 30 cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
```

    ## 3 x 30 sparse Matrix of class "dgCMatrix"

    ##   [[ suppressing 30 column names 'AAACATACAACCAC-1', 'AAACATTGAGCTAC-1', 'AAACATTGATCAGC-1' ... ]]

    ##                                                                    
    ## CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
    ## TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
    ## MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .

``` r
# Sparse size
object.size(as.matrix(pbmc.data))
```

    ## 709591472 bytes

``` r
# Dense size
object.size(pbmc.data)
```

    ## 29905192 bytes

# Pre-processing

## QC and selecting cells for further analysis

Common QC metrics include

- No. of unique genes detected in each cell
  - low-quality cells or empty droplets often have very few genes
  - cell doublets or multiplets may have aberrantly high gene count
- total no. molecules detected in a cell
  - correlates strongly with unique genes
- percentage of reads that map to the mitochondrial contamination
  - low-quality/dying cells often have extensive mt contamination
  - calculates percentage of counts originates from a set of features
    with PercentageFeatureSet()
  - set of all genes starting with MT- used as set of mt genes

Find no. unique genes and total molecules in the object metadata as they
are automatically calculated in CreateSeuratObject()

``` r
 # For first 5 cells
head(pbmc@meta.data, 5)
```

    ##                  orig.ident nCount_RNA nFeature_RNA
    ## AAACATACAACCAC-1     pbmc3k       2419          779
    ## AAACATTGAGCTAC-1     pbmc3k       4903         1352
    ## AAACATTGATCAGC-1     pbmc3k       3147         1129
    ## AAACCGTGCTTCCG-1     pbmc3k       2639          960
    ## AAACCGTGTATGCG-1     pbmc3k        980          521

Use PercentageFeatureSet() to find mt contamination. Add column to
object metadata using \[\[

``` r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

## Visualise QC metrics

- Violin plot

``` r
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha = 0.5)
```

    ## Warning: Default search for "data" layer in "RNA" assay yielded no results;
    ## utilizing "counts" layer instead.

![](seurat_tutorial_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

- Feature scatter

``` r
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

![](seurat_tutorial_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Filter cells

Filter cells based on QC metrics e.g. here keep cells with 200-2,500
unique features and with \<5% mitochondrial counts.

``` r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & nFeature_RNA < 2500 & percent.mt < 5)
```

# Normalise the data

By default, use a global-scaling normalisation method “LogNormalize”
that normalises feature expression measurements for each cell by total
expression, multiplies by a scale factor (default 10,000) and
log-transforms the result. Normalised values are stored in
pbmc\[\[“RNA”\]\]\$data in Seurat v5.

``` r
# Extra arguments here are not needed as using default
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

    ## Normalizing layer: counts

Assumes each cell originally contains the same number of RNA molecules.
Alternative workflows do not assume this e.g. SCTransform()

# Feature selection

Next, we calculate a subset of features that exhibit high cell-to-cell
variation in in the dataset. Have found that focusing on these
downstream helps to highlight biological signal in sc datasets.

FindVariableFeatures() function has improved on previous versions by
directly modelling the mean-variance relationship inherent in sc data.

- Variance-stabilising transformation to correct for the mean-variance
  relationship inherent to sc RNA-seq. If only choose genes based on
  log-norm sc variance would not account for this
  - mean and variance of each gene in unnorm data, apply
    lg10-transformation, fit a curve to predict variance of each gene as
    a function of its mean (local fitting of polynomials of degree 2).
    Provides a regularised estimator of variance given then mean of a
    feature so can use it to standardise feature counts without removing
    higher-than-expected variation.
- Transformation:
  - where F = feature, V = value, sd = standard deviation, C = cell
    standardised V of Fi in Cj = (raw V of Fi in Cj - mean raw V of
    Fi)/expected sd of Fi from global mean variance fit
  - clipped standardised values to max root of total number of cells
  - found variance of standardised values across all cells, this
    variance is a measure of sc dispersion after controlling for mean
    exp, used to directly rank the features
  - by default, returns 2,000 features per dataset with the highest
    standardised variance
- If multiple datasets integrated
  - feature selection on individual datasets as above
  - give priority to features that are highly variable in multiple
    experiments
  - take top 2,000 for downstream analysis
  - broke ties by examining ranks of tied features in each original
    dataset and take those with the highest median rank
    SelectIntegrationFeatures()

``` r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
```

    ## Finding variable features for layer counts

``` r
# Identify 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
```

    ## When using repel, set xnudge and ynudge to 0 for optimal results

``` r
plot1 + plot2
```

    ## Warning in scale_x_log10(): log-10 transformation introduced infinite values.
    ## log-10 transformation introduced infinite values.

![](seurat_tutorial_files/figure-gfm/unnamed-chunk-10-1.png)<!-- --> \#
Scaling

Apply a linear transformation (scaling) is a standard pre-processing
step to dimensionality reduction tehcnqiues

ScaleData() does:

- shifts exp of each gene so mean across all is 0 and variance is 1
  - gives equal weight to genes in downstream analyses so highly
    expressed genes do not dominate
- stores results in pbmc\[\[“RNA”\]\]\$scale.data
- only variable features are scaled by default
- specify features argument to scale additional features

``` r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

    ## Centering and scaling data matrix

Can also remove unwanted sources of variation from a dataset. E.g. could
‘regress out’ heterogeneity associated with e.g. cell cycle stage or mt
contamination

But, recommend SCTransform() for advanced users

``` r
#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
```

# Linear dimensional reduction

Perform PCA on scaled data. By default, use previously determined
variable features but can define different subset using features (must
pass through ScaleData first)

For first PCs, outputs a list of genes with most positive and most
negative loadings, representing modules of genes with correation or
anti-correlation across scs

``` r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

    ## PC_ 1 
    ## Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP 
    ##     FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP 
    ##     PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD 
    ## Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW 
    ##     CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A 
    ##     MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3 
    ## PC_ 2 
    ## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74 
    ##     HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB 
    ##     BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74 
    ## Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2 
    ##     CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX 
    ##     TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC 
    ## PC_ 3 
    ## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1, HLA-DPA1, CD74, MS4A1, HLA-DRB1, HLA-DRA 
    ##     HLA-DRB5, HLA-DQA2, TCL1A, LINC00926, HLA-DMB, HLA-DMA, CD37, HVCN1, FCRLA, IRF8 
    ##     PLAC8, BLNK, MALAT1, SMIM14, PLD4, P2RX5, IGLL5, LAT2, SWAP70, FCGR2B 
    ## Negative:  PPBP, PF4, SDPR, SPARC, GNG11, NRGN, GP9, RGS18, TUBB1, CLU 
    ##     HIST1H2AC, AP001189.4, ITGA2B, CD9, TMEM40, PTCRA, CA2, ACRBP, MMD, TREML1 
    ##     NGFRAP1, F13A1, SEPT5, RUFY1, TSC22D1, MPP1, CMTM5, RP11-367G6.3, MYL9, GP1BA 
    ## PC_ 4 
    ## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1, CD74, HIST1H2AC, HLA-DPB1, PF4, SDPR 
    ##     TCL1A, HLA-DRB1, HLA-DPA1, HLA-DQA2, PPBP, HLA-DRA, LINC00926, GNG11, SPARC, HLA-DRB5 
    ##     GP9, AP001189.4, CA2, PTCRA, CD9, NRGN, RGS18, CLU, TUBB1, GZMB 
    ## Negative:  VIM, IL7R, S100A6, IL32, S100A8, S100A4, GIMAP7, S100A10, S100A9, MAL 
    ##     AQP3, CD2, CD14, FYB, LGALS2, GIMAP4, ANXA1, CD27, FCN1, RBP7 
    ##     LYZ, S100A11, GIMAP5, MS4A6A, S100A12, FOLR3, TRABD2A, AIF1, IL8, IFI6 
    ## PC_ 5 
    ## Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY, CCL4, CST7, PRF1, GZMA, SPON2 
    ##     GZMH, S100A9, LGALS2, CCL3, CTSW, XCL2, CD14, CLIC3, S100A12, RBP7 
    ##     CCL5, MS4A6A, GSTP1, FOLR3, IGFBP7, TYROBP, TTC38, AKR1C3, XCL1, HOPX 
    ## Negative:  LTB, IL7R, CKB, VIM, MS4A7, AQP3, CYTIP, RP11-290F20.3, SIGLEC10, HMOX1 
    ##     LILRB2, PTGES3, MAL, CD27, HN1, CD2, GDI2, CORO1B, ANXA5, TUBA1B 
    ##     FAM110A, ATP1A1, TRADD, PPA1, CCDC109B, ABRACL, CTD-2006K23.1, WARS, VMO1, FYB
