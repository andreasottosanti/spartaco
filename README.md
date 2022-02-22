
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SPAtially Resolved TrAnscriptomics CO-clustering (SpaRTaCo)

by A. Sottosanti, D. Righelli and D. Risso

<!-- badges: start -->

<!-- badges: end -->

`spartaco` implements a novel statistical approach for detecting spatially expressed genes in sub-groups of spots in spatial transcriptomic experiments. Details of the statistical model appear in the article *Co-clustering of Spatially Resolved Transcriptomic Data* by A. Sottosanti and D. Risso available [here](https://arxiv.org/abs/2110.04872).

## Installation instructions

Install the development version from
[GitHub](https://github.com/andreasottosanti/spartaco) with:

``` r
BiocManager::install("andreasottosanti/spartaco")
```

## Run the model

Let `x` be the spatial experiment matrix containing the expression of `nrow(x)` genes measured over `ncol(x)` spots. The spatial coordinates of the spots are contained into the matrix `coordinates`. You can run SpaRTaCo searching for `K` clusters of genes and `R` clusters of spots by running the following code:

``` r
library(spartaco)
spartaco(x = x, coordinates = coordinates, K = K, R = R) 
```

### Starting points
The estimation algorithm is automatically run using random starting points. However, it is also possible to start from a previous output of the `spartaco` function:

```r 
output1 <- spartaco(x = x, coordinates = coordinates, K = K, R = R, max.iter = 10^3)
output2 <- spartaco(x = x, coordinates = coordinates, K = K, R = R, max.iter = 10^3, input.values = output1)
```
