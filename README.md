
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FindIT2:Find influential TF and influential Target

<!-- badges: start -->
<!-- badges: end -->

FindIT2 is a package that implements functions to find influential TF
and target based on different input type.It has five module

-   Multi-peak multi-gene annotaion(mmPeakAnno module)
-   Calculate regulation potential(calcRP module)
-   Find influential Target based on ChIP-Seq and RNA-Seq data(Find
    influential Target module)
-   Find influential TF based on different input(Find influential TF
    module)
-   Calculate peak-gene or peak-peak correlation(peakGeneCor module)

And there are also some other useful functions like integrate different
source information, calculate jaccard similarity for your TF. I will
introduce all these function in vignettes.

## Installation instructions

`FindIT2` is available on
[Bioconductor](https://bioconductor.org/packages/devel/bioc/html/FindIT2.html),
you can install it by:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("FindIT2")
```

> For the packages on bioconductor, please make sure you download the
> latest stable `R` release from [CRAN](http://cran.r-project.org/)

If you want the development version, install it directly from
[GitHub](https://github.com/shangguandong1996/FindIT2):

``` r
BiocManager::install("shangguandong1996/FindIT2")
```

## Document

If you want to download development version and view vignettes using
`browseVignettes(FindIT2)`, your R version should be 4.0 or greater
according to this
[issue](https://github.com/Bioconductor/BiocStyle/issues/78) because I
use the BiocStyle to output rmarkdown. Then you have to firstly download
the below packages

``` r
packages <- c("BiocStyle", "knitr", "rmarkdown", "sessioninfo", "TxDb.Athaliana.BioMart.plantsmart28")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
    BiocManager::install(packages[!installed_packages])
}
```

Meanwhile, adding the `build_vignettes=TRUE` when downloading
development version

``` r
BiocManager::install("shangguandong1996/FindIT2", build_vignettes=TRUE)
```

To view documentation of FindIT2, start R and enter:

``` r
browseVignettes("FindIT2")
```
