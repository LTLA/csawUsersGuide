---
title: TMM normalization with low counts
author: Aaron Lun
date: 19 August 2017
output: 
   html_document:
     fig.caption: false
---



# Background 

This script demonstrates that TMM normalization fails with low counts.
Namely, does it accurately estimate the normalization factors?
We set up a simulation where some genes are upregulated by `mult` in the second sample.


```r
simulator <- function(ngenes, genes.spiked, mu.back, mult, disp) {
    mu.spike <- mu.back*mult
    x1 <- rnbinom(ngenes, mu=mu.back, size=1/disp)
    normed <- mu.back*ngenes/(mu.back*(ngenes-genes.spiked)+mu.spike*genes.spiked)
    x2 <- rnbinom(ngenes, mu=mu.back*normed, size=1/disp)
    spiked <- sample(ngenes, genes.spiked)
    x2[spiked] <- rnbinom(genes.spiked, mu=mu.spike*normed, size=1/disp)
    return(list(counts=cbind(x1, x2), factor=normed, spiked=spiked))
}
```

# Simulating with various count sizes

Simulating over three times for various count sizes:


```r
library(edgeR)
set.seed(1000)
ngenes <- 10000
lapply(1:3, FUN=function(i) {
    x <- simulator(ngenes, 200, 2, 5, 0.05) # low count
    calcNormFactors(x$counts)
})
```

```
## [[1]]
## [1] 1.0013521 0.9986498
## 
## [[2]]
## [1] 1.0069616 0.9930865
## 
## [[3]]
## [1] 1.0052563 0.9947712
```

```r
lapply(1:3, FUN=function(i) {
    x <- simulator(ngenes, 200, 10, 5, 0.05) # middle count
    calcNormFactors(x$counts)
})
```

```
## [[1]]
## [1] 1.0291659 0.9716607
## 
## [[2]]
## [1] 1.029547 0.971301
## 
## [[3]]
## [1] 1.0341057 0.9670192
```

```r
lapply(1:3, FUN=function(i) {
    x <- simulator(ngenes, 200, 50, 5, 0.05) # high count
    calcNormFactors(x$counts)
})
```

```
## [[1]]
## [1] 1.0356454 0.9655815
## 
## [[2]]
## [1] 1.0345765 0.9665791
## 
## [[3]]
## [1] 1.0356693 0.9655591
```

We then compare these values to the truth.
The lower counts do not perform well, due to the low precision for trimming when M-values are discrete.
The shift of the median with unbalanced DE is also more pronounced when the non-DE M-values are more variable.


```r
x <- simulator(ngenes, 200, 50, 5, 0.05) 
c(1/sqrt(x$factor), sqrt(x$factor)) # Truth.
```

```
## [1] 1.0392305 0.9622504
```

# Failure even without understampling

Consider these simulations where there is no undersampling at all, just differences in library size.
True normalization factors should be 1, but this is not the case, corresponding to loss of precision in trimming.


```r
lapply(1:3, FUN=function(i) {
    x <- matrix(rnbinom(ngenes*2, mu=c(1, 5), size=20), nrow=ngenes, ncol=2, byrow=TRUE)
    calcNormFactors(x)
})
```

```
## [[1]]
## [1] 1.229705 0.813203
## 
## [[2]]
## [1] 1.2368873 0.8084811
## 
## [[3]]
## [1] 1.2328770 0.8111109
```

```r
lapply(1:3, FUN=function(i) {
    x <- matrix(rnbinom(ngenes*2, mu=c(10, 50), size=20), nrow=ngenes, ncol=2, byrow=TRUE)
    calcNormFactors(x)
})
```

```
## [[1]]
## [1] 0.9935526 1.0064892
## 
## [[2]]
## [1] 0.9947951 1.0052322
## 
## [[3]]
## [1] 0.9936823 1.0063579
```

# Wrapping up


```r
sessionInfo()
```

```
## R version 3.4.0 Patched (2017-04-24 r72627)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.5 LTS
## 
## Matrix products: default
## BLAS: /home/cri.camres.org/lun01/Software/R/R-3-4-branch_devel/lib/libRblas.so
## LAPACK: /home/cri.camres.org/lun01/Software/R/R-3-4-branch_devel/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] edgeR_3.19.3 limma_3.33.7
## 
## loaded via a namespace (and not attached):
##  [1] compiler_3.4.0  magrittr_1.5    tools_3.4.0     stringi_1.1.5  
##  [5] grid_3.4.0      locfit_1.5-9.1  knitr_1.17      stringr_1.2.0  
##  [9] lattice_0.20-35 evaluate_0.10.1
```

sessionInfo()

