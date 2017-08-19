---
title: Anticonservativeness in peak selection strategies
author: Aaron Lun
output: 
   html_document:
     fig.caption: false
---



# Background 

This simulation shows that ad hoc peak selection strategies not only result in conservatives (as shown in the NAR 2014 paper), but also loss of type I error control.
Consider an experimental design with two replicates in each of two conditions.


```r
group <- rep(c("A", "B"), each=2)
nlibs <- length(group)
design <- model.matrix(~group)
```

We'll simulate some counts where 10% of the sites are DB.


```r
library(edgeR)
set.seed(90000)
P <- 1/0.1
nsites <- 1e5
is.null <- seq_len(nsites)
counts <- rbind(matrix(rnbinom(nsites*nlibs, mu=50, size=P), ncol=nlibs, byrow=TRUE),
                matrix(rnbinom(nsites*0.1*nlibs, mu=c(100, 100, 0, 0), size=P), ncol=nlibs, byrow=TRUE))
```

We'll also set up a function to assess type I error control.


```r
plotAlpha <- function(pvals, ylab="Observed/specified", xlab="Specified", xlim=NULL, ...) {
    exp <- (seq_along(pvals) - 0.5)/length(pvals)
    n <- findInterval(exp, sort(pvals))
    obs <- n/length(pvals)
    if (is.null(xlim)) { # Stable at 20 observations.
        xlim <- c(exp[which(n >= 20)[1]], 1)
    }
    plot(exp, obs/exp, log="xy", xlim=xlim, ...)
}
```


# Applying the "at least 2" filter

Let's see what happens when we apply the "at least two" filter to retain the top proportion of sites.


```r
top.al2 <- apply(counts, 1, FUN=function(x) { sort(x)[nlibs-1] })
keep.al2 <- length(top.al2) - rank(top.al2, ties.method="random") +1 <= 20000
sum(which(keep.al2) %in% is.null)
```

```
## [1] 12729
```

Running these counts through _edgeR_, using standardized library sizes.
There is some inflation, but the presence of correct dispersion estimates for the DB sites keeps the dispersion low.


```r
y.al2 <- DGEList(counts[keep.al2,], lib.size=rep(1e6, nlibs))
y.al2 <- estimateDisp(y.al2, design)
y.al2$common.dispersion
```

```
## [1] 0.1138291
```

We observe loss of type I error control at low _p_-values.
This is because the dispersion inflation is minimized _and_ the "at least two" filter selects for spurious DB sites.


```r
fit.al2 <- glmQLFit(y.al2, design, robust=TRUE)
res.al2 <- glmQLFTest(fit.al2)
mean(res.al2$table$PValue[which(keep.al2) %in% is.null] <= 0.001)
```

```
## [1] 0.003456674
```

```r
plotAlpha(res.al2$table$PValue[which(keep.al2) %in% is.null])
```

![plot of chunk unnamed-chunk-7](figures-peak/unnamed-chunk-7-1.png)

# Applying a union filter.

Repeating the dose with a union filter.
Here we retain fewer sites, which ensures that the DB percentage in the retained set is higher (see below).


```r
top.u <- apply(counts, 1, FUN=function(x) { max(x) })
keep.u <- length(top.u) - rank(top.u, ties.method="random") +1 <= 5000
sum(which(keep.u) %in% is.null)
```

```
## [1] 410
```

Running these counts through _edgeR_.
The higher DB percentage keeps the dispersion inflation low.


```r
y.u <- DGEList(counts[keep.u,], lib.size=rep(1e6, nlibs))
y.u <- estimateDisp(y.u, design)
y.u$common.dispersion
```

```
## [1] 0.1630225
```

Testing again results in the loss of type I error control.
Normally, the union approach enriches outliers and inflates the dispersion.
However, enough DB sites ensures that the inflation is minimized, encouraging spurious rejection of the null.


```r
fit.u <- glmQLFit(y.u, design, robust=TRUE)
res.u <- glmQLFTest(fit.u)
mean(res.u$table$PValue[which(keep.u) %in% is.null] <= 0.05)
```

```
## [1] 0.1731707
```

```r
plotAlpha(res.u$table$PValue[which(keep.u) %in% is.null])
```

![plot of chunk unnamed-chunk-10](figures-peak/unnamed-chunk-10-1.png)

# Applying the mean filter

Now, to demonstrate the correct way of doing it, we use a filter on the mean count.


```r
top.m <- rowMeans(counts)
keep.m <- length(top.m) - rank(top.m, ties.method="random") +1 <= 10000
sum(which(keep.m) %in% is.null)
```

```
## [1] 8468
```

Running these counts through _edgeR_.
The higher DB percentage keeps the dispersion inflation low.


```r
y.m <- DGEList(counts[keep.m,], lib.size=rep(1e6, nlibs))
y.m <- estimateDisp(y.m, design)
y.m$common.dispersion
```

```
## [1] 0.0979567
```

Testing indicates that type I error control is mostly maintained.


```r
fit.m <- glmQLFit(y.m, design, robust=TRUE)
res.m <- glmQLFTest(fit.m)
mean(res.m$table$PValue[which(keep.m) %in% is.null] <= 0.01)
```

```
## [1] 0.01121871
```

```r
plotAlpha(res.m$table$PValue[which(keep.m) %in% is.null])
```

![plot of chunk unnamed-chunk-13](figures-peak/unnamed-chunk-13-1.png)

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
##  [1] Rcpp_0.12.12    locfit_1.5-9.1  lattice_0.20-35 rprojroot_1.2  
##  [5] digest_0.6.12   grid_3.4.0      backports_1.1.0 magrittr_1.5   
##  [9] evaluate_0.10.1 highr_0.6       stringi_1.1.5   rmarkdown_1.6  
## [13] splines_3.4.0   statmod_1.4.30  tools_3.4.0     stringr_1.2.0  
## [17] yaml_2.1.14     compiler_3.4.0  htmltools_0.3.6 knitr_1.17
```

