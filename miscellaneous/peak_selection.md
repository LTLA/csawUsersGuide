---
title: Anticonservativeness in peak selection strategies
author: Aaron Lun
date: 19 August 2017
output: 
   html_document:
     fig_caption: false
---



# Background 

The 2014 NAR paper (https://doi.org/10.1093/nar/gku351) showed that many _ad hoc_ peak selection strategies result in conservativeness.
One could consider erring on the side of conservativeness to be acceptable, especially if more DB sites pass the filter.
However, this simulation demonstrates that the same strategies can also result in loss of type I error control.
The possibility of anticonservativeness means that any increased detection from _ad hoc_ strategies cannot be trusted.
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
n.sites <- 1e5
db.sites <- n.sites*0.1
is.null <- seq_len(n.sites)
counts <- rbind(matrix(rnbinom(n.sites*nlibs, mu=50, size=P), ncol=nlibs, byrow=TRUE),
                matrix(rnbinom(db.sites*nlibs, mu=c(100, 100, 0, 0), size=P), ncol=nlibs, byrow=TRUE))
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
null.al2 <- which(keep.al2) %in% is.null
summary(null.al2)
```

```
##    Mode   FALSE    TRUE 
## logical    7271   12729
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
mean(res.al2$table$PValue[null.al2] <= 0.001)
```

```
## [1] 0.003456674
```

```r
plotAlpha(res.al2$table$PValue[null.al2])
```

![plot of chunk unnamed-chunk-7](figures-peak/unnamed-chunk-7-1.png)

One could argue that this is not a problem in practice, because the enrichment for DB sites ensures that FDR control is still preserved.
In the most extreme case, if you enrich for enough DB sites, even complete loss of type I error control among the true nulls will not breach the FDR threshold.
However, this assumes that you have enough power to detect all of the DB sites.
This is true in this particular simulation, where the DB is very strong - try setting `mu=c(70, 70, 30, 30)` for comparison.


```r
sig <- p.adjust(res.al2$table$PValue, method="BH") <= 0.05
table(sig, null.al2)
```

```
##        null.al2
## sig     FALSE  TRUE
##   FALSE     0 12353
##   TRUE   7271   376
```

```r
sum(sig & null.al2)/sum(sig)
```

```
## [1] 0.04916961
```

There are also other ways of mitigating the variance inflation that don't involve introducing DB sites between the two conditions.
For example, in an experimental design with multiple groups, you could add DB sites in the third group.
These would negate variance inflation but not contribute DB sites to the contrast between the first two groups.

# Applying a union filter.

Repeating the dose with a union filter.
Here we retain fewer sites, which ensures that the DB percentage in the retained set is higher (see below).


```r
top.u <- apply(counts, 1, FUN=function(x) { max(x) })
keep.u <- length(top.u) - rank(top.u, ties.method="random") +1 <= 5000
null.u <- which(keep.u) %in% is.null
summary(null.u)
```

```
##    Mode   FALSE    TRUE 
## logical    4590     410
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
mean(res.u$table$PValue[null.u] <= 0.05)
```

```
## [1] 0.1731707
```

```r
plotAlpha(res.u$table$PValue[null.u])
```

![plot of chunk unnamed-chunk-11](figures-peak/unnamed-chunk-11-1.png)

# Applying the mean filter

Now, to demonstrate the correct way of doing it, we use a filter on the mean count.


```r
top.m <- rowMeans(counts)
keep.m <- length(top.m) - rank(top.m, ties.method="random") +1 <= 10000
null.m <- which(keep.m) %in% is.null
summary(null.m)
```

```
##    Mode   FALSE    TRUE 
## logical    1532    8468
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
mean(res.m$table$PValue[null.m] <= 0.01)
```

```
## [1] 0.01121871
```

```r
plotAlpha(res.m$table$PValue[null.m])
```

![plot of chunk unnamed-chunk-14](figures-peak/unnamed-chunk-14-1.png)

# Using the mean filter with variable dispersions

The sample mean is approxiamtely independent of the dispersion estimate and p-value for Poisson and NB-distributed counts.
However, this only applies to a single distribution.
Consider a mixture of NB distributions with different dispersions but the same mean.


```r
set.seed(100000)
P.n <- rchisq(n.sites, df=10)
P.db <- rchisq(db.sites, df=10)
counts <- rbind(matrix(rnbinom(n.sites*nlibs, mu=50, size=P.n), ncol=nlibs, byrow=TRUE),
                matrix(rnbinom(db.sites*nlibs, mu=c(100, 100, 0, 0), size=P.db), ncol=nlibs, byrow=TRUE))
```

Running on all the sites doesn't cause any issues other than some anticonservativeness at low _p_-values.
This is relatively modest and acceptable given that the inverse-chi-squared distribution for the NB dispersions only mimics the true QL model.
(In contrast, the _ad hoc_ filters were tested on purely NB counts, so they should not have any problems.)


```r
y.all <- DGEList(counts, lib.size=rep(1e6, nlibs))
y.all <- estimateDisp(y.all, design)
y.all$common.dispersion
```

```
## [1] 0.1234853
```

```r
fit.all <- glmQLFit(y.all, design, robust=TRUE)
res.all <- glmQLFTest(fit.all)
mean(res.all$table$PValue[is.null] <= 0.01)
```

```
## [1] 0.01121
```

```r
plotAlpha(res.all$table$PValue[is.null])
```

![plot of chunk unnamed-chunk-16](figures-peak/unnamed-chunk-16-1.png)

However, applying a stringent mean filter will select for higher dispersions.
This is because high-dispersion features are more likely to achieve large sample means. 
The result is to encourage inflation of the dispersion estimate.


```r
top.m2 <- rowMeans(counts)
keep.m2 <- length(top.m2) - rank(top.m2, ties.method="random") +1 <= 1000
summary(which(keep.m2) %in% is.null)
```

```
##    Mode   FALSE    TRUE 
## logical     289     711
```

```r
y.m2 <- DGEList(counts[keep.m2,], lib.size=rep(1e6, nlibs))
y.m2 <- estimateDisp(y.m2, design)
y.m2$common.dispersion
```

```
## [1] 0.1907132
```

Funnily enough, this doesn't actually hurt the type I error rate.
This is possibly because the same issue affects both DB and non-DB sites, i.e., DB sites with high dispersions will be similarly preferred.
This ensures that there is no enrichment for low-dispersion DB sites to suppress the variance inflation and lead to anticonservativeness.


```r
fit.m2 <- glmQLFit(y.m2, design, robust=TRUE)
res.m2 <- glmQLFTest(fit.m2)
mean(res.m2$table$PValue[which(keep.m2) %in% is.null] <= 0.01)
```

```
## [1] 0.0140647
```

```r
plotAlpha(res.m2$table$PValue[which(keep.m2) %in% is.null])
```

![plot of chunk unnamed-chunk-18](figures-peak/unnamed-chunk-18-1.png)

In practice, this is not much of an issue.
For ChIP-seq data, the prior degrees of freedom is usually quite high such that there is not much variability in the dispersions.
For RNA-seq data, the filter boundary is not dense so the specifics of filtering doesn't matter.
The same effect is present in the other filters anyway (in addition to their poor performance on the constant dispersion case),
    as high dispersions make it more likely to get one or two peaks above the threshold.

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
##  [1] compiler_3.4.0  magrittr_1.5    tools_3.4.0     splines_3.4.0  
##  [5] stringi_1.1.5   highr_0.6       grid_3.4.0      locfit_1.5-9.1 
##  [9] knitr_1.17      stringr_1.2.0   statmod_1.4.30  lattice_0.20-35
## [13] evaluate_0.10.1
```

