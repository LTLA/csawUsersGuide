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
One could consider erring on the side of conservativeness to be acceptable, especially if more DB sites can pass the filter.
However, this document demonstrates that the same strategies can also result in loss of type I error control.
The possibility of anticonservativeness means that any increased detection from _ad hoc_ strategies cannot be trusted.

# Setting up the experimental design

Consider an experimental design with two replicates in each of two conditions.


```r
group <- rep(c("A", "B"), each=2)
nlibs <- length(group)
design <- model.matrix(~group)
```

We set up a simulator for a count matrix where a certain proportion (10% by default) of the sites are DB.
Note that the grand mean for DB and non-DB sites are the same, otherwise the filtering problem would be trivial.


```r
simulateCounts <- function(nlibs, n.sites=1e5, prop.db=0.1, 
                           dispfun=function(x) { 0.1 }, 
                           n.mu=50, db.mu=rep(c(100, 0), each=2)) {
    P.n <- 1/dispfun(n.sites)
    db.sites <- n.sites*prop.db
    P.db <- 1/dispfun(db.sites)
    is.null <- seq_len(n.sites)
    counts <- rbind(matrix(rnbinom(n.sites*nlibs, mu=n.mu, size=P.n), ncol=nlibs, byrow=TRUE),
        matrix(rnbinom(db.sites*nlibs, mu=db.mu, size=P.db), ncol=nlibs, byrow=TRUE))
    return(list(counts=counts, null=is.null))
}
```

We set up a function to perform the differential analysis with _edgeR_ to obtain _p_-values.
We use equal library sizes here, assuming that normalization has already been performed to correct for composition biases.


```r
library(edgeR)
detectDiff <- function(counts, design, coef=ncol(design), lib.size=rep(1e6, ncol(counts))) {
    y <- DGEList(counts, lib.size=lib.size)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit, coef=coef)
    return(list(common.dispersion=y$common.dispersion, PValue=res$table$PValue))
}
```

We'll also set up a function to assess type I error control.


```r
plotAlpha <- function(pvals, ylab="Observed/specified", xlab="Specified", xlim=NULL, ...) {
    for (i in seq_along(pvals)) { 
        cur.p <- pvals[[i]]
        exp <- (seq_along(cur.p) - 0.5)/length(cur.p)
        n <- findInterval(exp, sort(cur.p))
        obs <- n/length(cur.p)
        if (is.null(xlim)) { # Stable at 20 observations.
            xlim <- c(exp[which(n >= 20)[1]], 1)
        }
        if (i==1L) {
            plot(exp, obs/exp, log="xy", xlim=xlim, type="l", ...)
        } else {
            lines(exp, obs/exp, ...)
        }
    }
}
```

# Applying the "at least 2" filter

Let's see what happens when we apply the "at least two" filter to retain the top proportion of sites.
First we set up a filtering function.


```r
set.seed(10001)
AL2Filter <- function(counts) { 
    top.al2 <- apply(counts, 1, FUN=function(x) { sort(x, decreasing=TRUE)[2] })
    rank(-top.al2, ties.method="random")
}
```

Now we run through repeated simulations and collect the results.


```r
retained <- dispersions <- numeric(10)
null.p <- vector("list", 10)
for (it in 1:10) {
    out <- simulateCounts(nrow(design))
    keep <- AL2Filter(out$counts) <= 20000
    kept.null <- which(keep) %in% out$null
    retained[it] <- sum(kept.null)
    res <- detectDiff(out$counts[keep,], design)
    dispersions[it] <- res$common.dispersion
    null.p[[it]] <- res$PValue[kept.null]
}
```

There is some inflation, but the presence of correct dispersion estimates for the DB sites keeps the common dispersion low.


```r
summary(retained)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   12730   12764   12819   12797   12824   12838
```

```r
summary(dispersions)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1140  0.1150  0.1159  0.1157  0.1165  0.1170
```

We observe loss of type I error control at low _p_-values.
This is because the dispersion inflation is minimized _and_ the "at least two" filter selects for spurious DB sites.


```r
summary(sapply(null.p, FUN=function(x) { mean(x <= 1e-3) }))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.002260 0.002920 0.003284 0.003298 0.003729 0.004057
```

```r
plotAlpha(null.p)
```

![plot of chunk unnamed-chunk-9](figures-peak/unnamed-chunk-9-1.png)

# Applying a union filter.

Repeating the dose with a union filter.
Here we retain fewer sites, which ensures that the DB percentage in the retained set is higher (see below).


```r
set.seed(20002)
UnionFilter <- function(counts) {
    top.u <- apply(counts, 1, FUN=max)
    rank(-top.u, ties.method="random") 
}
```

Now we run through repeated simulations and collect the results.
We use a more stringent filter to obtain a higher DB percentage, which keeps the dispersion inflation low.


```r
retained <- dispersions <- numeric(10)
null.p <- vector("list", 10)
for (it in 1:10) {
    out <- simulateCounts(nrow(design))
    keep <- UnionFilter(out$counts) <= 5000
    kept.null <- which(keep) %in% out$null
    retained[it] <- sum(kept.null)
    res <- detectDiff(out$counts[keep,], design)
    dispersions[it] <- res$common.dispersion
    null.p[[it]] <- res$PValue[kept.null]
}
```

Some inflation occurs, but all in all, the dispersions are kept reasonably low.


```r
summary(retained)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   388.0   440.8   448.5   443.1   454.8   473.0
```

```r
summary(dispersions)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1666  0.1686  0.1691  0.1696  0.1711  0.1724
```

Testing again results in the loss of type I error control.
Normally, the union approach enriches outliers and inflates the dispersion.
However, enough DB sites ensures that the inflation is minimized, encouraging spurious rejection of the null.


```r
summary(sapply(null.p, FUN=function(x) { mean(x <= 1e-2) }))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.01418 0.01770 0.02385 0.02235 0.02559 0.02908
```

```r
plotAlpha(null.p)
```

![plot of chunk unnamed-chunk-13](figures-peak/unnamed-chunk-13-1.png)

# Applying the mean filter

Now, to demonstrate the correct way of doing it, we use a filter on the mean count.


```r
set.seed(30003)
MeanFilter <- function(counts) { 
    top.m <- rowMeans(counts)
    rank(-top.m, ties.method="random")
}
```

Running these counts through _edgeR_.


```r
retained <- dispersions <- numeric(10)
null.p <- vector("list", 10)
for (it in 1:10) {
    out <- simulateCounts(nrow(design))
    keep <- MeanFilter(out$counts) <= 10000
    kept.null <- which(keep) %in% out$null
    retained[it] <- sum(kept.null)
    res <- detectDiff(out$counts[keep,], design)
    dispersions[it] <- res$common.dispersion
    null.p[[it]] <- res$PValue[kept.null]
}
```

Dispersions are equal to their expected value.


```r
summary(retained)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    8479    8494    8509    8515    8537    8554
```

```r
summary(dispersions)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.09713 0.09754 0.09917 0.09898 0.10021 0.10100
```

Testing indicates that type I error control is mostly maintained.


```r
summary(sapply(null.p, FUN=function(x) { mean(x <= 1e-2) }))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.007963 0.008709 0.009973 0.009548 0.010281 0.010653
```

```r
plotAlpha(null.p, ylim=c(0.5, 2))
```

![plot of chunk unnamed-chunk-17](figures-peak/unnamed-chunk-17-1.png)

# Saving _ad hoc_ filters with FDR control

One could argue that anticonservativeness with _ad hoc_ filters is not a problem in practice, because the enrichment for DB sites ensures that FDR control is still preserved.
If you enrich for enough DB sites, even complete loss of type I error control among the true nulls will not breach the FDR threshold.
However, this assumes that you have enough power to detect all of the DB sites.
If you don't have enough power, FDR control is lost:


```r
set.seed(40004)
out <- simulateCounts(nrow(design), db.mu=c(70, 70, 30, 30))
keep <- AL2Filter(out$counts) <= 20000
res <- detectDiff(out$counts[keep,], design)
sig <- p.adjust(res$PValue, method="BH") <= 0.05
sum(sig & which(keep) %in% out$null)/sum(sig)
```

```
## [1] 0.08518296
```

There are also other ways of mitigating the variance inflation that don't involve introducing DB sites between two conditions.
For example, if an orthogonal batch effect increases binding in some sites, those sites will get preferentially enriched by an "at least 2" filter.
This will suppress the dispersion inflation in the other sites and encourage loss of type I error control.
FDR control is irrelevant here as there are no DB sites to the contrast between the first two groups.


```r
design.b <- model.matrix(~c(0,1,0,1)+c(0,0,1,1))
null.p <- vector("list", 10)
for (it in 1:10) {
    out <- simulateCounts(nrow(design.b), db.mu=rep(c(20, 80), 2)) # not really DB.
    keep <- AL2Filter(out$counts) <= 10000
    res <- detectDiff(out$counts[keep,], design.b)
    null.p[[it]] <- res$PValue # look at all p-values, as all nulls are true.
}
plotAlpha(null.p)
```

![plot of chunk unnamed-chunk-19](figures-peak/unnamed-chunk-19-1.png)

Alternatively, you could imagine a situation with three groups, and DB sites in the third group can suppress dispersion inflation.
This is done without contributing differences to the contrast between first two groups.


```r
group.3 <- factor(rep(LETTERS[1:3], each=2))
design.3 <- model.matrix(~group.3)
null.p <- vector("list", 10)
for (it in 1:10) {
    out <- simulateCounts(nrow(design.3), db.mu=rep(c(25, 100), c(4, 2))) # not really DB.
    keep <- AL2Filter(out$counts) <= 10000
    kept.null <- which(keep) %in% out$null
    res <- detectDiff(out$counts[keep,], design.3, coef=2)
    null.p[[it]] <- res$PValue # look at all p-values, as all nulls are true.
}
plotAlpha(null.p)
```

![plot of chunk unnamed-chunk-20](figures-peak/unnamed-chunk-20-1.png)

# Using the mean filter with variable dispersions

The sample mean is approxiamtely independent of the dispersion estimate and p-value for Poisson and NB-distributed counts.
However, this only applies to a single distribution.
Consider a mixture of NB distributions with different dispersions following an inverse-chi-squared distribution but the same mean.


```r
set.seed(50005)
PFUN <- function(n) { 1/rchisq(n, df=10) }
mean(PFUN(10000))
```

```
## [1] 0.1229193
```

Running on a randomly selected subset of sites (i.e., completely unbiased filtering) doesn't cause any issues other than some anticonservativeness at low _p_-values.
This is relatively modest and acceptable given that the inverse-chi-squared distribution for the NB dispersions only mimics the true QL model.
(In contrast, the _ad hoc_ filters were tested on purely NB counts, so they should not have any problems.)


```r
dispersions <- numeric(10)
null.p <- vector("list", 10)
for (it in 1:10) {
    out <- simulateCounts(nrow(design), dispfun=PFUN)
    keep <- sample(nrow(out$counts), 1000)
    kept.null <- keep %in% out$null
    res <- detectDiff(out$counts[keep,], design)
    dispersions[it] <- res$common.dispersion
    null.p[[it]] <- res$PValue[kept.null]
}
summary(dispersions)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1163  0.1183  0.1207  0.1215  0.1247  0.1289
```

```r
plotAlpha(null.p, ylim=c(0.5, 2))
```

![plot of chunk unnamed-chunk-22](figures-peak/unnamed-chunk-22-1.png)

However, applying a stringent mean filter will select for higher dispersions.
This is because high-dispersion features are more likely to achieve large sample means. 
The result is to encourage inflation of the dispersion estimate.


```r
dispersions <- numeric(10)
null.p <- vector("list", 10)
for (it in 1:10) {
    out <- simulateCounts(nrow(design), dispfun=PFUN)
    keep <- MeanFilter(out$counts) <= 1000
    kept.null <- which(keep) %in% out$null
    res <- detectDiff(out$counts[keep,], design)
    dispersions[it] <- res$common.dispersion
    null.p[[it]] <- res$PValue[kept.null]
}
summary(dispersions)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1798  0.1906  0.1922  0.1917  0.1945  0.1967
```

Funnily enough, this doesn't actually hurt the type I error rate.
This is possibly because the same issue affects both DB and non-DB sites, i.e., DB sites with high dispersions will be similarly preferred.
This ensures that there is no enrichment for low-dispersion DB sites to suppress the variance inflation and lead to anticonservativeness.


```r
plotAlpha(null.p, ylim=c(0.5, 2))
```

![plot of chunk unnamed-chunk-24](figures-peak/unnamed-chunk-24-1.png)

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

