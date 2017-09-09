---
title: Subtracting input counts from ChIP samples
author: Aaron Lun
date: 19 August 2017
output: 
   html_document:
     fig_caption: false
---



# Background

This script performs a simulation to demonstrate the effects of subtracting control counts from ChIP counts in a simple case with equal baseline coverage.
We set up a simulation with two ChIP replicates in each of two groups and matching input samples.


```r
exp.type <- rep(c("ChIP", "Con", "ChIP", "Con"), each=2)
group.no <- rep(c("A", "B"), each=4)
groupings <- paste0(exp.type, group.no)
nlibs <- length(groupings)
```

We generate the mean vectors for DB and non-DB sites.
The background is the same between groups in all cases, and the only difference is that there is genuine binding in group A for DB sites.


```r
library(edgeR)
baseline <- 50
binding <- 50
mu.nodb <- rep(baseline, nlibs)
mu.nodb[exp.type=="ChIP"] <- baseline+binding
mu.db <- rep(baseline, nlibs)
mu.db[exp.type=="ChIP" & group.no=="A"] <- baseline+binding
```

Simulating counts, with an equal number of DB and non-DB sites:


```r
set.seed(1000)
P <- 1/0.1
is.null <- 1:10000
counts <- rbind(matrix(rnbinom(10000*nlibs, mu=mu.nodb, size=P), ncol=nlibs, byrow=TRUE),
                matrix(rnbinom(10000*nlibs, mu=mu.db, size=P), ncol=nlibs, byrow=TRUE))
```

# Running without subtraction

As a control, we do a vanilla analysis between the two groups directly.


```r
g <- factor(groupings)
design <- model.matrix(~0 + g)
colnames(design) <- levels(g)
```

Using _edgeR_:


```r
y.d <- DGEList(counts, lib.size=rep(1e6, nlibs))
y.d <- estimateDisp(y.d, design)
fit.d <- glmQLFit(y.d, design, robust=TRUE)
res.d <- glmQLFTest(fit.d, contrast=makeContrasts(ChIPA - ChIPB, levels=design))
summary(y.d$trended.dispersion)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.09733 0.09822 0.09960 0.09923 0.10008 0.10071
```

You can see that type I error is controlled for the true nulls.


```r
nullp.d <- res.d$table$PValue[is.null]
sum(nullp.d <= 0.01)/length(nullp.d) 
```

```
## [1] 0.01
```

```r
sum(nullp.d <= 0.05)/length(nullp.d)
```

```
## [1] 0.0461
```
    
At the same thresholds, there are more DB sites that get detected than non-DB sites.
This indicates that power is good.


```r
altp.d <- res.d$table$PValue[-is.null]
sum(altp.d <= 0.01)/length(altp.d)
```

```
## [1] 0.2761
```

```r
sum(altp.d <= 0.05)/length(altp.d)
```

```
## [1] 0.5178
```

# Running with subtraction

Now seeing what happens if we subtract counts before testing.


```r
subcounts <- counts
is.chip <- exp.type=="ChIP"
is.A <- group.no=="A"
subcounts[,is.chip & is.A] <- subcounts[,is.chip & is.A] - subcounts[,!is.chip & is.A]
subcounts[,is.chip & !is.A] <- subcounts[,is.chip & !is.A] - subcounts[,!is.chip & !is.A]
subcounts[subcounts < 0] <- 0
subcounts <- subcounts[,is.chip]
```

Setting up the new design matrix.


```r
g2 <- factor(groupings[is.chip])
design2 <- model.matrix(~0 + g2)
colnames(design2) <- levels(g2)
```

Running through _edgeR_:


```r
y.s <- DGEList(subcounts, lib.size=rep(1e6, length(g2)))
y.s <- estimateDisp(y.s, design2)
fit.s <- glmQLFit(y.s, design2, robust=TRUE)
res.s <- glmQLFTest(fit.s, contrast=makeContrasts(ChIPA - ChIPB, levels=design2))
summary(y.s$trended.dispersion)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.5194  0.9881  1.4736  1.5486  1.9523  3.6892
```

The results are now way too conservative, due to inflation of the dispersions.


```r
nullp.s <- res.s$table$PValue[is.null]
sum(nullp.s <= 0.01)/length(nullp.s) 
```

```
## [1] 0.0089
```

```r
sum(nullp.s <= 0.05)/length(nullp.s)
```

```
## [1] 0.0258
```

We see a concomitant reduction in power relative to the no-subtraction case.


```r
altp.s <- res.s$table$PValue[-is.null]
sum(altp.s <= 0.01)/length(altp.s)
```

```
## [1] 0.1579
```

```r
sum(altp.s <= 0.05)/length(altp.s)
```

```
## [1] 0.2689
```

Can the conservativeness upon subtraction be offset by simply increasing the threshold (notwithstanding the loss of interpretability of the error rates)?
No, based on AUC curves.


```r
thresholds <- 1:100/1000
tp.s <- findInterval(thresholds, sort(altp.s))/length(altp.s)
fp.s <- findInterval(thresholds, sort(nullp.s))/length(nullp.s)
tp.d <- findInterval(thresholds, sort(altp.d))/length(altp.d)
fp.d <- findInterval(thresholds, sort(nullp.d))/length(nullp.d)
plot(fp.s, tp.s, col="red", type="l", xlab="FPR", ylab="TPR", 
    xlim=c(0, 0.1), ylim=c(0, 1))
lines(fp.d, tp.d, col="blue")
```

![plot of chunk unnamed-chunk-14](figures-subtract/unnamed-chunk-14-1.png)

# Anticonservativeness with buffering

Buffering with lots of entries that are high-abundance and did not require much subtraction.


```r
others <- 1001:nrow(subcounts)
bufcounts <- subcounts
bufcounts[others,] <- matrix(rnbinom(length(others)*ncol(subcounts), mu=binding, size=P), length(others))
```

Running these through _edgeR_:


```r
y.b <- DGEList(bufcounts, lib.size=rep(1e6, length(g2)))
y.b <- estimateDisp(y.b, design2)
fit.b <- glmQLFit(y.b, design2, robust=TRUE)
res.b <- glmQLFTest(fit.b, contrast=makeContrasts(ChIPA - ChIPB, levels=design2))
summary(y.b$trended.dispersion)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1116  0.1131  0.1171  0.1221  0.1237  0.6110
```

We see loss of type I error control, because the buffering removes the protection from variance inflation.


```r
nullp.b <- res.b$table$PValue[-others]
sum(nullp.b <= 0.01)/length(nullp.b) 
```

```
## [1] 0.102
```

```r
sum(nullp.b <= 0.05)/length(nullp.b)
```

```
## [1] 0.187
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
##  [1] Rcpp_0.12.12    locfit_1.5-9.1  lattice_0.20-35 rprojroot_1.2  
##  [5] digest_0.6.12   grid_3.4.0      backports_1.1.0 magrittr_1.5   
##  [9] evaluate_0.10.1 highr_0.6       stringi_1.1.5   rmarkdown_1.6  
## [13] splines_3.4.0   statmod_1.4.30  tools_3.4.0     stringr_1.2.0  
## [17] yaml_2.1.14     compiler_3.4.0  htmltools_0.3.6 knitr_1.17
```
