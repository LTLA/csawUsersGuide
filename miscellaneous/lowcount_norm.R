# This script demonstrates that TMM normalization fails with low counts.
# Namely, does it accurately estaimte the normalization factors?

library(edgeR)

#########################################################
# Setting up a simulator function.

simulator <- function(ngenes, genes.spiked, mu.back, mult, disp) {
    mu.spike <- mu.back*mult
    x1 <- rnbinom(ngenes, mu=mu.back, size=1/disp)
    normed <- mu.back*ngenes/(mu.back*(ngenes-genes.spiked)+mu.spike*genes.spiked)
    x2 <- rnbinom(ngenes, mu=mu.back*normed, size=1/disp)
    spiked <- sample(ngenes, genes.spiked)
    x2[spiked] <- rnbinom(genes.spiked, mu=mu.spike*normed, size=1/disp)
    return(list(counts=cbind(x1, x2), factor=normed, spiked=spiked))
}

set.seed(10000)
ngenes <- 10000

#########################################################
# Simulating whether non-unity factors are accurately estimated.

lapply(1:3, FUN=function(i) {
    x <- simulator(ngenes, 200, 2, 5, 0.05) # low count
    calcNormFactors(x$counts)
})

lapply(1:3, FUN=function(i) {
    x <- simulator(ngenes, 200, 10, 5, 0.05) # low count
    calcNormFactors(x$counts)
})

lapply(1:3, FUN=function(i) {
    x <- simulator(ngenes, 200, 50, 5, 0.05) # low count
    calcNormFactors(x$counts)
})

c(1/sqrt(x$factor), sqrt(x$factor)) # Truth.

#########################################################
# No undersampling at all, just differences in library size.
# True values should be 1, 1.

lapply(1:3, FUN=function(i) {
    x <- matrix(rnbinom(ngenes*2, mu=c(1, 5), size=20), nrow=ngenes, ncol=2, byrow=TRUE)
    calcNormFactors(x)
})

lapply(1:3, FUN=function(i) {
    x <- matrix(rnbinom(ngenes*2, mu=c(10, 50), size=20), nrow=ngenes, ncol=2, byrow=TRUE)
    calcNormFactors(x)
})

#########################################################
# Wrapping up

sessionInfo()

