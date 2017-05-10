# This script performs a simulation to demonstrate the effects of subtracting
# control counts from ChIP counts in a simple case with equal baseline coverage.

library(edgeR)
set.seed(1000)

baseline <- 50
binding <- 50
exp.type <- rep(c("ChIP", "Con", "ChIP", "Con"), each=2)
group.no <- rep(c("A", "B"), each=4)
groupings <- paste0(exp.type, group.no)
nlibs <- length(groupings)

mu.nodb <- rep(baseline, nlibs)
mu.nodb[exp.type=="ChIP"] <- baseline+binding
mu.db <- rep(baseline, nlibs)
mu.db[exp.type=="ChIP" & group.no=="A"] <- baseline+binding

P <- 1/0.1
counts <- rbind(matrix(rnbinom(10000*nlibs, mu=mu.nodb, size=P), ncol=nlibs, byrow=TRUE),
                matrix(rnbinom(10000*nlibs, mu=mu.db, size=P), ncol=nlibs, byrow=TRUE))
 
#######################################                   
# Running without subtraction

g <- factor(groupings)
design <- model.matrix(~0 + g)
colnames(design) <- levels(g)

y.d <- DGEList(counts, lib.size=rep(1e6, nlibs))
y.d <- estimateDisp(y.d, design)
fit.d <- glmQLFit(y.d, design, robust=TRUE)
res.d <- glmQLFTest(fit.d, contrast=makeContrasts(ChIPA - ChIPB, levels=design))

# Type I error controlled.
is.null <- 1:10000
nullp.d <- res.d$table$PValue[is.null]
sum(nullp.d <= 0.01)/length(nullp.d) 
sum(nullp.d <= 0.05)/length(nullp.d)
summary(fit.d$dispersion)
summary(fit.d$df.prior)
    
# Alternative p-values mostly below null p-values at the same threshold.
altp.d <- res.d$table$PValue[-is.null]
sum(altp.d <= 0.01)/length(altp.d)
sum(altp.d <= 0.05)/length(altp.d)
mean(findInterval(altp.d[altp.d <= 0.01], sort(nullp.d[nullp.d <= 0.01]))/sum(nullp.d <= 0.01))
mean(findInterval(altp.d[altp.d <= 0.05], sort(nullp.d[nullp.d <= 0.05]))/sum(nullp.d <= 0.05))

#######################################                   
# Running with subtraction

subcounts <- counts
is.chip <- exp.type=="ChIP"
is.A <- group.no=="A"
subcounts[,is.chip & is.A] <- subcounts[,is.chip & is.A] - subcounts[,!is.chip & is.A]
subcounts[,is.chip & !is.A] <- subcounts[,is.chip & !is.A] - subcounts[,!is.chip & !is.A]
subcounts[subcounts < 0] <- 0
subcounts <- subcounts[,is.chip]

g2 <- factor(groupings[is.chip])
design2 <- model.matrix(~0 + g2)
colnames(design2) <- levels(g2)

y.s <- DGEList(subcounts, lib.size=rep(1e6, length(g2)))
y.s <- estimateDisp(y.s, design2)
fit.s <- glmQLFit(y.s, design2, robust=TRUE)
res.s <- glmQLFTest(fit.s, contrast=makeContrasts(ChIPA - ChIPB, levels=design2))

# Way too conservative, due to inflation of variances.
nullp.s <- res.s$table$PValue[is.null]
sum(nullp.s <= 0.01)/length(nullp.s) 
sum(nullp.s <= 0.05)/length(nullp.s)
summary(fit.s$dispersion)
summary(fit.s$df.prior)
    
# More alternative p-values above the null p-values.
# (i.e., more false positives, if threshold is increased to offset conservativeness).
altp.s <- res.s$table$PValue[-is.null]
sum(altp.s <= 0.01)/length(altp.s)
sum(altp.s <= 0.05)/length(altp.s)
mean(findInterval(altp.s[altp.s <= 0.01], sort(nullp.s[nullp.s <= 0.01]))/sum(nullp.s <= 0.01))
mean(findInterval(altp.s[altp.s <= 0.05], sort(nullp.s[nullp.s <= 0.05]))/sum(nullp.s <= 0.05))

#######################################
# Buffering with lots of entries that are high-abundance and did not require much subtraction.

others <- 1001:nrow(subcounts)
bufcounts <- subcounts
bufcounts[others,] <- matrix(rnbinom(length(others)*ncol(subcounts), mu=binding, size=P), length(others))

y.b <- DGEList(bufcounts, lib.size=rep(1e6, length(g2)))
y.b <- estimateDisp(y.b, design2)
fit.b <- glmQLFit(y.b, design2, robust=TRUE)
res.b <- glmQLFTest(fit.b, contrast=makeContrasts(ChIPA - ChIPB, levels=design2))

nullp.b <- res.b$table$PValue[-others]
sum(nullp.b <= 0.01)/length(nullp.b) 
sum(nullp.b <= 0.05)/length(nullp.b)
summary(fit.b$dispersion)
summary(fit.b$df.prior)

