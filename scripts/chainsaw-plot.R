setwd("~/git/fluclades")

chainsaw <- read.csv("results/chainsaw-nsubtrees.csv")
chainsaw <- chainsaw[chainsaw$cutoff > 0.06,]

# manual runs of chainsaw.py
pdf("results/chainsaw.pdf", width=4.5, height=4.5)
par(mar=c(5,5,1,1))
#x <- c(0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25, 0.30)
#y <- c(  64,   43,   29,   25,   22,   21,   18,   16,   15,   13,   13,   12,    9)
plot(chainsaw$cutoff, chainsaw$nsubtrees, type='s', 
     xlab="Maximum internal branch length", 
     ylab="Number of subtrees", bty='n')
points(chainsaw$cutoff, chainsaw$nsubtrees, pch=19, cex=0.5)
abline(h=18, lty=2)
dev.off()

# examine distribution of serotype labels among subtrees
#labels <- read.csv("results/chainsaw-0.15.full.csv")
labels <- read.csv("results/chainsaw-0.18.labels.csv")
pat <- ".+_(H[0-9]+)N*[0-9]*_.+"
labels$serotype <- gsub(pat, "\\1", labels$tip.label)
labels$serotype[!grepl(pat, labels$tip.label)] <- NA

#require(xtable)
tab <- table(labels$subtree, labels$serotype)
#xtable(tab)

# generate a plot
x <- tab / apply(tab, 1, sum)
xval <- as.integer(gsub("H([0-9]+)", "\\1", colnames(x)))
xsum <- apply(x, 1, function(row) sum(xval*row))
io <- order(xsum)
jo <- order(xval)  # column index

#hc <- hclust(dist(tab))
pdf("results/chainsaw-table.pdf", width=4.5, height=4.5)
par(mar=c(5,5,1,2))
plot(NA, xlim=c(0, ncol(x)), ylim=c(0, nrow(x)), xaxt='n', yaxt='n',
     xlab="Serotype labels", ylab="Subtree", bty='n')
for (i in 1:nrow(x)) {
  for (j in 1:ncol(x)) {
    xx <- 1-x[io[i], jo[j]]
    if (xx < 1) xx <- min(0.9, xx) 
    rect(j-1, i-1, j, i, col=rgb(xx, xx, xx), border=NA)
    
    if (xx == 0.9) text(j-0.5, i-0.5, adj=0.5, label=tab[io[i], jo[j]], cex=0.5)
  }
}
for (i in 0:ncol(x)-1) {
  abline(v=i, col='grey80')
  abline(h=i, col='grey80')
}
axis(side=1, at=1:ncol(x)-0.5, label=paste("H", 1:18, sep=""), 
     cex.axis=0.8, las=2)
axis(side=2, at=1:nrow(x)-0.5, label=io-1, las=2, cex.axis=0.8)
#axis(side=4, at=1:nrow(x)-0.5, label=io-1, las=2, cex.axis=0.8)
axis(side=4, at=1:nrow(x)-0.5, label=apply(tab[io,], 1, sum),
     cex.axis=0.6, las=2, lwd=0, line=-1)
dev.off()

