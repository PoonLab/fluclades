setwd("~/git/fluclades")

rescale <- function(x, from, to) {
  diff(range(to)) * (x-min(from)) / diff(range(from)) + min(to)
}

chainsaw <- read.csv("results/chainsaw-nsubtrees.csv")
chainsaw <- chainsaw[chainsaw$cutoff > 0.06,]
#chainsaw <- read.csv("results/chainsaw-nsubtrees-na.csv")

# manual runs of chainsaw.py
pdf("results/chainsaw-ha.pdf", width=5, height=4)
par(mar=c(5,5,1,5))
#x <- c(0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25, 0.30)
#y <- c(  64,   43,   29,   25,   22,   21,   18,   16,   15,   13,   13,   12,    9)
plot(chainsaw$cutoff, chainsaw$nsubtrees, type='n', 
     xlab="Maximum internal branch length", ylab=NA,
     bty='n', yaxt='n', cex.axis=0.8)
axis(side=2, col='royalblue3', col.axis='royalblue3', las=2, cex.axis=0.8)
mtext(side=2, text="Number of subtrees", line=3, col='royalblue3')

# rescale 
y <- rescale(chainsaw$normalized, chainsaw$normalized, chainsaw$nsubtrees)
polygon(c(min(chainsaw$cutoff), chainsaw$cutoff, max(chainsaw$cutoff)), 
        c(0, y, 0), col=rgb(0.9,0,0,0.2), border=NA)
p <- pretty(chainsaw$normalized)
axis(side=4, at=rescale(p, chainsaw$normalized, chainsaw$nsubtrees), 
     label=p, cex.axis=0.7, las=2,
     col='firebrick3', col.axis='firebrick3')
u <- par('usr')
text(x=u[2]+0.06, y=mean(u[3:4]), 
     label="Normalized mutual information", col='firebrick3', 
     srt=-90, xpd=NA)

lines(chainsaw$cutoff, chainsaw$nsubtrees, type='s', col='royalblue')
points(chainsaw$cutoff, chainsaw$nsubtrees, pch=19, cex=0.5, col='royalblue')

abline(h=18, lty=2)
dev.off()


######  NA  ######
chainsaw <- read.csv("results/chainsaw-nsubtrees-na.csv")

pdf("results/chainsaw-na.pdf", width=5, height=4)
par(mar=c(5,5,1,5))
#x <- c(0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25, 0.30)
#y <- c(  64,   43,   29,   25,   22,   21,   18,   16,   15,   13,   13,   12,    9)
plot(chainsaw$cutoff, chainsaw$nsubtrees, type='n', 
     xlab="Maximum internal branch length", ylab=NA,
     bty='n', yaxt='n', cex.axis=0.8)
axis(side=2, col='royalblue3', col.axis='royalblue3', las=2, cex.axis=0.8)
mtext(side=2, text="Number of subtrees", line=3, col='royalblue3')

# rescale 
foo <- chainsaw$normalized[chainsaw$normalized > 0.1]
y <- rescale(chainsaw$normalized, foo, chainsaw$nsubtrees)
polygon(c(min(chainsaw$cutoff), chainsaw$cutoff, max(chainsaw$cutoff)), 
        c(0, y, 0), col=rgb(0.9,0,0,0.2), border=NA)
p <- pretty(chainsaw$normalized)
axis(side=4, at=rescale(p, chainsaw$normalized, chainsaw$nsubtrees), 
     label=p, cex.axis=0.7, las=2,
     col='firebrick3', col.axis='firebrick3')
u <- par('usr')
text(x=u[2]+0.09, y=mean(u[3:4]), 
     label="Normalized mutual information", col='firebrick3', 
     srt=-90, xpd=NA)

lines(chainsaw$cutoff, chainsaw$nsubtrees, type='s', col='royalblue')
points(chainsaw$cutoff, chainsaw$nsubtrees, pch=19, cex=0.5, col='royalblue')

abline(h=11, lty=2)
dev.off()



# examine distribution of serotype labels among subtrees
#labels <- read.csv("results/chainsaw-0.15.full.csv")
labels <- read.csv("results/chainsaw-HA-0.2.csv")
labels <- read.csv("results/chainsaw-NA-0.4.labels.csv")
pat <- ".+_(H[0-9]+)N*[0-9]*_.+" # HA
pat <- ".+_H[0-9]+(N[0-9]+)_.+"  # NA
labels$serotype <- gsub(pat, "\\1", labels$tip.label)
labels$serotype[!grepl(pat, labels$tip.label)] <- NA

#require(xtable)
tab <- table(labels$subtree, labels$serotype)
#xtable(tab)

# generate a plot
x <- tab / apply(tab, 1, sum)
#xval <- as.integer(gsub("H([0-9]+)", "\\1", colnames(x)))
xval <- as.integer(gsub("N([0-9]+)", "\\1", colnames(x)))
xsum <- apply(x, 1, function(row) sum(xval*row))
io <- order(xsum)
jo <- order(xval)  # column index

#hc <- hclust(dist(tab))
pdf("results/chainsaw-NA-table.pdf", width=4.5, height=4.5)
par(mar=c(5,5,1,2))
plot(NA, xlim=c(0, ncol(x)), ylim=c(0, nrow(x)), xaxt='n', yaxt='n',
     xlab="Serotype labels", ylab="Subtree", bty='n')
for (i in 1:nrow(x)) {
  for (j in 1:ncol(x)) {
    xx <- 1-x[io[i], jo[j]]
    if (xx < 1) xx <- min(0.9, xx) 
    rect(j-1, i-1, j, i, col=rgb(xx, xx, xx), border=NA)
    count <- tab[io[i], jo[j]]
    if (count > 0) text(j-0.5, i-0.5, adj=0.5, label=count, cex=0.5)
  }
}
for (i in 0:ncol(x)-1) {
  abline(v=i, col='grey80')
  abline(h=i, col='grey80')
}
#axis(side=1, at=1:ncol(x)-0.5, label=paste("H", 1:18, sep=""), 
#     cex.axis=0.8, las=2)
axis(side=1, at=1:ncol(x)-0.5, label=paste("N", 1:11, sep=""), 
     cex.axis=0.8, las=2)
axis(side=2, at=1:nrow(x)-0.5, label=io-1, las=2, cex.axis=0.8)
#axis(side=4, at=1:nrow(x)-0.5, label=io-1, las=2, cex.axis=0.8)
axis(side=4, at=1:nrow(x)-0.5, label=apply(tab[io,], 1, sum),
     cex.axis=0.6, las=2, lwd=0, line=-1)
dev.off()


# calculate mutual information
mutinfo <- function(labels) {
  n <- sum(!is.na(labels$serotype))
  pu <- table(labels$serotype) / n  # serotypes
  pv <- table(labels$subtree[!is.na(labels$serotype)]) / n  # subtrees
  mi <- 0
  for (sero in names(pu)) {
    for (st in names(pv)) {
      puv <- nrow(labels[labels$subtree==st & 
                           !is.na(labels$serotype) & 
                           labels$serotype==sero,]) / n
      if (puv == 0) next
      mi <- mi + puv * log(puv / (pu[sero] * pv[st]))
    }
  }
  return(mi)
}

mutinfo(labels)


# examine NA subtree 6
st <- 6
foo <- gsub(".+_(H[0-9]+)N.+", "\\1", labels$tip.label[labels$subtree==st])
foo[!grepl(".+_(H[0-9]+)N.+", labels$tip.label[labels$subtree==st])] <- NA
table(foo, labels$serotype[labels$subtree==st])

ha.sero <- gsub(".+_(H[0-9]+)N.+", "\\1", labels$tip.label)
ha.sero[!grepl(".+_(H[0-9]+)N.+", labels$tip.label)] <- NA
labels$ha.sero <- ha.sero

labels$accn <- sapply(labels$tip.label, function(x) strsplit(x, "_")[[1]][1])

meta <- read.csv("data/gb-metadata-na.csv")
idx <- match(labels$accn, meta$accn)
labels$host <- meta$host[idx]
labels$is.human <- grepl("Homo", labels$host)

tab <- table(labels$ha.sero[!labels$is.human], labels$serotype[!labels$is.human])
xtable(tab[c(1, 11:18, 2:10), c(1, 4:11, 2, 3)])

tab <- table(labels$ha.sero[labels$is.human], labels$serotype[labels$is.human])
xtable(tab[c(1, 3:7, 2), ])


phy <- read.tree("data/na.pruned.nwk")
phy <- ladderize(phy)
L <- tree.layout(phy)
idx <- match(L$nodes$label, labels$tip.label)
L$nodes$subtree <- labels$subtree[idx]
L$nodes$label <- gsub("(.+H[0-9]+N[0-9]+)_.+", "\\1", L$nodes$label)

tips <- L$nodes$label

set.seed(6)
pal <- sample(gg.rainbow(11), 11, replace=FALSE)

pdf("results/na-tree.pdf", width=7, height=9)
plot(L, mar=c(0,0,0,11), lwd=1, cex=0.4, col='grey')
for (i in 0:10) {
  draw.clade(L, tips=tips[which(L$nodes$subtree==i)], 
             col=pal[i+1], lwd=1.5)  
}
add.scalebar(L, len=0.2, x0=0.5, y0=10, cex=0.7)
dev.off()

############# PB2 ###############

others <- read.csv("results/chainsaw-nsubtrees-others.csv")
others <- others[order(others$gene, others$cutoff), ]

marks <- list(
  PB2=c(5, 0.046),
  PB1=c(4, 0.04),
  PA=c(5, 0.043),
  NP=c(4, 0.055),
  M1M2=c(6, 0.045),
  NS1NS2=c(5, 0.15)
)

pdf("~/papers/fluclades/chainsaw-others.pdf", width=11, height=7.3)
pal <- gg.rainbow(n=6, l=50)
par(mfrow=c(2,3), mar=c(5,5,1,1), cex=1)
i <- 1
for (g in c('PB2', 'PB1', 'PA', 'NP', 'M1M2', 'NS1NS2')) {
  rows <- others[others$gene==g, ]
  plot(NA, xlim=c(0.01, 0.2), ylim=c(2, 25), bty='n', 
       xlab="Maximum internal branch length",
       ylab="Number of subtrees", )
  lines(rows$cutoff, rows$nsubtrees, type='s', col=pal[i])
  title(main=g, font.main=1, adj=0.9, line=-2, cex.main=1.2)
  points(rows$cutoff, rows$nsubtrees, pch=19, cex=0.5, col=pal[i])
  points(marks[[g]][2], marks[[g]][1]+0.7, pch=25, bg='black', cex=0.9)
  i <- i+1
}
dev.off()


