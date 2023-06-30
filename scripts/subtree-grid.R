setwd("~/git/fluclades")
grid <- read.csv("results/subtree-grid.csv")
grid$st.diff <- log(grid$nsubtrees/18)

#plot(grid$minlen, grid$maxlen, cex=grid$mean.n.types)
pdf("~/papers/fluclades/subtree-grid.pdf", width=10, height=10)
par(mar=c(5,5,2,1), mfrow=c(2,2), cex=1)

plot(grid$minlen, grid$maxlen, pch=21,
     cex=sqrt(abs(grid$nsubtrees-18))/5,
     bg=ifelse(grid$nsubtrees > 18, 'pink', 'skyblue'),
     xlab="Minimum divergence", ylab="Maximum patristic distance")
title(main="Number of subtrees", adj=0, font.main=1)
text(x=-0.02, y=2.2, label="A", cex=2, xpd=NA)

plot(jitter(grid$maxlen), grid$mean.n.types, pch=19, cex=0.5,
     ylab="Mean number per subtree",
     xlab="Maximum patristic distance")
title(main="Number of different serotype labels", adj=0, font.main=1)
text(x=-0.35, y=9.1, label="B", cex=2, xpd=NA)

plot(grid$minlen, grid$maxlen, cex=2.5*sqrt(1-grid$nlabels), 
     pch=22, bg='red3', col=NA,
     xlab="Minimum divergence", ylab="Maximum patristic distance")
points(grid$minlen[grid$nlabels<1], grid$maxlen[grid$nlabels<1], 
       pch=22, cex=2.5, lwd=0.5, col='grey20')
title(main="Number of dropped labels", adj=0, font.main=1)
text(x=-0.02, y=2.2, label="C", cex=2, xpd=NA)

res <- read.csv("results/HA.mindiv0_08.maxpat1_2.subtypes.csv")
res <- res[res$serotype!='H192329',]
#res <- read.csv("results/HA.mindiv0.maxpat1_1.subtypes.csv")
res <- read.csv("results/HA.mindiv0_08.maxpat1.subtypes.csv")

par(mar=c(5,5,2,2))
plot(NA, xlim=c(0.75,17.25), ylim=c(-0.25, 16.25), xaxt='n', yaxt='n',
     xlab="Serotype labels", ylab="Subtree", bty='n')
title(main="Distribution of labels at y=1.2, v=0.08", adj=0, font.main=1)
text(x=-3.2, y=18., label="D", cex=2, xpd=NA)

#si <- c(2,0,3,4,1,6,7:16,17,5)  # yes this is ugly
for (i in 0:max(res$subtree)) {  # subtree
  #rows <- res[res$subtree==si[i+1],]
  rows <- res[res$subtree==i,]
  for (sero in rows$serotype) {  # serotype
    j <- as.integer(gsub("H([0-9]+)", "\\1", sero))
    count <- rows[rows$serotype==sero, 3]
    xx <- 1- (count / sum(rows$count))
    rect(j-1, i-1, j, i, col=rgb(xx, xx, xx), border=NA)
    if (xx >= 0.5) text(j-0.5, i-0.5, adj=0.5, label=count, cex=0.5)
  }
}
axis(side=1, at=1:18-0.5, label=paste("H", 1:18, sep=""), 
     cex.axis=0.8, las=2)
axis(side=2, at=0:17-0.5, label=si, las=2, cex.axis=0.8)
axis(side=4, at=0:17-0.5, 
     label=sapply(split(res$count, res$subtree), sum)[si+1],
     cex.axis=0.6, las=2, lwd=0, line=-0.5)
for (i in 0:18) {
  abline(v=i, col='grey80')
  abline(h=i, col='grey80')
}

dev.off()


if(FALSE) {
  require(ggfree)
  pdf("results/nsubtrees.pdf", width=5, height=5)
  par(mar=c(5,5,1,1))
  ridgeplot(split(log10(grid$nsubtrees), grid$minlen), fill=gg.rainbow(10, alpha=0.5),
            ylab="Minimum divergence", xlab="Log(number of subtrees)", xaxt='n',
            cex.axis=0.8, step=0.3)
  axis(side=1, at=0:4, labels=10^(0:4), cex.axis=0.8)
  abline(v=log10(18), lty=2)
  rug(log10(grid$nsubtrees))
  dev.off()
  
  boxplot(split(grid$nlabels, grid$minlen))  
}

