setwd("~/git/fluclades")
grid <- read.csv("results/subtree-grid.csv")
grid$st.diff <- log(grid$nsubtrees/18)

#plot(grid$minlen, grid$maxlen, cex=grid$mean.n.types)
pdf("results/subtree-grid.pdf", width=15, height=5)
par(mar=c(5,5,2,1), mfrow=c(1,3), cex=1)

plot(grid$minlen, grid$maxlen, pch=21,
     cex=sqrt(abs(grid$nsubtrees-18))/5,
     bg=ifelse(grid$nsubtrees > 18, 'pink', 'skyblue'),
     xlab="Minimum divergence", ylab="Maximum patristic distance")
title(main="Number of subtrees", adj=0, font.main=1)

plot(jitter(grid$maxlen), grid$mean.n.types, pch=19, cex=0.5,
     ylab="Mean number per subtree",
     xlab="Maximum patristic distance")
title(main="Number of different serotype labels", adj=0, font.main=1)

plot(grid$minlen, grid$maxlen, cex=2.5*sqrt(1-grid$nlabels), 
     pch=22, bg='red3', col=NA,
     xlab="Minimum divergence", ylab="Maximum patristic distance")
points(grid$minlen[grid$nlabels<1], grid$maxlen[grid$nlabels<1], 
       pch=22, cex=2.5, lwd=0.5, col='grey20')
title(main="Number of dropped labels", adj=0, font.main=1)

dev.off()


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
