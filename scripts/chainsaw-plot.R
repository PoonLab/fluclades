setwd("~/git/fluclades")

# manual runs of chainsaw.py
pdf("results/chainsaw.pdf", width=4.5, height=4.5)
par(mar=c(5,5,1,1))
x <- c(0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25, 0.30)
y <- c(  64,   43,   29,   25,   22,   21,   18,   16,   15,   13,   13,   12,    9)
plot(x, y, type='s', xlab="Maximum internal branch length", 
     ylab="Number of subtrees", bty='n')
abline(h=18, lty=2)
dev.off()


labels <- read.csv("results/chainsaw-0.15.full.csv")
pat <- ".+_(H[0-9]+)N*[0-9]*_.+"
labels$serotype <- gsub(pat, "\\1", labels$tip.label)
labels$serotype[!grepl(pat, labels$tip.label)] <- NA

table(labels$subtree, labels$serotype)
