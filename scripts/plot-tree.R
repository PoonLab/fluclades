require(ggfree)
require(phangorn)

small <- read.tree("~/git/flu/data/removeX.n1000.nwk")
small <- ladderize(midpoint(small))

pdf("~/git/flu/small.pdf", height=100, width=30)
plot(tree.layout(small), cex=0.5, mar=c(0,0,0,20))
dev.off()

# all sequences

# repair metadata
phy0 <- read.tree("~/git/flu/data/removeX.ft2.nwk")
metadata <- read.csv("~/git/flu/data/genbank-metadata.csv")
names(metadata)[1] <- "locus"
metadata$accn <- gsub("^'([A-Z]+[0-9]+)\\..+", "\\1", phy0$tip.label)

phy <- read.tree("~/git/flu/data/removeX2.ft2.nwk")
phy <- midpoint(phy)
phy <- ladderize(phy)

pal <- gg.rainbow(18)

# it takes a couple of minutes to generate this layout
L <- tree.layout(phy)

# augment data frame with metadata
#metadata <- read.csv("~/git/flu/data/genbank-metadata.csv")
L$nodes$accn <- gsub("^'([A-Z]+[0-9]+)\\..+", "\\1", L$nodes$label)
idx <- match(L$nodes$accn, metadata$accn)
# record.name stores LOCUS, not ACCESSION
#idx <- which(as.character(L$nodes$accn[1:nrow(metadata)]) != as.character(metadata$accn))
#names(metadata)[1] <- "locus"
#metadata$accn <- L$nodes$accn[1:Ntip(phy)]

L$nodes$strain <- metadata$strain[idx]
L$nodes$serotype <- metadata$serotype[idx]
L$nodes$hsero <- gsub("^.*([Hh][0-9]+).*$", "\\1", L$nodes$serotype)
L$nodes$hsero <- toupper(L$nodes$hsero)
L$nodes$hsero[!grepl("^H[0-9]", L$nodes$hsero)] <- NA

L$nodes$host <- metadata$host[idx]

L$nodes$country <- metadata$country[idx]

L$nodes$coldate <- metadata$coldate[idx]
L$nodes$coldate <- as.Date(L$nodes$coldate, format="%d-%b-%Y")

rup <- function(x, n=2) {
  # runner up
  sort(x, decreasing=TRUE)[min(n, length(x))]
}

res <- 600
png("~/git/flu/bigtree2.png", width=8*res, height=11*res, res=res)
plot(L, label='n', lwd=1, mar=c(0,0,0,0))
add.scalebar(L, len=0.1)

for (i in 1:18) {
  print(i)
  # note phy$tip.label == L2$nodes$label[1:Ntip(phy)]
  #idx2 <- which(grepl(paste("\\(H", i, "[N)]", sep=""), phy$tip.label))
  idx <- which(L$nodes$hsero == paste("H", i, sep=""))
  root.idx <- draw.clade(L, tips=phy$tip.label[idx], col=pal[i], lwd=1)

  # label clades using tip farthest from origin
  last.tip <- idx[which.max(L$nodes$x[idx])]
  max.x <- rup(L$nodes$x[idx], 5) + (0.03*max(L$nodes$x))
  mid.y <- median(L$nodes$y[idx])
  #subtype <- gsub(".+\\((H[0-9]+).+", "\\1", L$nodes$label[last.tip])
  text(x=max.x, y=mid.y, label=L$nodes$hsero[last.tip], col=pal[i], xpd=NA)
  
  # mark tips from humans
  hum <- grepl("[Hh]omo|[Hh]uman", L$nodes$host)
  points(L$nodes$x[hum], L$nodes$y[hum], cex=0.15, pch=19)
}
dev.off()

# visualize/browse with Taxonium
write.tree(phy, file="temp.nwk")


# what are the nearest neighbours for unclassified sequences?
nearest.sero <- function(cidx, L, skip=0) {
  parent <- L$edges$parent[which(L$edges$child==cidx)]
  children <- setdiff(get.tips(parent, L), cidx)
  tab <- table(L$nodes$hsero[children])
  tries <- 0
  while (length(tab) == 0 || tries <= skip) {
    # all descendants of immediate ancestor are unknown
    if (tries > 20) {
      print(paste("Failed to retrieve labelled neighbour for ", cidx))
      break  # failsafe
    }
    parent <- L$edges$parent[which(L$edges$child==parent)]
    children <- setdiff(get.tips(parent, L), cidx)
    tab <- table(L$nodes$hsero[children])
    tries <- tries + 1
  }
  names(sort(tab, decreasing=TRUE))[1]
}

unknown <- which(L$nodes$n.tips==0 & is.na(L$nodes$hsero))
res <- sapply(unknown, function(x) nearest.sero(x, L, skip=1))


# which sequences might be mis-classified?
par(mar=c(5,5,1,1))
boxplot(split(L$nodes$y, L$nodes$hsero)) -> bx
miss <- which(L$nodes$y %in% bx$out)
predicted <- sapply(miss, function(cidx) {
  nearest.sero(cidx, L, skip=1)
})

toprint <- L$nodes[miss, c('accn', 'strain', 'hsero')]
toprint$predicted <- predicted
print(xtable(toprint), include.rownames=FALSE)

####################################
# unrooted layout

L2 <- tree.layout(phy, type='u')

res <- 600
png("~/git/flu/unrooted.png", width=10*res, height=10*res, res=res)
plot(L2, label='n')
for (i in 1:18) {
  print(i)
  # note phy$tip.label == L2$nodes$label[1:Ntip(phy)]
  idx <- which(grepl(paste("\\(H", i, "[N)]", sep=""), phy$tip.label))
  draw.clade(L2, tips=phy$tip.label[idx], col=pal[i])
  # label clades using tip farthest from origin
  last.tip <- idx[which.max(L2$nodes$r[idx])]
}
dev.off()

