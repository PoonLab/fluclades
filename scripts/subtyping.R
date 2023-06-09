require(phangorn)
setwd("~/git/fluclades/data")

# filtered = removed patent and vaccine sequences
tfiles <- Sys.glob("gb-filtered-*.fa.nwk")
for (tf in tfiles) {
  phy <- read.tree(tf)
  of <- gsub("\\.fa\\.nwk", ".mid.nwk", tf)
  write.tree(midpoint(phy), of)
}

phy <- read.tree("HA.nwk")
#phy <- midpoint(phy)
phy <- reorder(phy, "postorder")

get.path <- function(phy, child, res=c()) {
  # returns sequence of node indices from tip to root
  parent <- phy$edge[which(phy$edge[,2] == child), 1]
  #res <- c(res, parent)
  res <- c(child, res)
  if (parent == Ntip(phy)+1) {
    return(res)  # root
  }
  return(get.path(phy, parent, res))
}

df <- data.frame(parent=phy$edge[,1], child=phy$edge[,2], 
                 len=phy$edge.length, label=NA, 
                 total.len=0, mean.pl=0, n.tips=0)
df$label[match(1:Ntip(phy), df$child)] <- gsub("^'|'$", "", phy$tip.label)

for (i in 1:nrow(df)) {
  if (i%%100 == 0) { print(i) }
  if (df$child[i] <= Ntip(phy)) {
    df$n.tips[i] <- 1
  }
  p.row <- which(df$child==df$parent[i])
  df$n.tips[p.row] <- df$n.tips[p.row] + df$n.tips[i]
  df$total.len[p.row] <- df$total.len[p.row] + df$len[i] + df$total.len[i]
  df$mean.pl[p.row] <- df$mean.pl[p.row] + df$len[i]*df$n.tips[i] + df$mean.pl[i]
}
df$mean.pl <- df$mean.pl / df$n.tips
df$mean.bl <- df$total.len / df$n.tips

#write.csv(df, file="HA.subtyping.csv")
df <- read.csv("HA.subtyping.csv")

# visualize to select criteria?
plot(df$len[df$n.tips>1], df$mean.bl[df$n.tips>1])

path <- get.path(phy, 1)
plot(1:length(path), df$len[match(path, df$child)], type='l', xlim=c(1, 30))
for (i in sample(1:Ntip(phy), 10)) {
  path <- get.path(phy, i)
  lines(1:length(path), df$len[match(path, df$child)])
}

plot(NA, xlim=c(0, 1), ylim=c(0, 1), xlab="Normalized depth", 
     ylab="Scaled mean path length")
for (i in sample(1:Ntip(phy), 100)) {
  path <- get.path(phy, i)
  tot.len <- sum(df$len[match(path, df$child)])
  lines(1:length(path)/length(path), df$mean.pl[match(path, df$child)] / tot.len)
}


# apply search criteria to extract subtrees
subtree.search <- function(df, child=Ntip(phy)+1, cutoff=0.15, min.len=0.01, 
                           depth=0, max.depth=10, res=c()) {
  if (child %in% df$child && 
      df$mean.pl[df$child==child] < cutoff && 
      df$len[df$child==child] >= min.len) {
    # mean path length in clade above this node falls below cutoff
    return(c(res, child))
  }
  else {
    children <- df$child[df$parent==child]
    depth <- depth+1
    if (depth > max.depth) {
      return(res)
    }
    for (child in children) {
      res <- subtree.search(df, child, cutoff=cutoff, min.len=min.len, 
                            depth=depth, max.depth=max.depth, res=res)
    }
    return(res)
  }
}

# returns indices of internal nodes selected by search criteria
res <- subtree.search(df, cutoff=0.3, min.len=-1, max.depth=100)

get.tips <- function(df, parent, res=c()) {
  children <- df$child[df$parent==parent]
  for (child in children) {
    if (df$n.tips[df$child==child] == 1) {
      res <- c(res, child)
    }
    else {
      res <- get.tips(df, child, res)  
    }
  }
  return(res)
}

gsub2 <- function(pattern, replacement, x) {
  # gsub returns full string if pattern does not match
  idx <- !grepl(pattern, x)
  x <- gsub(pattern, "\\1", x)
  x[idx] <- NA
  x
}
df$gb.sero <- gsub2("^.+[^A-Z](H[0-9]{1,2}).+$", "\\1", df$label)

# this step takes too long - Python would be faster
subtypes <- lapply(res, function(i) {
  tips <- get.tips(df, i)
  return(df$gb.sero[df$child %in% tips])
})
lapply(subtypes, table)


# add subtype annotation from Genbank metadata?
metadata <- read.csv("gb-filtered-HA.metadata.csv")
metadata$subtype <- gsub("(H[X0-9]+).*", "\\1", toupper(metadata$serotype))
df$accn <- sapply(df$label, function(x) strsplit(x, "_")[[1]][1])
df$accn <- gsub("\\.[1-9]$", "", df$accn)
df$meta.subtype <- metadata$subtype[match(df$accn, metadata$accn)]




