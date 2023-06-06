require(phangorn)
setwd("~/git/flu/data")

tfiles <- Sys.glob("gb-filtered-*.fa.nwk")
for (tf in tfiles) {
  phy <- read.tree(tf)
  of <- gsub("\\.fa\\.nwk", ".mid.nwk", tf)
  write.tree(midpoint(phy), of)
}

phy <- read.tree("gb-filtered-HA.mid.nwk")
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

#write.csv(df, file="gb-filtered-HA.subtyping.csv")
df <- read.csv("gb-filtered-HA.subtyping.csv")

# visualize to select criteria?
plot(df$len[df$n.tips>1], df$mean.pl[df$n.tips>1], 
     cex=sqrt(df$n.tips[df$n.tips>1])/10)

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
subtree.search <- function(df, child=Ntip(phy)+1, cutoff=0.15, min.len=0.05, res=c()) {
  if (child %in% df$child && df$mean.pl[df$child==child] < cutoff) {
    return(c(res, child))
  }
  else {
    children <- df$child[df$parent==child]
    for (child in children) {
      res <- subtree.search(df, child, cutoff=cutoff, res=res)
    }
    return(res)
  }
}

res <- subtree.search(df)

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

subtypes <- lapply(res, function(i) {
  tips <- get.tips(df, i)
  return(df$gb.sero[df$child %in% tips])
})
lapply(subtypes, table)

