#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) != 2) {
    stop("Usage: midpoint.R <input NWK> <output NWK>")
}
require(phangorn)
phy <- read.tree(args[1])
phy <- midpoint(phy)
write.tree(phy, args[2])

