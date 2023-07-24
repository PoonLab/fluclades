setwd("~/papers/fluclades")

# these metadata cannot be released to public domain
gisaid <- read.csv("gisaid.csv")

## WHO FluNet database 
#fnt <- read.csv("VIW_FNT.csv")
#viw <- read.csv("VIW_FID_EPI.csv")
#viw <- viw[viw$AGEGROUP_CODE=='All',]
#ncases <- sapply(split(viw$ILI_CASES, viw$MMWR_YEAR), function(x) sum(x, na.rm=T))

## number of influenza A detections (all subtypes)
#flua <- sapply(split(fnt$INF_A, fnt$MMWR_YEAR), function(x) sum(x, na.rm=T))
#nsamp <- sapply(split(fnt$SPEC_RECEIVED_NB, fnt$MMWR_YEAR), function(x) sum(x, na.rm=T))
#idx <- match(gisaid$year, names(flua))
#gisaid$ndetect <- flua[idx]

pdf("~/papers/fluclades/gisaid-nseqs.pdf", width=4.5, height=4.5)

par(mar=c(5,5,1,1))
barplot(gisaid$nseq[order(gisaid$year)]/1000, space=0, 
        border=NA, col='grey50', cex.axis=0.8, las=2,
        ylab="Number of IAV sequences (thousands)")
p <- pretty(1:nrow(gisaid))
years <- sort(unique(gisaid$year))
axis(side=1, at=p+0.5, labels=years[p+1], cex.axis=0.8)
title(xlab="Year of sample collection")

dev.off()

# evidence of exponential trend
fit <- lm(log(nseq)~year, data=gisaid)
summary(fit)  # displays R-squared

