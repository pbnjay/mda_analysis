#!/usr/bin/env Rscript --vanilla --no-restore
#
# USAGE: get_genotypes.R celfiles/ genotypes.txt
#
# this is pretty much straight out of the MouseDivGeno vignette:
#  http://cgd.jax.org/tools/mousedivgeno/MouseDivGeno.pdf

args <- commandArgs(TRUE)

library("MouseDivGeno")
load("MouseDivData.RData")

genoVinoResult <- mouseDivGenotypeCEL(
		snpProbeInfo = snpProbeInfo,
		snpInfo = snpInfo,
		referenceDistribution = snpReferenceDistribution,
		chromosomes = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, "X", "Y", "M"),
		celFiles = args[1],
		confScoreThreshold = 0)

write.csv(genoVinoResult$geno, args[2], quote=FALSE)

# get confidence values for every call
#write.csv(genoVinoResult$conf, args[2], quote=FALSE)

# get annotated results
#infoMatched <- snpInfo[match(rownames(genoVinoResult), snpInfo$snpId), ]
#annoResults <- cbind(infoMatched, annoResults)
#write.csv(annoResults[1, , drop=FALSE], args[2], quote=FALSE)
