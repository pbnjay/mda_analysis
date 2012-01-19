#!/usr/bin/env Rscript --vanilla --no-restore
#
# USAGE: chr_plot.R datafile.txt output.pdf
#
# Feb. 14, 2011
#  - Original script by Daniel Gatti <Dan.Gatti@jax.org>
#
# Jan. 17, 2012
#  - Modified by Jeremy Jay <Jeremy.Jay@jax.org>
#    added ability to work on the commandline and accept generic arguments
#    simplified input file formats also
#
################################################################################
 
options(stringsAsFactors = F)
library(org.Mm.eg.db)
library(RColorBrewer)
 
args <- commandArgs(TRUE)

dataFile = args[1]
outputPDF = args[2]

################################################################################
founderSplit=" // "
#founderSplit=""

# format is SNPNAME001 CHR POS (founders-of-A) (founders-of-B) (founders-of-C) ...
thecalls = read.delim(dataFile)
# remove Y,M, relabel X, and scale positions
thecalls = thecalls[thecalls[,2] != "Y" & thecalls[,2] != "M",]
numchr = length(unique(thecalls[,2]))
thecalls[thecalls[,2] == "X",2] = numchr
thecalls[,3] = thecalls[,3] * 1e-6

chrnames = sort(as.numeric(unique(thecalls[,2])))
chrnames[ length(chrnames) ] = "X"

# Get the chromosome lengths.
chrlen = org.Mm.egCHRLENGTHS * 1e-6
chrlen = chrlen[names(chrlen) %in% chrnames]

# determine the array source names
arrayNames = colnames(thecalls)
arrayNames = arrayNames[4:length(arrayNames)]

# determine all the founder names
foundernames = c("-")
for(i in 1:length(arrayNames)) {
	mf = unlist(strsplit(thecalls[,3+i],founderSplit))
	foundernames = unique(c(foundernames,mf))
}

# assign colors to the founder names using a qualitative Brewer color scheme
states = aperm(as.array(brewer.pal(length(foundernames),"Paired")))
rownames(states) = foundernames

################################################################################
# Plot the genotypes using only the maximum probability at each SNP and arrange
# the colors to minimize jumping from one strand to the other.
# Arguments: thecalls: SNPID, chrN, position, call1, call2, ...
#						 whichcall: an integer 1-N where N is the number of calls in thecalls
plot.genotype.max2 = function(thecalls, whichcall, title, onScreen) {
  # Plot the chromosome skeletons.
  par(font = 2, font.axis = 2, font.lab = 2, cex.main=2, las = 1, mar=c(5,4,4,10)+0.1)
  plot(1, 1, col = 0, xlim = c(1, (2 * length(chrlen) + 1)),
       ylim = c(-3, 1.1 + max(chrlen)), xaxt = "n", xlab = "", ylab = "Mb", frame.plot=F)
  # Draw chr skeletons
  for(i in 1:length(chrlen)) {
    lines(2 * c(i, i), c(0, chrlen[i]))
  } # for(i)
  # Draw Mb lines.
  abline(h = 0:20 * 10, col = rgb(0.7,0.7,0.7))
  abline(h = 0:4  * 50, col = rgb(0.5,0.5,0.5), lwd = 1.2)
  text(2 * 1:length(chrlen), -5, names(chrlen), cex=2)
  title(title)

	par(xpd=T, font=1)
	toshow = which(rownames(states)!="-")
	legend((2 * length(chrlen) + 1)+2,1.1+max(chrlen),legend=rownames(states)[toshow], fill=states[toshow])
	par(xpd=F)
 
  #############

  offset = 0.6
  for(c in 1:length(chrlen)) {
    ss = which(thecalls[,2] == c)
    if(length(ss) > 0) {
			ms = t(matrix(unlist(strsplit(thecalls[ss,3+whichcall],founderSplit)),nr=2))
			leftfounders=ms[,1]
			rightfounders=ms[,2]
      y = thecalls[ss,3]  # chr position

      # xm is the middle of the chromosome.
      xm = 2 * c
      # xl is the left side.
      xl = xm - offset
      # xr is the right side.
      xr = xm + offset
      # Get the left and right colors.
      lcol = states[leftfounders]
      rcol = states[rightfounders]
      # Get the left and right side breaks points.
      lbreaks = c(1, which(lcol[1:(length(lcol)-1)] != lcol[2:length(lcol)]) +
                  1, length(lcol))
      rbreaks = c(1, which(rcol[1:(length(rcol)-1)] != rcol[2:length(rcol)]) +
                  1, length(rcol))
      lcol = lcol[lbreaks]
      rcol = rcol[rbreaks]
 
      # Make a data.frame of all breaks and try to minimize swapping colors
      # from one chr to the other.
      all.breaks = sort(union(lbreaks, rbreaks))
      all.breaks = data.frame(breaks = all.breaks,
									left = rep(NA, length(all.breaks)), 
									right = rep(NA, length(all.breaks)))
      all.breaks$left[match(lbreaks, all.breaks$breaks)]  = lcol
      all.breaks$right[match(rbreaks, all.breaks$breaks)] = rcol
 
			# Fill in any NA values.
			for(i in 1:nrow(all.breaks)) {
					 if(is.na(all.breaks$left[i])) {
					all.breaks$left[i] = all.breaks$left[i-1]
				} # if(is.na(all.breaks$left[i]))
 
				if(is.na(all.breaks$right[i])) {
					all.breaks$right[i] = all.breaks$right[i-1]
				} # if(is.na(all.breaks$right[i]))
			} #for(i)
 
			rect(xl, y[all.breaks$breaks[1]], xm, y[all.breaks$breaks[2]],
					col = all.breaks$left[1], density = NA, border =
					all.breaks$left[1])
 
			# Plot the right side of the chromosome.
			rect(xm, y[all.breaks$breaks[1]], xr, y[all.breaks$breaks[2]],
					 col = all.breaks$right[1], density = NA, border =
					 all.breaks$right[1])
 
			# Go through each break and see if the next color should be swapped.
			for(i in 2:nrow(all.breaks)) {
					 if(all.breaks$left[i] != all.breaks$right[i]) {
					if(all.breaks$left[i] == all.breaks$right[i-1] |
						 all.breaks$left[i-1] == all.breaks$right[i]) {
							 tmp = all.breaks$left[i]
							 all.breaks$left[i] = all.breaks$right[i]
							 all.breaks$right[i] = tmp
					} # if(all.breaks$left[i] == all.breaks$right[i-1] | ...
					 } # if(all.breaks$left[i] != all.breaks$right[i])
					
			rect(xl, y[all.breaks$breaks[i]], xm, y[all.breaks$breaks[i+1]],
					col = all.breaks$left[i], density = NA, border =
					all.breaks$left[i])
 
			# Plot the right side of the chromosome.
					 rect(xm, y[all.breaks$breaks[i]], xr, y[all.breaks$breaks[i+1]],
								col = all.breaks$right[i], density = NA, border =
								all.breaks$right[i])
			} # for(i)
    } # if(length(ss) > 0)
  } # for(c)
} # plot.genotype.max2()

# draw the desired plot into a multi-page PDF
pdf(file=outputPDF, height=10, width=(length(chrlen)*0.75 + 3))
for(i in 1:length(arrayNames)) {
	plot.genotype.max2(thecalls, i, arrayNames[i], FALSE)
}
dev.off()
