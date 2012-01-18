#!/usr/bin/env Rscript --vanilla --no-restore
#
# USAGE: chr_plot.R color_defs.txt datafile.txt whichArray output.pdf "plot title"
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
 
args <- commandArgs(TRUE)

colorsFile = args[1]
dataFile = args[2]

whichArray = as.numeric(args[3])
outputPDF = args[4]

################################################################################
# Get the chromosome lengths.
chrlen = org.Mm.egCHRLENGTHS
remove = grep("random", names(chrlen))
if(length(remove) > 0) {
  chrlen = chrlen[-remove]
} # if(length(remove) > 0)
chrlen = chrlen * 1e-6
 
# format is AB #000000 #000000
states = read.delim(colorsFile)
rownames(states) = states[,1]

# format is SNPNAME001 CHR POS XY AB CD ...
thecalls = read.delim(dataFile, header=FALSE)
# remove Y,M, relabel X, and scale positions
thecalls = thecalls[thecalls[,2] != "Y" & thecalls[,2] != "M",]
thecalls[thecalls[,2] == "X",2] = "20"
thecalls[,3] = thecalls[,3] * 1e-6

#############################################################################
# Draw the skeleton of the chromosomes.
# Arguments: chr: vector with chr names to plot.
#            chrlen: vector with chr lengths, named with chr names.
#            onScreen: logical, T if plot is on screen, false if to file.
plot.chr.skeletons = function(chr, chrlen, onScreen = T) {
  chrlen = chrlen[names(chrlen) %in% chr]
  # Draw the plotting area.
  if(!onScreen) {
    par(cex = 2)
  } # if(!onScreen)
  par(font = 2, font.axis = 2, font.lab = 2, las = 1)
  plot(1, 1, col = 0, xlim = c(1, (2 * length(chrlen) + 1)),
       ylim = c(-3, 1.1 + max(chrlen)), xaxt = "n", xlab = "", ylab = "Mb")
  # Draw chr skeletons
  for(i in 1:length(chrlen)) {
    lines(2 * c(i, i), c(0, chrlen[i]))
  } # for(i)
  # Draw Mb lines.
  abline(h = 0:20 * 10, col = rgb(0.7,0.7,0.7))
  abline(h = 0:4  * 50, col = rgb(0.5,0.5,0.5), lwd = 1.2)
  text(2 * 1:length(chrlen), -5, names(chrlen))
} # plot.chr.skeletons()
 
 
################################################################################
# Plot the genotypes using only the maximum probability at each SNP and arrange
# the colors to minimize jumping from one strand to the other.
# Arguments: thecalls: SNPID, chrN, position, call1, call2, ...
#						 whichcall: an integer 1-N where N is the number of calls in thecalls
plot.genotype.max2 = function(thecalls, whichcall, title, onScreen) {
  # Plot the chromosome skeletons.
  chr = sort(as.numeric(unique(thecalls[,2])))
  chr[20] = "X"
  plot.chr.skeletons(chr, chrlen, onScreen)
  title(title)
 
  offset = 0.6
  for(c in 1:length(chrlen)) {
    ss = which(thecalls[,2] == c)
    if(length(ss) > 0) {
      ms = thecalls[ss,3+whichcall]  # should be the founder call ex: "CD"
      y = thecalls[ss,3]  # chr position

      # xm is the middle of the chromosome.
      xm = 2 * c
      # xl is the left side.
      xl = xm - offset
      # xr is the right side.
      xr = xm + offset
      # Get the left and right colors.
      lcol = states[ms,2]
      rcol = states[ms,3]
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

# draw the desired plot
pdf(file=outputPDF, height=8.5, width=11)
plot.genotype.max2(thecalls, whichArray, args[5], FALSE)
