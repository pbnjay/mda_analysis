#!/usr/bin/env Rscript --vanilla --no-restore
#
# USAGE: chr_plot.R datafile.txt output.pdf [qtls.txt]
#
# Feb. 14, 2011
#  - Based on plot.genotype.max2 by Daniel Gatti <Dan.Gatti@jax.org>
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

qtlFile = NA
if(length(args)>=3) {
  qtlFile = args[3]
}

################################################################################
ignoreChrs=c("Y","M")
onlyQTLChrs=TRUE
fillChrEnds=TRUE
mergeUncalled=TRUE
chromoSpacing=6
founderSplit=" // "
#founderSplit=""

# format is SNPNAME001 CHR POS (founders-of-A) (founders-of-B) (founders-of-C) ...
thecalls = read.delim(dataFile)
thecalls[,3] = thecalls[,3] * 1e-6

# remove ignored chromosomes
thecalls = thecalls[thecalls[,2] %in% ignoreChrs == FALSE,]

# Get the chromosome lengths.
chrlen = org.Mm.egCHRLENGTHS * 1e-6
remove = grep("random", names(chrlen))
chrlen = chrlen[-remove]
chrlen = chrlen[names(chrlen) %in% ignoreChrs==FALSE]

# determine the array source names
arrayNames = colnames(thecalls)
arrayNames = arrayNames[4:length(arrayNames)]

# determine all the founder names
foundernames = c("-")
for(i in 1:length(arrayNames)) {
  mf = unlist(strsplit(thecalls[,3+i],founderSplit))
  foundernames = unique(c(foundernames,mf))
}
foundernames = sort(foundernames)

# assign colors to the founder names using a qualitative Brewer color scheme
if( mergeUncalled ) {
  states = aperm(as.array(c(NA, brewer.pal(length(foundernames)-1,"Paired"))))
} else {
  states = aperm(as.array(brewer.pal(length(foundernames),"Paired")))
}
rownames(states) = foundernames

if( !is.na(qtlFile) ) {
  qtls = read.delim(qtlFile)
  qtls = qtls[qtls[,2] %in% ignoreChrs == FALSE,]
  qtls[,3:5] = qtls[,3:5] * 1e-6
}

################################################################################
# Plot the genotypes using only the maximum probability at each SNP and arrange
# the colors to minimize jumping from one strand to the other.
# Arguments: thecalls: SNPID, chrN, position, call1, call2, ...
#             whichcall: an integer 1-N where N is the number of calls in thecalls
plot.genotype.max2 = function(thecalls, arrayPair) {
  par(font = 2, font.axis = 2, font.lab = 2, cex.main=1.5, las = 1, mar=c(2,5,4,10)+0.1)
  offset = chromoSpacing/4
  chrs = chrlen

  if( length(arrayPair)>1 ) {
    offset = chromoSpacing/10
    if( onlyQTLChrs && !is.na(qtlFile) ) {
      chrs = chrs[names(chrs) %in% qtls[,2]]
    }
    plot(1, 1, col = 0, xlim = c(chromoSpacing-(offset*2), (chromoSpacing * length(chrs))+(offset*2)),
         ylim = c(-3, 1.1 + max(chrs)), xaxt = "n", xlab = "", ylab = "Mb", frame.plot=F)
  } else {
    plot(1, 1, col = 0, xlim = c(chromoSpacing, (chromoSpacing * length(chrs))),
         ylim = c(-3, 1.1 + max(chrs)), xaxt = "n", xlab = "", ylab = "Mb", frame.plot=F)
  }

  # Draw chr skeletons
  #for(i in 1:length(chrs)) {
  #  lines(chromoSpacing * c(i, i), c(0, chrs[i]))
  #} # for(i)

  # Draw Mb lines.
  abline(h = 0:20 * 10, col = rgb(0.9,0.9,0.9))
  abline(h = 0:4  * 50, col = rgb(0.7,0.7,0.7), lwd = 1.2)
  text(chromoSpacing * 1:length(chrs), -5, names(chrs), cex=2)

  # title the plot
  if( length(arrayPair)>1 ) {
    plot_title <- paste(arrayNames[arrayPair][1], "(left) vs", arrayNames[arrayPair][2], "(right)")
  } else {
    plot_title <- arrayNames[arrayPair[1]]
  }
  #title(main="Founder Map", sub=plot_title, col.sub="black")
  title(main=plot_title)

  pb <- par("usr")
  par(xpd=NA, font=1)
  toshow = which(rownames(states)!="-")
  legend( pb[2]+1, pb[4], legend=rownames(states)[toshow], fill=states[toshow])
  par(xpd=F)
 
  #############
  if( length(arrayPair)>1 && !is.na(qtlFile) ) {
    for(c in 1:length(chrs)) {
      cn = names(chrs)[c]

      # draw the qtls
      qs = which(qtls[,2] == cn)
      if(length(qs)>0) {
        qoff <- offset*3.5
        # draw span in red (brightness depends on % length)
        rect(chromoSpacing*c-qoff, qtls[qs,4], chromoSpacing*c+qoff, qtls[qs,5], col = rgb(1.0-((qtls[qs,5]-qtls[qs,4])/chrs[c]),0,0), border=NA)
        # draw peak in green
        rect(chromoSpacing*c-qoff, qtls[qs,3]-0.25, chromoSpacing*c+qoff, qtls[qs,3]+0.25, col = rgb(0,1.0,0), border=NA)
      }
    }
  }

  for(wap in 1:length(arrayPair)) {
    whichcall <- arrayPair[wap]

    for(c in 1:length(chrs)) {
      cn = names(chrs)[c]

      # draw the founder haplotypes
      ss = which(thecalls[,2] == cn)
      if(length(ss) > 0) {
        ms = t(matrix(unlist(strsplit(thecalls[ss,3+whichcall],founderSplit)),nr=2))
        leftfounders=ms[,1]
        rightfounders=ms[,2]
        y = thecalls[ss,3]  # chr position

        # xm is the middle of the chromosome.
        xm = chromoSpacing * c
        if( length(arrayPair)>1 && wap==1 ) {
          xm = xm-(offset*2)
        } else if( length(arrayPair)>1 ) {
          xm = xm+(offset*2)
        }
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
        for(i in 2:nrow(all.breaks)) {
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
          if( is.na(all.breaks$left[i]) || is.na(all.breaks$right[i]) || is.na(all.breaks$left[i-1]) || is.na(all.breaks$right[i-1])) {}
          else if( all.breaks$left[i] != all.breaks$right[i]) {
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

        if( fillChrEnds ) {
          rect(xl, 1, xm, y[all.breaks$breaks[1]],
            col = all.breaks$left[1], density = NA, border =
            all.breaks$left[1])

          # Plot the right side of the chromosome.
          rect(xm, 1, xr, y[all.breaks$breaks[1]],
            col = all.breaks$right[1], density = NA, border =
            all.breaks$right[1])

          nr = nrow(all.breaks)
          rect(xl, y[all.breaks$breaks[nr]], xm, chrs[c],
            col = all.breaks$left[nr], density = NA, border =
            all.breaks$left[nr])

          # Plot the right side of the chromosome.
          rect(xm, y[all.breaks$breaks[nr]], xr, chrs[c],
            col = all.breaks$right[nr], density = NA, border =
            all.breaks$right[nr])
        }
      } # if(length(ss) > 0)
    } # for(c)
  } # for(wap)
} # plot.genotype.max2()

# draw the desired plot into a multi-page PDF
pdf(file=outputPDF, height=8.5, width=length(chrlen)+3)
for(i in 1:length(arrayNames)) {
  for(j in 1:length(arrayNames)) {
    if( i<j ) {
      arrayPair=c(i,j)
      plot.genotype.max2(thecalls, arrayPair)
    }
  }
}
for(i in 1:length(arrayNames)) {
  plot.genotype.max2(thecalls, c(i))
}
dev.off()
