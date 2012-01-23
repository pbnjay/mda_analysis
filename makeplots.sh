#!/bin/bash
#
# this script will make a founder haplotype plot given the proper inputs:
#   genotypes.txt - SNP A/B calls in the format produced by get_genotypes.R
#   snpdownload.csv - SNP-to-founder mapping (see README for source)
#   celfiledir - if provided, will call get_genotypes.R for you
#   threshold - an optional confidence threshold to pass along to get_genotypes.R

# determine the right dir for this script

ODIR=`pwd`
CURDIR=`dirname $0`
cd $CURDIR
CURDIR=`pwd`
cd $ODIR

if [ -d "$3" ]; then
	echo get_genotypes.R "$3" "$1" $4
	$CURDIR/get_genotypes.R "$3" "$1" $4
fi

if [ ! -f "$1" ]; then
	echo "USAGE: makeplots.sh genotypes.txt snpdownload.csv [celfiledir [threshold]]"
	exit
elif [ ! -f "$2" ]; then
	echo "USAGE: makeplots.sh genotypes.txt snpdownload.csv [celfiledir [threshold]]"
	exit
fi

GENOS=`basename "$1" .txt`
echo determine_founders.py "$2" "$GENOS.txt" "$GENOS.possible.txt"
$CURDIR/determine_founders.py "$2" "$GENOS.txt" "$GENOS.possible.txt"

echo call_founders6.py "$GENOS.possible.txt" "$GENOS.founders.txt"
$CURDIR/call_founders.py "$GENOS.possible.txt" "$GENOS.founders.txt" "$GENOS.accuracy.txt"

echo chr_plot.R "color_defs.txt" "$GENOS.founders.txt" "$GENOS"
$CURDIR/chr_plot.R "$GENOS.founders.txt" "$GENOS"

echo "You should have a founder haplotype plot in $GENOS.pdf"
