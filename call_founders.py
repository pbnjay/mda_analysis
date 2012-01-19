#!/usr/bin/env python
#
# USAGE: python call_founders.py possible_founders.txt founders.txt
#
# where possible_founders.txt is the results file from determine_founders.py
#
# this script will aggregate calls in a region to guess at the most
# likely founder strain. for example, if there are 3 possible founders
# (ABC) at a few SNPs in sequence, followed by a single SNP for founder
# B, it will name the preceding calls as B since that is the most
# parsimonious call to make.

import sys

founders=set()

# read in the file, split on tabs, and sort by chr/pos
all_lines = open(sys.argv[1]).read().split('\n')
header = all_lines[0].split('\t')
all_lines = [line.split('\t') for line in all_lines[1:-1]]
all_lines.sort(key=lambda x:(x[-4],int(x[-3])))

num_arrays = len(all_lines[0])-5
out_distributions =[{} for x in range(0,num_arrays+1)]
out_sets = [[] for x in range(0,num_arrays+1)]
lastcross=[[1,1] for x in range(0,num_arrays+1)]
lastchr=''
canbe=None
possible=None

fout = open(sys.argv[2], 'w')
print >>fout, "%s\t%s\t%s\t%s" % ( header[0], header[-4], header[-3], '\t'.join(header[1:num_arrays+1]))
for ln in range(1,len(all_lines)):
    d=all_lines[ln]
    chr = d[-4][3:]
    pos = d[-3]
    a = d[-2].split(': ')
    a[1]=set(a[1].split(' // '))
    b = d[-1].split(': ')
    b[1]=set(b[1].split(' // '))

    if '' in a[1]:
        a[1].remove('')
    if '' in b[1]:
        b[1].remove('')
    if len(a[1])==0 or len(b[1])==0:
        continue

    founders.update( a[1] )
    founders.update( b[1] )

    if lastchr!=chr:
        possible=['-']
        for idx in range(1,num_arrays+1):
            possible.append( [set(),set()] )

    out_sets[0].append( [d[0], chr, pos] )
    for idx in range(1,num_arrays+1):
        left=set()
        right=set()
        left_could_be=0
        right_could_be=0
        if d[idx][0]==a[0][0]: # A allele
            left=a[1]

            if d[idx][1]==a[0][0]: # homozygous A
                right=a[1]
            elif d[idx][1]==b[0][0]: # het
                right=b[1]

        elif d[idx][0]==b[0][0]: # B allele
            left=b[1]

            if d[idx][1]==b[0][0]: # homozygous B
                right=b[1]
            elif d[idx][1]==a[0][0]: # het
                right=a[1]

        crossover=None
        if not possible[idx][0].isdisjoint(left) and not possible[idx][1].isdisjoint(right):
            possible[idx][0].intersection_update(left)
            possible[idx][1].intersection_update(right)
        elif not possible[idx][1].isdisjoint(left):
            possible[idx][1].intersection_update(left)
            possible[idx][0]=set(right)
            crossover=[0]
        elif not possible[idx][0].isdisjoint(right):
            possible[idx][0].intersection_update(right)
            possible[idx][1]=set(left)
            crossover=[1]
        else:
            possible[idx][0]=set(left)
            possible[idx][1]=set(right)
            crossover=[0,1]

        if crossover is not None:
            for i in crossover:
                start = lastcross[idx][i]
                lastcross[idx][i] = len(out_sets[0])
                # copy this set up to the last crossover point
                for li in range(start, lastcross[idx][i]):
                    out_sets[idx][li-1][i] = out_sets[idx][-1][i]

        out_sets[idx].append( [set(possible[idx][0]),set(possible[idx][1])] )

        if len(possible[idx][0])==1:
            for pf in possible[idx][0]:
                if pf not in out_distributions[idx]:
                    out_distributions[idx][ pf ]=1
                else:
                    out_distributions[idx][ pf ]+=1
        if len(possible[idx][1])==1:
            for pf in possible[idx][1]:
                if pf not in out_distributions[idx]:
                    out_distributions[idx][ pf ]=1
                else:
                    out_distributions[idx][ pf ]+=1

    lastchr=chr

for ln in range(0,len(out_sets[0])):
    no_out=False
    founderstates=[]
    for xi in range(1,num_arrays+1):
        ostr=''
        for i in [0,1]:
            if len(out_sets[xi][ln][i])==0:
                no_out=True
                break

            if len(out_sets[xi][ln][i])==1:
                for oset in out_sets[xi][ln][i]:
                    ostr+=oset+" // "
                continue

            ostr='- // - // '
            break

        if no_out:
            break

        founderstates.append(ostr[:-4])

    if no_out:
        continue

    print >>fout, "%s\t%s" % ('\t'.join(out_sets[0][ln]), '\t'.join(founderstates))
fout.close()
