#!/usr/bin/env python
#
# USAGE: python call_founders.py results.txt founders.txt
#
# where results.txt is the results file from determine_founders.py
#
# this script will aggregate calls in a region to guess at the most
# likely founder. for example, if there are 3 possible founders (ABC)
# at a few SNPs in sequence, followed by a single SNP for founder B, it
# will name the preceding 3 calls as B since that is the most
# parsimonious call to make.

import sys

founders=[]
canbe=None
possible=None
lastchr=''
lastcross=[1,1,1,1,1,]

all_lines = open(sys.argv[1]).read().split('\n')
all_lines = [line.split('\t') for line in all_lines[1:-1]]
all_lines.sort(key=lambda x:(x[5],int(x[6])))
out_sets = [[],[],[],[],[]]

fout = open(sys.argv[2], 'w')

for ln in range(1,len(all_lines)):
    d=all_lines[ln]
    pos = int(d[6])
    chr = d[5][3:]
    a = d[7].split(': ')
    a[1]=a[1].split(' // ')
    b = d[8].split(': ')
    if len(b)!=2:
        b.append([])
    else:
        b[1]=b[1].split(' // ')

    if len(founders)==0:
        founders=frozenset(a[1]+b[1])

    if lastchr!=chr:
        possible=['------',[],[],[],[]]

    for idx in range(1,5):
        if d[idx][0]==d[idx][1]: # homozygous
            if a[0][0]==d[idx][0]: # A allele
                canbe=frozenset(a[1])
            elif b[0][0]==d[idx][0]: # B allele
                canbe=frozenset(b[1])
        else: # heterozygous
            out_sets[idx].append('---')
            continue

        if len(possible[idx])==0 or possible[idx].isdisjoint(canbe):
            if lastcross[idx]!=0:
                for li in range(lastcross[idx], ln-1):
                    out_sets[idx][li-1] = out_sets[idx][ln-2]

            lastcross[idx]=ln
            possible[idx]=set(canbe)
        else:
            possible[idx].intersection_update(canbe)

        out_sets[idx].append(' // '.join(possible[idx]))

    lastchr=chr

for ln in range(0,len(all_lines)):
    d = all_lines[ln]
    print >>fout, '\t'.join( d[:7]+[out_sets[xi][ln-1] for xi in range(1,5)]+d[7:])
fout.close()
