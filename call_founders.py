#!/usr/bin/env python
#
# USAGE: python call_founders.py results.txt founders.txt
#
# where results.txt is the results file from determine_founders.py
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
all_lines = [line.split('\t') for line in all_lines[1:-1]]
all_lines.sort(key=lambda x:(x[5],int(x[6])))

num_arrays = len(all_lines[0])-5
out_distributions =[{} for x in range(0,num_arrays+1)]
out_sets = [[] for x in range(0,num_arrays+1)]
lastcross=[1 for x in range(0,num_arrays+1)]
lastchr=''
canbe=None
possible=None

fout = open(sys.argv[2], 'w')

for ln in range(1,len(all_lines)):
    d=all_lines[ln]
    chr = d[-4][3:]
    #pos = int(d[-3])
    a = d[-2].split(': ')
    a[1]=a[1].split(' // ')
    b = d[-1].split(': ')
    b[1]=b[1].split(' // ')

    founders.update( a[1] )
    founders.update( b[1] )

    if lastchr!=chr:
        possible=['-']
        for idx in range(1,num_arrays+1):
            possible.append( set() )

    for idx in range(1,num_arrays+1):
        if d[idx][0]==d[idx][1]: # homozygous
            if a[0][0]==d[idx][0]: # A allele
                canbe=frozenset(a[1])
            elif b[0][0]==d[idx][0]: # B allele
                canbe=frozenset(b[1])

        else: # heterozygous
            canbe=frozenset(a[1]+b[1])
            #out_sets[idx].append(set())
            #continue

        if len(possible[idx])==0 or possible[idx].isdisjoint(canbe):
            if lastcross[idx]!=0:
                for li in range(lastcross[idx], ln-1):
                    out_sets[idx][li-1] = out_sets[idx][ln-2]

            lastcross[idx]=ln
            possible[idx]=set(canbe)
        else:
            possible[idx].intersection_update(canbe)

        out_sets[idx].append(possible[idx])

        if len(possible[idx])==1:
            for pf in possible[idx]:
                if pf not in out_distributions[idx]:
                    out_distributions[idx][ pf ]=1
                else:
                    out_distributions[idx][ pf ]+=1

    lastchr=chr

foundermap={}
fstr="ABCDEFGHIJKLMNOPQRSTUV"
for f in founders:
    if f=='':
        continue
    foundermap[f]=fstr[len(foundermap)]

print foundermap

for ln in range(0,len(all_lines)):
    d = all_lines[ln]

    no_out=False
    founderstates=[]
    for xi in range(1,num_arrays+1):
        if '-' in out_sets[xi][ln-1]:
            out_sets[xi][ln-1].remove('-')
        if '' in out_sets[xi][ln-1]:
            out_sets[xi][ln-1].remove('')

        if len(out_sets[xi][ln-1])==0:
            no_out=True
            break

        if len(out_sets[xi][ln-1])==1:
            ostr=''
            for oset in out_sets[xi][ln-1]:
                ostr+=foundermap[ oset ]
            ostr+=ostr
            founderstates.append( ostr )
        else:
            no_out=True
            continue

            osets = list(out_sets[xi][ln-1])
            osets.sort(key=lambda x:out_distributions[xi][x])

            ostr = foundermap[ osets[0] ]
            ostr+= foundermap[ osets[1] ]
            if ostr[0]>ostr[1]:
                ostr = ostr[::-1]
            founderstates.append( ostr )

    if no_out:
        continue

    #print >>fout, "%s\t%s\t%s\t%s" % (d[0],d[5][3:], d[6], '\t'.join([' // '.join(out_sets[xi][ln-1]) for xi in range(1,5)]))
    print >>fout, "%s\t%s\t%s\t%s" % (d[0],d[5][3:], d[6], '\t'.join(founderstates))
fout.close()
