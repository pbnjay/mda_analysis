#!/usr/bin/env python
#
# USAGE: python call_founders.py possible_founders.txt founders.txt founder_accuracy.txt
#
# where possible_founders.txt is the results file from determine_founders.py
#
# basic idea here is to walk down the chromosome, and if this snp and
# the previous snp share possible founders, intersection them and save
# the intersection to both snps. after doing this first in the forward
# direction, repeat it in the backward direction. (if you ever come to a
# heterozygous call, skip it).
#
# In a second phase, repeat the forward walk again, this time counting
# SNPs with a single founder, aggregated by founder. (for this we use a
# sliding window of N=1500 SNPs) When the walk reaches a SNP which still
# have an ambiguous founder (ie 3 possible founders), find the founder
# which has been seen the most in the last N singleton-founder SNPs, and
# save it. Repeat this process in the backwards direction, replacing
# "best" founders if they have been seen in more of the last N SNPs than
# the forward direction. Then, walk through the saved maps and replace
# the founder calls.
#
# Once these is completed, walk down the chromosome again, and output any
# (homozygous) snps that have a single founder attributed. Also
# calculate the error rate of the determinations.
#

from collections import deque
import sys

windowsize=25
previous_max=1500

# read in the file and parse header
f = open(sys.argv[1])
header = f.next().strip().split('\t')
num_arrays = len(header)-5

chrsnps = {}
allsnps = {}
for line in f:
    d = line.strip().split('\t')
    if d[0]=='':
        continue
    chr = d[-4][3:]
    pos = int(d[-3])
    a = d[-2].split(': ')
    b = d[-1].split(': ')
    if len(a)==2:
        a[1]=set(a[1].split(' // '))
        if '' in a[1]:
            a[1].remove('')
        a[1] = frozenset(a[1])
    else:
        a.append(frozenset())

    if len(b)==2:
        b[1]=set(b[1].split(' // '))
        if '' in b[1]:
            b[1].remove('')
        b[1] = frozenset(b[1])
    else:
        b.append(frozenset())

    hasmatch = False
    snpinfo = {'snp': d[0], 'chr': chr, 'pos': pos, 'calls': [[] for x in range(0,num_arrays)] }
    for idx in range(1,num_arrays+1):
        if d[idx][0]==d[idx][1]: # homozygous
            if d[idx][0]==a[0][0]: # homozygous A
                snpinfo['calls'][idx-1] = (1, set(a[1]), a[1])
                hasmatch=True
            elif d[idx][0]==b[0][0]: # homozygous B
                snpinfo['calls'][idx-1] = (1, set(b[1]), b[1])
                hasmatch=True
            elif d[idx][0]=='-': # no call made
                snpinfo['calls'][idx-1] = (0,)
            else:
                print 'Homozygous, but not AA/BB : %s -- %s, %s' %(d[idx], a,b)
                # wtf?
                snpinfo['calls'][idx-1] = (0, set(),set(), frozenset(),frozenset())

        # verify it is a het call
        elif (d[idx][0]==a[0][0] and d[idx][1]==b[0][0]) or (d[idx][1]==a[0][0] and d[idx][0]==b[0][0]):
            snpinfo['calls'][idx-1] = (2, set(a[1]), set(b[1]), a[1], b[1])
            hasmatch=True

        else:
            print 'Unknown call : %s -- %s, %s' %(d[idx], a,b)
            snpinfo['calls'][idx-1] = (0, set(), set(), frozenset(),frozenset())
            continue

    if hasmatch:
        if chr not in chrsnps:
            chrsnps[chr] = set()
        chrsnps[chr].add( (pos, d[0]) )
        allsnps[snpinfo['snp']] = snpinfo['calls']
f.close()
print "%d SNPs loaded." % (len(allsnps),)

def _process_chr(thesnps, idx):
    for snpset in [sorted(thesnps), sorted(thesnps, reverse=True)]:
        last = None
        for pos,snp in snpset:
            if allsnps[snp][idx][0]!=1:
                continue

            if last is not None:
                if not allsnps[last[1]][idx][1].isdisjoint( allsnps[snp][idx][1] ):
                    allsnps[last[1]][idx][1].intersection_update( allsnps[snp][idx][1] )
                    allsnps[snp][idx][1].intersection_update( allsnps[last[1]][idx][1] )

            last = pos,snp

def _process_chr2(thesnps, idx):
    best_founders = {}

    for snpset in [sorted(thesnps), sorted(thesnps, reverse=True)]:
        previous_founders = deque()
        recent_counts = {}

        for pos,snp in snpset:
            if allsnps[snp][idx][0]!=1:
                continue

            if len(allsnps[snp][idx][1])==1:
                for fdr in allsnps[snp][idx][1]:
                    if fdr not in recent_counts:
                        recent_counts[ fdr ] = 1
                    else:
                        recent_counts[ fdr ] += 1
                    previous_founders.append( fdr )

                if len(previous_founders)==previous_max:
                    recent_counts[ previous_founders.popleft() ] -= 1

            else:
                best_num = 0
                best = None
                for fdr in allsnps[snp][idx][1]:
                    if fdr in recent_counts and recent_counts[fdr]>best_num:
                        best=fdr
                        best_num=recent_counts[fdr]

                if best is not None:
                    if snp not in best_founders:
                        best_founders[snp] = (best, best_num)
                    elif best_founders[snp][1]<best_num:
                        best_founders[snp] = (best, best_num)

    for snp in best_founders:
        allsnps[snp][idx][1].intersection_update( [best_founders[snp][0]] )


fout = open(sys.argv[2], 'w')
fout2 = open(sys.argv[3], 'w')
print >>fout, "SNP ID\tchr\tposition\t%s" % ("\t".join(header[1:num_arrays+1]), )
print >>fout2, "SNP ID\tchr\tposition\t%s" % ("\t".join(header[1:num_arrays+1]), )

for chr,snpset in chrsnps.items():
    chrscore = [0.0,0]
    print "  Chr"+chr
    for idx in range(0,num_arrays):
        _process_chr(snpset,idx)
        _process_chr2(snpset,idx)

    # output our founder calls and accuracy metric
    lastfoundercall=None
    scores = [deque() for x in range(0,num_arrays)]
    for pos,snp in sorted(snpset):
        chrfounders=[]
        goods=0
        for idx in range(0,num_arrays):
            if allsnps[snp][idx][0]==1 and len(allsnps[snp][idx][1])==1:
                f = ''.join(allsnps[snp][idx][1])
                chrfounders.append(' // '.join([f,f]))
                goods+=1

            #elif allsnps[snp][idx][0]==1 or allsnps[snp][idx][0]==2:
            #    chrfounders.append( ' // '.join(allsnps[snp][idx][1]) +' /// '+ ' // '.join(allsnps[snp][idx][2]) )
            #    goods+=1
            #else:
            #    chrfounders.append( '--' )

        if goods==num_arrays:
            lastfoundercall = chrfounders
            print >>fout, "%s\t%s\t%s\t%s" % (snp,chr,pos, '\t'.join(chrfounders))

        # calculate accuracy of calls so far...
        if lastfoundercall is None:
            # beginning of chromosome...
            continue

        curscores = []
        for idx in range(0,num_arrays):

            # do not score no-calls
            if allsnps[snp][idx][0]>0:
                score=0.0
                tc=[-1,-1]

                if allsnps[snp][idx][0]==2:
                    tc[0]=-2
                fdrs = lastfoundercall[idx].split(' // ')
                if fdrs[0] in allsnps[snp][idx][ tc[0] ]:
                    score+=0.5
                if fdrs[1] in allsnps[snp][idx][ tc[1] ]:
                    score+=0.5

                scores[idx].append(score)
                if len(scores[idx])>windowsize:
                    scores[idx].popleft()
                chrscore[0]+=score
                chrscore[1]+=1

            sum=0.0
            for score in scores[idx]:
                sum+=score
            sum = float(sum) / float(len(scores[idx]))
            curscores.append("%4.3f // %4.3f" % (sum,sum))

        print >>fout2, "%s\t%s\t%s\t%s" % (snp,chr,pos, '\t'.join(curscores))

    try:
        print "    Chr %s Accuracy: %4.3f" % (chr, chrscore[0]/float(chrscore[1]),)
    except:
        print "    Chr %s Accuracy: 0.000" % (chr,)

fout.close()
fout2.close()

