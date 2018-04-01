import sys
import collections
from Bio import SeqIO

#fasta file with seqs of interest
ff=sys.argv[1]
#file with blast results
blastf=sys.argv[2]

r1=open(blastf, "r")

seq_d={}
count_d={}
coord_d={}

for seq in SeqIO.parse(ff, "fasta"):
    id=seq.id
   # print(str(id))
    length=len(str(seq.seq))
    cnt=0
    ar=[]
    while(cnt < length):
        ar.append(0)
        cnt=cnt+1
    seq_d[id]=ar
    count_d[id]=[]
    coord_d[id]=set()

for l in r1:
    l_arr=l.rstrip().split("\t")
    rname=l_arr[0]
    name=l_arr[1]
    start=int(l_arr[8])
    end=int(l_arr[9])
    count_d[name].append(rname)
    coord_d[name].add(start)
    coord_d[name].add(end)
    if (start < end):
#takes into account indexing from 0
        c=start-1
        while(c < end):
            pos=seq_d[name][c]
            npos=pos+1
            seq_d[name][c]=npos
            c=c+1
    else:
        c=start-1
        while(c > end-2):
            pos=seq_d[name][c]
            npos=pos+1
            seq_d[name][c]=npos
            c=c-1

for g in seq_d:
    cov_p=0
    for p in seq_d[g]:
        if(p > 0):
            cov_p=cov_p+1
    frac=float(cov_p)/float(len(seq_d[g]))
    print(g+"\t"+str(frac)+"\t"+str(len(count_d[g]))+"\t"+str(len(coord_d[g])))
