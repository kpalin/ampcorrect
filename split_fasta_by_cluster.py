import sys
from Bio import SeqIO
from collections import defaultdict
import re
from random import shuffle

MAXsequences = 100

pat = re.compile("cluster=([^; ]+);")

outLines = defaultdict(list)

for inFasta in SeqIO.parse(sys.argv[1],"fasta"):
    cName = pat.search(inFasta.description).group(1)
    outLines[cName].append(inFasta)


for cName,seqs in outLines.iteritems():
    shuffle(seqs)
    for seqIdx in range(0,len(seqs),MAXsequences):
        fname = "%s.%d.fasta"%(cName,seqIdx)
        o = open(fname,"w")
        seqSlice = seqs[seqIdx:seqIdx+MAXsequences]
        SeqIO.write(seqSlice, o, "fasta")
        o.close()
        print fname, len(seqSlice)
