"""Split reads in a fasta file according to reference sequences they have been mappend in a SAM file"""
import sys
from Bio import SeqIO
from collections import Counter

MAXsequences = 100

SAMish, FASTA = sys.argv[1:]

partFile = open(SAMish)

groups = {}
for line in partFile:
    if line.startswith("@"):
        continue
    line = line.split("\t")
    if line[2] in groups:
        exclude.add(line[2])
        groups[line[2]] = None
    else:
        groups[line[2]] = line[0]

for r in exclude:
    print "Excluding ambiguous read", r
    del (groups[r])

outIO = {}

group_counter = Counter()

for grp in set(groups.values() + ["rest"]):
    fname = "%s.fasta" % (grp)
    fname = re.sub("[^a-zA-Z0-9.-]", "_", fname)
    o = open(fname, "w")
    outIO[grp] = o

for inFasta in SeqIO.parse(FASTA, "fasta"):
    g = groups.get(inFasta.id, "rest")
    o = outIO[g]
    SeqIO.write(inFasta, o, "fasta")
    group_counter[g] += 1

for o in outIO.values():
    o.close()

print group_counter
