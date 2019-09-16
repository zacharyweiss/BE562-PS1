#!/usr/bin/env python

####
# BE562: Problem Set 1 - string hashing/dotplots
#
# INSTRUCTIONS FOR USE:
# call program as follows:
#   ./ps1-dotplot.py <FASTA 1> <FASTA 2> <PLOTFILE>
#   e.g. ./ps1-dotplot.py human-hoxa-region.fa mouse-hoxa-region.fa dotplot.jpg
#
# Make sure the ps1-dotplot.py is marked as executable:
#   chmod +x ps1-dotplot.py
# or in windows with:
#   python ps1-dotplot.py human-hoxa-region.fa mouse-hoxa-region.fa dotplot.jpg
# once you have put python in your path
#
# matplotlib
# matplotlib is used to generate plots for this problem. It is a widely used
# python package. Installation instructions can be found at
# http://matplotlib.org/users/installing.html


import sys
from collections import defaultdict


def readSeq(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seq = []

    for line in stream:
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())

    return "".join(seq)


def quality(hits, slope1, slope2, offset1, offset2):
    """determines the quality of a list of hits"""

    goodhits = []

    for hit in hits:
        upper = slope1 * hit[0] + offset1
        lower = slope2 * hit[0] + offset2

        if lower < hit[1] < upper:
            goodhits.append(hit)

    return goodhits


def makeDotplot(filename, hits):
    """generate a dotplot from a list of hits
       filename may end in the following file extensions:
         *.ps, *.png, *.jpg
    """
    import matplotlib.pyplot as plt
    x, y = zip(*hits)

    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000)
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000

    hits2 = quality(hits, slope1, slope2, offset1, offset2)
    print("%.5f%% hits on diagonal" % (100 * len(hits2) / float(len(hits))))

    x2 = [0, 1e6]
    y2 = [offset1, slope1 * 1e6 + offset1]
    y3 = [offset2, slope2 * 1e6 + offset2]

    # create plot
    plt.plot(x2, y2, 'b', x2, y3, 'g')
    plt.scatter(x, y, s=1, c='r', edgecolors='none', marker=',')
    plt.title("dotplot (%d hits, %.5f%% hits on diagonal)" %
              (len(hits), 100 * len(hits2) / float(len(hits))))
    plt.xlabel("sequence 1")
    plt.ylabel("sequence 2")
    plt.xlim(x2)
    plt.ylim(x2)
    plt.tight_layout()

    # output plot
    plt.savefig(filename)


def main():

    # parse command-line arguments
    if len(sys.argv) < 4:
        print("you must call program as:  ")
        print("   ./ps1-2.py <FASTA 1> <FASTA 2> <PLOT FILE>")
        print("   PLOT FILE may be *.ps, *.png, *.jpg")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    plotfile = sys.argv[3]

    # read sequences
    print("reading sequences")
    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    #
    # You will need to modify the code below to

    # length of hash key
    kmerlen = 30

    # hash table for finding hits
    lookup = defaultdict(list)

    # store sequence hashes in hash table
    print("hashing seq1...")
    for i in range(len(seq1) - kmerlen + 1):
        key = seq1[i:i+kmerlen]
        lookup[key].append(i)

    # look up hashes in hash table
    print("hashing seq2...")
    hits = []
    for i in range(len(seq2) - kmerlen + 1):
        key = seq2[i:i+kmerlen]

        # store hits to hits list
        for hit in lookup.get(key, []):
            hits.append((i, hit))

    # hits should be a list of tuples
    # [(index1_in_seq2, index1_in_seq1),
    #  (index2_in_seq2, index2_in_seq1),
    #  ...]

    print("%d hits found" % len(hits))
    print("making plot...")

    makeDotplot(plotfile, hits)

if __name__ == "__main__":
    main()
