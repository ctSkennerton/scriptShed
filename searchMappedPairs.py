#!/usr/bin/env python
###############################################################################
#
#    searchMappedPairs.py
#
#    Given a bam file, this srcipt will calculate where the mate of the read
#    is mapped or unmappedand in which contig that mate is mapped to.  The aim
#    being to generate a graphviz image and report file of all of the pairs that
#    map to a different contigs and to no contigs at all.  Hopefully this will
#    help with assemblies and stuff like that to piece together contigs.
#
#    Copyright (C) 2011, 2012, 2014 Connor Skennerton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

import argparse
import sys
import pysam
import networkx as nx
from operator import itemgetter


def getSubList(subFile):
    subList = set()
    for line in subFile:
        subList.add(line.rstrip())
    return subList

def findEndLinks(G, bamFile, contig, length, endLength=500):
    # get the links at the start of the reference
    for read in bamFile.fetch(contig, 1, endLength):
        if checkLink(bamFile, read, length, contig):
            mate_contig = bamFile.getrname(read.mrnm)
            G.add_node(contig, length=length)
            G.add_node(mate_contig, length=bamFile.lengths[read.mrnm])
            addLink(G, contig, mate_contig)

    # now get oll of the links at the end of the reference
    for read in bamFile.fetch(contig, length - endLength, length):
        if checkLink(bamFile, read, length, contig):
            mate_contig = bamFile.getrname(read.mrnm)
            G.add_node(contig, length=length)
            G.add_node(mate_contig, length=bamFile.lengths[read.mrnm])
            addLink(G, contig, mate_contig)


def addLink(G, contig, mate_contig):

    if contig < mate_contig:
        G.add_edge(contig, mate_contig)
        try:
            G[contig][mate_contig]['weight'] += 1
        except:
            G[contig][mate_contig]['weight'] = 1
    else:
        G.add_edge(mate_contig, contig)
        try:
            G[mate_contig][contig]['weight'] += 1
        except:
            G[mate_contig][contig]['weight'] = 1


def checkLink(bamFile, read, length, contig):
    if isMated(read):
        if hasMissingMates(bamFile, read, contig):
            # mate is on a different contig
            return True

# checks for a number of features for each aligned read.  If a read's mate is
# on a different contig then it returns that contig name.  For all other
# possibilities returns None
def hasMissingMates(bamFile, alignedRead, contig):
    mate_contig = bamFile.getrname(alignedRead.mrnm)
    if (mate_contig != contig):
        return True
    return False

# checks the position of the read and it's mate to see if they are on oposite
# ends of a contig. Returns True if they are, False if not
def isCircularLink(alignedReadPos, matePos, contigLength, endLength=500):
    if ((alignedReadPos < endLength) and (matePos > contigLength - endLength)):
        return True
    elif (alignedReadPos > (contigLength - endLength) and (matePos < endLength)):
        return True
    else:
        return False

def isMated(alignedRead):
    if alignedRead.is_paired:
        if not alignedRead.mate_is_unmapped:
            return True
    return False


if __name__ =='__main__':

    # intialise the options parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("bam", help="the name of the input bam file")
    parser.add_argument('outfile', help='Name of output file of graph in GML format')
    parser.add_argument('-f', '--wanted-contigs', dest='wantedContigs', type=argparse.FileType(),
            help='A file of contig names to be considered')
    parser.add_argument("-n","--numberLinks", type=int, dest="numOfLinks", default=3,
            help="the number of links that two contigs must share for the links to even be considered 'real'")
    parser.add_argument('-m', '--min-contig-len', type=int, dest='minContigLen', default=500,
            help='The minimum length of the contig to be considered for adding links')

    # get and check options
    args = parser.parse_args()

    endLength = 500
    sublist = None
    if args.wantedContigs is not None:
        sublist = getSubList(args.wantedContigs)

    try:
        bamFile = pysam.Samfile(args.bam, 'rb')
    except:
        print "The input file must be in bam format"
        sys.exit(1)

    G = nx.Graph()
    for contig, length in zip(bamFile.references, bamFile.lengths):

        if length < args.minContigLen:
            continue

        if contig not in sublist:
            continue

        findEndLinks(G, bamFile, contig, length)

    # now subset the graph based on the edge conditions

    SG = nx.Graph( [ (u,v,d) for u,v,d in G.edges(data=True) if d ['weight'] >= args.numOfLinks] )

    nx.write_gml(SG, args.outfile)
