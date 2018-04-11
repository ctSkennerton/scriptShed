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
#    Copyright (C) 2011, 2012 Connor Skennerton
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

from optparse import OptionParser
import sys
import pysam
from sets import Set
from operator import itemgetter

# create an edge object between two contigs
# an edge needs to contain both of the contig names and the list of reads that have the edges
# the position of the reads in each contig and the number of reads
class ContigLinks:
    def __init__(self,link):
        self.links = []
        self.links.append(link)  
    
    # add on a two element tuple to the links list
    def append(self, link):
        self.links.append(link)
    
    # return the number of links between the two contigs    
    def size(self):
        return len(self.links)
    
    # convert the set of links into a list
    def to_list(self):
        return list(self.links)
    
    # convert back to a set
    def to_set(self, list):
        self.links = set(list)
    
    def unique(self):
        self.links = list(set(self.links))
        
    # sort by the first value of the tuple within the set
    # requires that the set be first converted into a list
    # since sets are unordered 
    def sortByFirst(self):
        self.links.sort(key = itemgetter(0))
        
    def sortBySecond(self):
        self.links.sort(key = itemgetter(1))
    
    def findClusteredLinkPositions(self, z, n):
        self.sortByFirst()
        master_tmp_list = []
        tmp_list = []
        total_span = 0
        #print "-------------------------------"
        for i in range(len(self.links) - 1):
            if total_span < z:
                tmp_list.append(self.links[i])
                total_span = self.links[i][0] - tmp_list[0][0]

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
def checkLink(bamFile, read, length, contig):
    if isMated(read):
        if hasMissingMates(bamFile, read, contig):
            # mate is on a different contig
            return True


def getSubList(subFile):
    f = open(subFile, 'r')
    subList = []
    for line in f:
        subList.append(line.rstrip())
    return subList
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
    args = parser.parse_args()

    endLength = 500
    sublist = None
    if args.wantedContigs is not None:
        sublist = getSubList(args.wantedContigs)

    try:
        bamFile = pysam.Samfile(args.bam, 'rb')
    except:
        print "The input file must be in bam format"
    if (opts.clusterWidth is not None):
        cluster_size = opts.clusterWidth

    if (opts.numOfLinks is not None):
        min_links = opts.numOfLinks

    finder = LinkFinder(bamFile, cluster_size, min_links, opts.endOnly,
            opts.prefix, doSubset, subList)
    finder.iterateThroughContigs()

    finder.clusterLinks()

    if opts.cytoscape is True:
        finder.printCytoscape()
    else:
        finder.printGraph()

