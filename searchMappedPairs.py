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

import subprocess
import os
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

            elif len(tmp_list) >= n:
                tmp_list.pop()
                #print tmp_list
                master_tmp_list.append(list(tmp_list))
                del tmp_list[:]
                total_span = 0
            else:
                del tmp_list[:]
                total_span = 0

        self.links = master_tmp_list

    def findRangeOfLinks(self, links):
        tmp_list = [x[0] for x in links]
        min1 = min(tmp_list)
        max1 = max(tmp_list)
        tmp_list = [x[1] for x in links]
        min2 = min(tmp_list)
        max2 = max(tmp_list)
        return (min1, max1, min2, max2)
    
    def generateGraphvizEdgeLabel(self, link):
        rangeLink = self.findRangeOfLinks(link)
        return '"(' +str(len(link))+') '+ str(rangeLink[0]) + '-' + str(rangeLink[1]) + ':' + str(rangeLink[2]) + '-' + str(rangeLink[3]) + '"'
    
    def to_graphviz(self, key):
        if len(self.links) > 0:
            for links in self.links:
                edge_label = self.generateGraphvizEdgeLabel(links)
                contig = key.split(':')
                parts1 = contig[0].split('_')
                parts2 = contig[1].split('_')
                label = 'N'+parts1[1]+'_L'+parts1[3]+' -> N'+parts2[1]+'_L'+parts2[3]+ '[len = 3 label = '+edge_label+'];\n'
                return label

    def to_cytoscapeConnections(self, key):
         if len(self.links) > 0:
            for links in self.links:
                contigs = key.split(':')
                return contigs[0]+'\t0\t'+contigs[1]+'\n'


            
class LinkFinder:
    def __init__(self, bamFile, clusterSize, mniNumLinks, endSize,
            outputPrefix):
        self.bamFile = bamFile
        self.clusterSize = clusterSize
        self.minNumLinks = mniNumLinks
        self.endLength = endSize
        self.hashOfLinks = {}
        self.outputPrefix = outputPrefix
    
    def linkList(self):
        return self.hashOfLinks

    # remove any contig pairings that have fewer links than the specified number    
    def filterLinksOnNumber(self):
        for k in self.hashOfLinks.keys():
            if (self.hashOfLinks[k].size() < i):
                del self.hashOfLinks[k]

    def addOrAppend(self,mate_contig, contig, read):
        if(mate_contig < contig):
            key = mate_contig +':'+contig
        else:
            key = contig+':'+mate_contig
        
        link = (read.mpos, read.pos)
                    
        if(key in self.hashOfLinks):
            self.hashOfLinks[key].append(link)
        else:
            s = ContigLinks( link)
            self.hashOfLinks[key] = s

    def findAllPossibleLinks(self):
        references = self.bamFile.references
        for contig in references:
            for read in self.bamFile.fetch(contig):
                if self.hasMissingMates(read, contig):
                    self.addOrAppend(self.bamFile.getrname(read.mrnm), contig, read )

    def findEndLinks(self):

        for contig,length in zip(self.bamFile.references, self.bamFile.lengths):

            # get the links at the start of the reference
            for read in self.bamFile.fetch(contig, 1, self.endLength):
                self.checkLink(read, length, contig)

            # now get oll of the links at the end of the reference
            for read in self.bamFile.fetch(contig, length - self.endLength, length):
                self.checkLink(read, length, contig)

    def checkLink(self, read, length, contig):
        if self.isMated(read):
            if self.hasMissingMates(read, contig):
                # mate is on a different contig
                self.addOrAppend(self.bamFile.getrname(read.mrnm), contig, read )
            elif self.isCircularLink(read.pos, read.mpos, length):
                # mates are at either end of the same contig
                self.addOrAppend(self.bamFile.getrname(read.mrnm), contig, read )

    # checks for a number of features for each aligned read.  If a read's mate is
    # on a different contig then it returns that contig name.  For all other
    # possibilities returns None
    def hasMissingMates(self, alignedRead, contig):
        mate_contig = self.bamFile.getrname(alignedRead.mrnm)
        if (mate_contig != contig):
            return True
        return False

    # checks the position of the read and it's mate to see if they are on oposite
    # ends of a contig. Returns True if they are, False if not
    def isCircularLink(self, alignedReadPos, matePos, contigLength):
        if ((alignedReadPos < self.endLength) and (matePos > contigLength - self.endLength)):
            return True
        elif (alignedReadPos > (contigLength - self.endLength) and (matePos < self.endLength)):
            return True
        else:
            return False

    def isMated(self, alignedRead):
        if alignedRead.is_paired:
            if not alignedRead.mate_is_unmapped:
                return True
        return False
    
    def clusterLinks(self):
        for k,v in self.hashOfLinks.iteritems():
           v.unique()
           v.findClusteredLinkPositions(self.clusterSize, self.minNumLinks)

    def printGraph(self):
        outFile = open(self.outputPrefix+'gv', 'w')
        outFile.write('digraph A {\n')
        for k,v in self.hashOfLinks.iteritems():
           buf =v.to_graphviz(k)
           if buf is not None:
               outFile.write(buf)
        outFile.write('}\n')

    def printCytoscape(self):
        connectFile = open(self.outputPrefix+"tab", 'w')
        for k,v in self.hashOfLinks.iteritems():
            # print two node per contig for beginning and end
            buf =v.to_cytoscapeConnections(k)
            if buf is not None:
                connectFile.write(buf)

def getSubList(subFile):
    f = open(subFile, 'r')
    subList = []
    for line in f:
        subList.append(line.rstrip())
    return subList

if __name__ =='__main__':
    
    # intialise the options parser
    parser = OptionParser("\n\n %prog -i <file.bam> ")
    parser.add_option("-i","--input",type="string",dest="inFile",help="the name of the input bam file")
    parser.add_option("-c","--clusterWidth", type="int", dest="clusterWidth", help="the maximum number of positions that two links can be apart to be considered part of the same cluster [default: 300]")
    parser.add_option("-n","--numberLinks", type="int", dest="numOfLinks", help="the number of links that two contigs must share for the links to even be considered 'real' [default: 3]")
    parser.add_option("-e","--end-length", type="int", dest="endOnly", help="Only consider links in the first or last x num of bases of a contig")
    parser.add_option("-t", "--cytoscape", action="store_true", dest="cytoscape", help="Output connections in cytoscape formatted tables. Default: false")
    parser.add_option("-p","--prefix", type="string", dest="prefix", default="connections.", help="Prefix for output files default: connections.")
    parser.add_option("-s","--subset", type="string",dest="subFile",help="file \
            containing a list of contig headers.  Only those contigs listed \
            will be evaluated for connections")
    # get and check options
    (opts, args) = parser.parse_args()
    if(opts.inFile is None):
        print"please specify the name of the input bam file with -i"
        parser.print_help()
        sys.exit(1)
    else:
        try:
            bamFile = pysam.Samfile(opts.inFile, 'rb')
        except:
            print"The input file must be in bam format"

    cluster_size = 300
    min_links = 3
    doSubset = False
    doEndOnly = False

    if opts.subFile is not None:
        subList = getSubList(opts.subFile)
        doSubset = True

    if (opts.clusterWidth is not None):
        cluster_size = opts.clusterWidth

    if (opts.numOfLinks is not None):
        min_links = opts.numOfLinks

    if opts.endOnly is None:
        finder = LinkFinder(bamFile, cluster_size, min_links, 0, opts.prefix)
        finder.findAllPossibleLinks()
    else:
        doEndOnly = True
        finder = LinkFinder(bamFile, cluster_size, min_links, opts.endOnly, opts.prefix)
        finder.findEndLinks()

    finder.clusterLinks()

    if opts.cytoscape is True:
        finder.printCytoscape()
    else:
        finder.printGraph()

