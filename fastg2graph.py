#!/usr/bin/env python
###############################################################################
#
# fastg2graph 
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "uqcskenn"
__copyright__ = "Copyright 2014"
__credits__ = ["uqcskenn"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "uqcskenn"
__email__ = "uqcskenn@rudd"
__status__ = "Development"

###############################################################################

import argparse
import sys
import networkx as nx

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')     

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def doWork( args ):
    G = nx.DiGraph()
    for name, seq, qual in readfq(args.fastg):
        name = name.rstrip()
        fields = name.split(':')
        if len(fields) == 1:
            continue
        
        name_fields = fields[0].split('_')
        # this is our 'key node' the other fields are the edges
        if fields[0][-1] == "'":
            fields[0] = fields[0][:-1]

        G.add_node(fields[0], length=int(name_fields[3]),
                coverage=float(name_fields[5]))
        for i in range(1,len(fields)):
            name_fields = fields[i].split('_')
            if fields[i][-1] == ';':
                fields[i] = fields[i][:-1]

            if fields[i][-1] == "'":
                # this is reverse complement therefore this node comes before
                # our current key node
                G.add_node(fields[i][:-1], length=int(name_fields[3]),
                        coverage=float(name_fields[5]))
                G.add_edge(fields[i][:-1], fields[0])
            else:
                G.add_node(fields[i], length=int(name_fields[3]),
                        coverage=float(name_fields[5]))
                G.add_edge(fields[0], fields[i])

    nx.write_gml(G, args.graphfile)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    
    parser = argparse.ArgumentParser()
    parser.add_argument('fastg', type=argparse.FileType(),
            help="Fastg input file")
    parser.add_argument('graphfile', help="output graph file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")
    
    # parse the arguments
    args = parser.parse_args()        

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

