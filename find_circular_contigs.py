#!/usr/bin/env python
###############################################################################
#
# find_circular_contigs 
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
from __future__ import division
__author__ = "uqcskenn"
__copyright__ = "Copyright 2013"
__credits__ = ["uqcskenn"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "uqcskenn"
__email__ = "uqcskenn@rudd"
__status__ = "Development"

###############################################################################

import argparse
import sys
import pysam

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
class BamCallback:
    def __init__(self, refLength, start=700, end=700):
        self.circularPairs = {}
        self.otherContigLinks = {}
        self.refLen = refLength
        self.startBoundary = start
        self.endBoundary = refLength - end
        self.Verbose= False

    def __call__(self, alignedRead):
        if alignedRead.is_paired and not alignedRead.is_unmapped and not alignedRead.mate_is_unmapped:
            if alignedRead.tid != alignedRead.rnext:
                if alignedRead.pos <= self.startBoundary:
                    try:
                        self.otherContigLinks[alignedRead.qname].append(alignedRead)
                    except KeyError:
                        self.otherContigLinks[alignedRead.qname] = []
                        self.otherContigLinks[alignedRead.qname].append(alignedRead)
                elif alignedRead.pos >= self.endBoundary:
                    try:
                        self.otherContigLinks[alignedRead.qname].append(alignedRead)
                    except KeyError:
                        self.otherContigLinks[alignedRead.qname] = []
                        self.otherContigLinks[alignedRead.qname].append(alignedRead)
            else:
                if alignedRead.pos <= self.startBoundary and alignedRead.pnext >= self.endBoundary:
                    try:
                        self.circularPairs[alignedRead.qname].append(alignedRead)
                    except KeyError:
                        self.circularPairs[alignedRead.qname] = []
                        self.circularPairs[alignedRead.qname].append(alignedRead)
                elif alignedRead.pos >= self.endBoundary and alignedRead.pnext <= self.startBoundary:
                    try:
                        self.circularPairs[alignedRead.qname].append(alignedRead)
                    except KeyError:
                        self.circularPairs[alignedRead.qname] = []
                        self.circularPairs[alignedRead.qname].append(alignedRead)

def add_circ_to_db(args, circ_contigs):
    import sqlite3
    conn=sqlite3.connect(args.database, detect_types=sqlite3.PARSE_DECLTYPES)
    cur=conn.cursor()
    for name in circ_contigs:
        cur.execute('''UPDATE contigs SET Circular=1 WHERE Name = ?''', (name,))
    conn.commit()

def filter_contigs(infile):
    all_set = set()
    bad_set = set()
    for line in infile:
        fields = line.rstrip().split()
        if len(fields) != 14:
            raise RuntimeError("provide tabular blast+ with -outfmt '6 std qlen slen'")
        if fields[0] == fields[1]:
            continue
        all_set.add(fields[0])
        if int(fields[3]) / int(fields[12]) >= 0.9:
            if int(fields[13]) / int(fields[12]) >= 1.5:
                #print '\t'.join(fields)
                bad_set.add(fields[0])
    return all_set - bad_set

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    filtered_contigs = set()
    if args.blast is not None:
        filtered_contigs = filter_contigs(args.blast)

    circ_contigs = set()
    for bamfile in args.bamfile:
        bam = pysam.Samfile(bamfile, 'rb')
        total_count = 0.0
        circ_count = 0.0
        for reference, length in zip(bam.references, bam.lengths):
            if args.blast is not None and reference not in filtered_contigs:
                continue
            total_count += 1.0
            if length < args.min_frag_length:
                continue
            rl = BamCallback(length)
            bam.fetch(reference, 0, length, callback = rl )
            if len(rl.circularPairs) >= args.links:
                circ_contigs.add(reference)
                if not args.quiet:
                    print reference, length, len(rl.circularPairs), len(rl.otherContigLinks)
                circ_count += 1.0
        if args.summary:
            print "total: %d circular: %d percentage: %f" % (total_count, circ_count,circ_count / total_count)
    if args.database is not None:
        add_circ_to_db(args, circ_contigs)
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    
    parser = argparse.ArgumentParser()
    parser.add_argument('bamfile', nargs='+', help="Input bam file to check")
    parser.add_argument('-q', '--quiet', dest='quiet', action='store_true',
            help="do not print individual contig stats")
    parser.add_argument('-S', '--summary', dest='summary', action='store_true',
            help="print summary at the end")
    parser.add_argument('-d', '--database', dest='database',
            help="sqlite database name for updating information")
    parser.add_argument('-m', '--min-contig-length', type=int, default=3000,
            help="The minimum contig length to consider for circularity",
            dest='min_frag_length')
    parser.add_argument('-s', '--start-boundary', default=700, type=int,
            help="only consider reads that begin within this boundary on the"\
            " reference", dest='start')
    parser.add_argument('-l', '--minimum-links', default=3, type=int, help="The\
            minimum number of pairs between the end of the contig for it to be\
            considered circular", dest='links')
    parser.add_argument('-e', '--end-boundary', default=700, type=int,
            dest='end', help="only consider reads that begin within this boundary on the"\
            " reference")
    parser.add_argument('-b', '--blast-output', type=argparse.FileType('r'),
            help='A tabular blast output file containing blast hits of contigs \
            to each other.  Will filter contigs out if they are fragments of \
            other contigs', dest='blast')
    
    # parse the arguments
    args = parser.parse_args()        

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

