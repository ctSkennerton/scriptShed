#!/usr/bin/env python
###############################################################################
#
# blast2sam.py - description!
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

from __future__ import print_function

__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__status__ = "Development"

###############################################################################

import argparse
import sys
import pysam
from Bio.Blast import NCBIXML as blastxml
from Bio.Seq import Seq
from Bio import SeqIO
#from Bio.Seq import reverse_complement
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
def makeCigar(hsp, align_length):
    ''' create a pysam CIGAR tuple, given an alignment from a blast file
    '''
    # cigar
    # The cigar alignment (None if not present)
    # The operation code comes first then the number of basses
    # The alignment is returned as a list of operations. The operations are:
    #  M	BAM_CMATCH	0
    #  I	BAM_CINS	1
    #  D	BAM_CDEL	2
    #  N	BAM_CREF_SKIP	3
    #  S	BAM_CSOFT_CLIP	4
    #  H	BAM_CHARD_CLIP	5
    #  P	BAM_CPAD	6
    #  =	BAM_CEQUAL	7
    #  X	BAM_CDIFF	8

    query = hsp.query
    match = hsp.match
    subject = hsp.sbjct
    qlen = len(query)
    
    current_operation = -1
    current_operation_count = 0
    operations = []
    if hsp.query_start != 1:
        operations.append((5, hsp.query_start - 1))
    for i in range(len(query)):
        if match[i] != '|':
            if query[i] != '-' and subject[i] != '-':  # straight mismatch
                if current_operation != 0:
                    if current_operation != -1:
                        operations.append((current_operation, current_operation_count))
                    current_operation = 0
                    current_operation_count = 1
                else:
                    current_operation_count += 1
            elif query[i] == '-' and subject[i] != '-':  # Deletion in query
                if current_operation != 2:
                    if current_operation != -1:
                        operations.append((current_operation, current_operation_count))
                    current_operation = 2
                    current_operation_count = 1
                else:
                    current_operation_count += 1
            elif query[i] != '-' and subject[i] == '-':  # Insertion in query
                if current_operation != 1:
                    if current_operation != -1:
                        operations.append((current_operation, current_operation_count))
                    current_operation = 1
                    current_operation_count = 1
                else:
                    current_operation_count += 1
            else:
                raise RuntimeError('Should not have got here')
        else:
            if current_operation != 7:
                if current_operation != -1:
                    operations.append((current_operation, current_operation_count))
                current_operation = 7
                current_operation_count = 1
            else:
                current_operation_count += 1
    operations.append((current_operation, current_operation_count))
    if hsp.query_end < align_length:
        #print(align_length, ' ', qlen, ' ', qlen + hsp.query_start - 1)
        operations.append((5, align_length - hsp.query_end))

    # sanity check
    # look at SAM documentation to see this calculation
    count = 0
    for op in operations:
        if op[0] != 2:
            count += op[1]
    if count != align_length:
        print("CIGAR does not match align length: %s\t%i\t%i\n%s" %
                (str(operations), count, align_length, str(hsp)))

    return tuple(operations)

def parseReferences(infile, informat='fasta'):
    headers = []
    header_lookup = {}
    counter = 0
    for record in SeqIO.parse(infile, informat):
        headers.append({'LN': len(record), 'SN': record.id})
        header_lookup[record.id] = counter
        counter += 1
    return headers, header_lookup

def doWork( args ):
    """ Main wrapper"""

    # make sam header
    header = {
                'HD': {'VN': '1.0'}
             }
    headers, header_lookup = parseReferences(args.ref)
    header['SQ'] = headers
    # open outfile
    outfile = pysam.Samfile( args.samfile, "wh", header = header )

    # parse in the blast file 
    blast_records = blastxml.parse(open(args.blast))
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                #print(alignment.title)
                read = pysam.AlignedRead()
                read.qname = blast_record.query
                read.flag = 0  # fix up reverse complementing
                dna = hsp.query.replace('-', '')
                #print(hsp.frame)
                cigar = makeCigar(hsp, blast_record.query_letters)  # represented as tuple of 2-tuples
                if hsp.frame[1] ^ hsp.frame[0]:
                    seq = Seq(dna)
                    rc = seq.reverse_complement()
                    read.seq = str(rc)
                    read.flag |= 0x10
                    read.pos = hsp.sbjct_end - 1
                    read.cigar = cigar[::-1]
                else:
                    read.seq = dna
                    read.pos = hsp.sbjct_start - 1
                    read.cigar = cigar
                read.rname = header_lookup[alignment.hit_def]  # index to list of headers
                read.mapq = 255  # phred scaled probability score
                read.mrnm = -1  # index of the mate
                read.mpos = -1  # position of the mate
                read.tlen = 0  # insert size of the mates
                outfile.write(read)
                #print(read)
    outfile.close()


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('blast', help="Blast XML file")
    parser.add_argument('ref', help='Fasta formatted file containing the' \
            ' referece sequences')
    parser.add_argument('samfile', help="Name of the output SAM file")
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
