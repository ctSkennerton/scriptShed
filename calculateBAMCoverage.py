#!/usr/bin/python

#
#  Take a BAM file and calculate the coverage of each base
#  Output some general statistics for each column

import subprocess
import os
from optparse import OptionParser
import sys
import pysam
from Bio import SeqIO


if __name__ == '__main__':
    #parse user options
    parser = OptionParser("\n\n %prog -b bamFile ")
    parser.add_option("-b", "--bam_file", type="string",
                    dest="samFileName", help="Specify a name for the sam file")

    # get and check options
    (opts, args) = parser.parse_args()
    if(opts.samFileName is None):
        print("please specify the name of the sam file with -s")
        parser.print_help()
        sys.exit(1)
    # open the SAM/BAM
    try:
        samFile = pysam.Samfile(opts.samFileName, 'rb')
    except:
        print "Unable to open BAM file -- did you supply a SAM file instead?"
        sys.exit(1)
    for reference, length in zip(samFile.references, samFile.lengths):
        #num_reads = samFile.count(reference, 0, length)
        #print num_reads / length
        tid = samFile.gettid(reference)
        if(tid != -1):
            tmp_cov = 0
            for base in samFile.pileup(samFile.header['SQ'][tid]['SN'],
                                    0, length):
                tmp_cov += base.n
            ave_coverage = float(tmp_cov) / float(length)
            print reference, ave_coverage
