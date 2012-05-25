#!/usr/bin/python

#
#  Take a BAM file and calculate the coverage of each base
#  Output some general statistics for each column

from optparse import OptionParser
import sys
import pysam
import csv
import os

if __name__ == '__main__':
    #parse user options
    parser = OptionParser()
    parser.add_option("-a", "--averages", dest="averages", 
            default=False, action="store_true")
    # get and check options
    (opts,args) = parser.parse_args()
    # open the SAM/BAM
    for bam_file in args:
        try:
            samFile = pysam.Samfile(bam_file, 'rb')
        except:
            print "Unable to open BAM file -- did you supply a SAM file instead?"
            sys.exit(1)
        basename = os.path.splitext(bam_file)[0]

        output_csv = csv.writer(open(basename + '.cov.csv', "wb"),
                quoting = csv.QUOTE_NONNUMERIC )
        for reference, length in zip(samFile.references, samFile.lengths):
            #num_reads = samFile.count(reference, 0, length)
            #print num_reads / length
            tid = samFile.gettid(reference)
            if(tid != -1):
                tmp_cov = 0
                cov_vec = []
                for base in samFile.pileup(samFile.header['SQ'][tid]['SN'],
                                        0, length):
                    if opts.averages:
                        tmp_cov += base.n
                    else:
                        cov_vec.append(base.n)
                if opts.averages:
                    ave_coverage = float(tmp_cov) / float(length)
                    print reference, ave_coverage
                else:
                    if len(cov_vec) == 0:
                        cov_vec = [0 for x in range(length)]
                    cov_vec.insert(0,reference)
                    output_csv.writerow(cov_vec)
