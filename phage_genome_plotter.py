#!/usr/bin/env python
###############################################################################
#
# phage_genome_plotter.py - Tell me all about a phage and the spacers and snps
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

__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################

import argparse
import sys
from collections import defaultdict
import os
#import errno

import numpy as np
np.seterr(all='raise')

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.ticker import NullFormatter, MaxNLocator
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure

import vcf
from datetime import date

###############################################################################
###############################################################################
###############################################################################
###############################################################################

sample_names_to_dates = {
        '2705': date(2011, 5, 27),
        '0108': date(2011, 8, 1),
        '0809': date(2011, 9, 8),
        '2511': date(2011, 11, 25),
        '1801': date(2012, 1, 18),
        '2206': date(2011, 6, 22),
        }
viral_samples = {
        'viral_mapping_set_V90102': date(2012, 2, 1),
        'viral_mapping_set_V90104': date(2011, 4, 1),
        'viral_mapping_set_V90106': date(2011, 6, 1),
        'viral_mapping_set_V90308': date(2011, 8, 3),
        'viral_mapping_set_V90401': date(2012, 1, 4),
        'viral_mapping_set_V90709': date(2011, 9, 7),
        'viral_mapping_set_V90806': date(2011, 6, 8),
        'viral_mapping_set_V90903': date(2011, 3, 9),
        'viral_mapping_set_V91105': date(2011, 5, 11),
        'viral_mapping_set_V91210': date(2011, 10, 12),
        'viral_mapping_set_V91307': date(2011, 7, 13),
        'viral_mapping_set_V91409': date(2011, 9, 14),
        'viral_mapping_set_V91412': date(2011, 12, 14),
        'viral_mapping_set_V91506': date(2011, 6, 15),
        'viral_mapping_set_V91801': date(2012, 1, 18),
        'viral_mapping_set_V91802': date(2011, 2, 18),
        'viral_mapping_set_V91805': date(2011, 5, 18),
        'viral_mapping_set_V92007': date(2011, 7, 20),
        'viral_mapping_set_V92206': date(2011, 6, 22),
        'viral_mapping_set_V92704': date(2011, 4, 27),
        'viral_mapping_set_V92809': date(2011, 9, 28),
        'viral_mapping_set_V92810': date(2011, 10, 28),
        'viral_mapping_set_V93011': date(2011, 11, 3),
        'viral_mapping_set_V93108': date(2011, 8, 31),
        }


class Snp(object):
    def __init__(self, pos, alt, timepoints=None):
        super(Snp, self).__init__()
        self.pos = pos
        self.alt = alt
        if timepoints is None:
            self.timepoints = []
        else:
            self.timepoints = timepoints


class Spacer(object):
    def __init__(self, name, start, end, contig, timepoint=None):
        self.start = int(start)
        self.end = int(end)
        self.contig = contig
        self.name = name
        self.positions = set(range(self.start, self.end))
        self.timepoint = timepoint


class Contig(object):
    def __init__(self, length, timepoints=24):
        super(Contig, self).__init__()
        self.length = length
        # list of 2-tuple of the positions of and timepoints of snps and
        # spacers
        self.snps = []
        self.spacers = []

def fx_include(name, seq, qual, headers=None):
    """Returns sequences that contain the headers
    """
    return name in headers

def fx_parse(fp, callback=None, **kwargs):
        last = None  # this is a buffer keeping the last unprocessed line
        while True:  # mimic closure; is it a bad idea?
            if not last:  # the first record or a record following a fastq
                for l in fp:  # search for the start of the next record
                    if l[0] in '>@':  # fasta/q header line
                        last = l[:-1]  # save this line
                        break
            if not last:
                break
            name, seqs, last = last[1:].split()[0], [], None
            for l in fp:  # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+':  # this is a fasta record
                if callback:
                    if callback(name, ''.join(seqs), None, **kwargs):
                        yield name, ''.join(seqs), None  # yield a fasta record
                else:
                    yield name, ''.join(seqs), None  # yield a fasta record
                if not last:
                    break
            else:  # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp:  # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq):  # have read enough quality
                        last = None
                        if callback:
                            if callback(name, seq, ''.join(seqs), **kwargs):
                                yield name, seq, ''.join(seqs);  # yield a fastq record
                        else:
                            yield name, seq, ''.join(seqs);  # yield a fastq record
                        break
                if last:  # reach EOF before reading enough quality
                    if callback:
                        if callback(name, seq, None, **kwargs):
                            yield name, seq, None  # yield a fasta record instead
                    else:
                        yield name, seq, None  # yield a fasta record instead
                    break


def get_spacer_positions(fields):
    ''' take a line from the blast file and return a list of bp where the spacer
        targets
    '''
    # final field is the length of the spacer
    #sp_length = int(fields[-1])
    al_len = int(fields[3])
    mismatch = int(fields[4])
    gap = int(fields[5])
    start = int(fields[8])
    end = int(fields[9])

    #if float(al_len + mismatch + gap) < float(sp_length * 0.9):
    #    return None

    if start > end:
        start = end
        end = int(fields[8])

    return Spacer(fields[0], start, end, fields[1],
            sample_names_to_dates[fields[0][2:6]])

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    """ Main wrapper"""

    sorted_dates = sorted(viral_samples.values())
    # open a file of contig stats and get the contig lengths, rel abundance
    contig_stats = defaultdict(dict)
    with open(args.contigs) as fp:
        callback = None
        if args.filter is not None:
            callback = fx_include

        for name, seq, qual in fx_parse(fp, callback=callback,
                headers=args.filter):
            contig_stats[name]['length'] = len(seq)
            contig_stats[name]['spacers'] = []
            contig_stats[name]['snps'] = []
            contig_stats[name]['rel_abs'] = []

    # open the vcf file containing SNPs get the timepoints that the SNP is
    # present
    vcf_reader = vcf.Reader(open(args.snps))
    for record in vcf_reader:
        if record.QUAL is not None and int(record.QUAL) < int(args.snp_quality):
            continue
        if record.CHROM not in contig_stats:
            continue
        for sample in record.samples:
            if sample['GT'] != '0/0' and sample['GT'] != './.':
                try:
                    contig_stats[record.CHROM]['snps'].append([record.POS, viral_samples[os.path.splitext(sample.sample)[0]]])
                except:
                    contig_stats[record.CHROM]['snps'] = []
                    contig_stats[record.CHROM]['snps'].append([record.POS, viral_samples[os.path.splitext(sample.sample)[0]]])

    # open the Spacer file get the positions and timepoints
    with open(args.spacers) as sp:
        for line in sp:
            line = line.rstrip()
            fields = line.split('\t')
            if fields[0].startswith('M8') or fields[1] not in contig_stats:
                continue
            ret = get_spacer_positions(fields)
            try:
                contig_stats[fields[1]]['spacers'].append([ret.start,
                    ret.timepoint])
            except:
                contig_stats[fields[1]]['spacers'] = []
                contig_stats[fields[1]]['spacers'].append([ret.start,
                    ret.timepoint])
                
    # open relative abundance file and add that data in
    with open(args.rel_abs) as ra:
        rel_abs_header = True
        header = []
        times = []
        for line in ra:
            line = line.rstrip()
            fields = line.split('\t')
            if rel_abs_header:
                times = [viral_samples[x] for x in fields[1:]]
                rel_abs_header = False
            elif fields[0] in contig_stats:
                contig_stats[fields[0]]['rel_abs'] = zip(fields[1:], times)

    #-----
    # make a 2d plot
    nullfmt = NullFormatter()
    for contig_name, data in contig_stats.items():

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        left_h = left+width+0.02

        rect_scatter = [left, bottom, width, height]
        rect_histy = [left_h, bottom, 0.2, height]
        
        fig = plt.figure()
        axScatter = fig.add_axes(rect_scatter)
        axHisty = fig.add_axes(rect_histy)
        #axScatter = plt.axes(rect_scatter)
        #axHisty = plt.axes(rect_histy)
        axHisty.set_xlabel('relative abundance')
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        axHisty.set_ylim((mpl.dates.date2num(sorted_dates[0]),
            mpl.dates.date2num(sorted_dates[-1])))
        axHisty.xaxis.set_major_locator(MaxNLocator(4))

        #phi_rel_abundance = np.random.randn(len(sorted_dates))
        #ra_point = [x[0] for x in data['rel_abs']]
        ra_time = [x[1] for x in data['rel_abs']]
        ra_time = mpl.dates.date2num(ra_time)  # sorted(ra_time)
        sorted_indexes = np.argsort(ra_time)
        ra_point = []
        for i in sorted_indexes:
            ra_point.append(data['rel_abs'][i][0])
        axHisty.plot_date(ra_point, sorted(ra_time), xdate=False,
                ydate=True, color='0.7', marker=' ', linestyle='--')

        #ax.xaxis.set_ticks_position('top')
        axScatter.set_xlim([0, data['length']])
        axScatter.set_ylim((mpl.dates.date2num(sorted_dates[0]),
            mpl.dates.date2num(sorted_dates[-1])))
        axScatter.set_xlabel('genome position (bp)')
        axScatter.set_title(contig_name)

        sp_dates = [x[1] for x in data['spacers']]
        sp_dates = mpl.dates.date2num(sp_dates)
        sp_points = [x[0] for x in data['spacers']]
        axScatter.plot_date(sp_points, sp_dates, color='r', alpha=0.5,
                xdate=False, ydate=True)

        sn_dates =[x[1] for x in data['snps']]
        sn_dates = mpl.dates.date2num(sn_dates)
        sn_points = [x[0] for x in data['snps']]
        axScatter.plot_date(sn_points, sn_dates, color='0.7', marker='.',
                xdate=False, ydate=True)
         
        axScatter.tick_params(axis='y', labelsize='small')
        axHisty.tick_params(axis='y', labelsize='small')       
        #
        # Change the formatting of the xlabels to make them pretty
        #
        labels = axScatter.get_xticklabels() 
        for label in labels: 
            label.set_rotation(30)
            label.set_horizontalalignment('right')
            label.set_size('small')
        
        labels = axHisty.get_xticklabels()
        for label in labels:
            label.set_rotation(30)
            label.set_horizontalalignment('right')
            label.set_size('small')
            
        plt.savefig(os.path.join(args.output, contig_name + ".png"), dpi=300,
            format='png')
        #-----
        # clean up!
        plt.close(fig)
        #plt.close(axHisty)
        del fig

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--contigs', dest='contigs',
            help="fasta file of contigs")
    parser.add_argument('-r', '--relative-abundance', dest='rel_abs',
            help="text file containing relative abundance of phage at each \
                    timepoint")
    parser.add_argument('-q', '--snp-quality', dest='snp_quality', default=20,
            help="Minimum quality for a SNP")
    parser.add_argument('-s', '--snp', dest='snps',
            help="vcf file containing SNPs")
    parser.add_argument('-S', '--spacers', dest='spacers',
            help="tabular blast file containing spacer hits")
    parser.add_argument('-f', '--filter', dest='filter', action='append',
            help="only make plots for the named contig")
    parser.add_argument('-o', '--output', dest='output', default='.',
            help="output directory for the image")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
