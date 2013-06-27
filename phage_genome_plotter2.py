#!/usr/bin/env python
###############################################################################
#
# phage_genome_plotter2.py - Tell me all about a phage and the spacers and snps
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
import sqlite3
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
names_to_times = {    
    'M92206': date(2011, 6, 22),
    'M92705': date(2011, 5, 27),
    'M90108': date(2011, 8, 1),
    'M80509': date(2011, 9, 5),
    'M81612': date(2011, 12, 16),
    'M81706': date(2011, 6, 17),
    'M91801': date(2012, 1, 18),
    'M90809': date(2011, 9, 8),
    'M92511': date(2011, 11, 25),
    'V90102': date(2012, 2, 1),
    'V90104': date(2011, 4, 1),
    'V90106': date(2011, 6, 1),
    'V90308': date(2011, 8, 3),
    'V90401': date(2012, 1, 4),
    'V90709': date(2011, 9, 7),
    'V90806': date(2011, 6, 8),
    'V90903': date(2011, 3, 9),
    'V91210': date(2011, 10, 12),
    'V91307': date(2011, 7, 13),
    'V91409': date(2011, 9, 14),
    'V91412': date(2011, 12, 14),
    'V91506': date(2011, 6, 15),
    'V91801': date(2012, 1, 18),
    'V91802': date(2011, 2, 18),
    'V91805': date(2011, 5, 18),
    'V92007': date(2011, 7, 20),
    'V92206': date(2011, 6, 22),
    'V92704': date(2011, 4, 27),
    'V92809': date(2011, 9, 28),
    'V92810': date(2011, 10, 28),
    'V93108': date(2011, 8, 31)
    }
names_to_times['ACSBR9'] = names_to_times['V91506']

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
    def __init__(self, start, end, contig, crispr, host=None, timepoint=None):
        if start > end:
            self.start = end
            self.end = start
        else:
            self.start = start
            self.end = end
        self.contig = contig
        self.timepoint = timepoint
        self.host = host
        self.crispr = crispr


class Contig(object):
    def __init__(self, Id, Name, length):
        super(Contig, self).__init__()
        self.id = Id
        self.Name = Name
        self.length = length
        # list of 2-tuple of the positions of and timepoints of snps and
        # spacers

        self.snps = []
        self.spacers = []
        self.relabs = {}
        
    def __len__(self):
        return self.length
    
    def relative_abundance(self):
        time = sorted(self.relabs.keys())
        points = [self.relabs[x] for x in time]
        return time, points
        
    def get_spacers(self):
        ret = {} #defaultdict(dict)
        for sp in self.spacers:
            if sp.crispr in ret:
                ret[sp.crispr]['host'] = sp.host
                ret[sp.crispr]['data'].append([sp.start, sp.timepoint])
            else:
                ret[sp.crispr] = dict()
                ret[sp.crispr]['host'] = sp.host
                ret[sp.crispr]['data'] = []
                ret[sp.crispr]['data'].append([sp.start, sp.timepoint])
                
        return ret

class ScaffoldFragment(Contig):
    def __init__(self, Id, Name, start, Length, complement):
        super(ScaffoldFragment, self).__init__(Id, Name, Length)
        self.start = start
        self.complement = complement
        
class Scaffold(object):
    def __init__(self, Id, Name, Length):
        super(Scaffold, self).__init__()
        self.id = Id
        self.name = Name
        self.length = Length
        self.fragments = {}
        
    def __len__(self):
        return self.length
    
    def relative_abundance(self):
        sums = defaultdict(float)
        for fragid, fragobj in self.fragments.items():
            frag_relabs = fragobj.relative_abundance()
            for t, r in frag_relabs.items():
                sums[t] += r
        frag_count = len(self.fragments)
        for t in sums.keys():
            sums[t] = sums[t] / frag_count
        time = sorted(sums.keys())
        points = [sums[x] for x in time]
        return time, points
    
    def get_spacers(self):
        ret = defaultdict(dict)
        for frag in self.fragments.values():
            for sp in frag.spacers:
                ret[sp.crispr]['host'] = sp.host
                try:
                    ret[sp.crispr]['data'].append([sp.start + frag.start, sp.timepoint])
                except TypeError:
                    ret[sp.crispr]['data'] = [[sp.start + frag.start, sp.timepoint]]  
        return ret      
            
def get_contig_length(cur, name):
    cur.execute('''SELECT Length FROM contigs WHERE Id = ?''', (name,))
    result = cur.fetchone()
    if not result:
        raise RuntimeError('could not find contigs name %s in the database' % name)
    return result[0]

def get_spacers_for_contig(cur, name):
    cur.execute(''' create temp view spacer_host_match_times as
        select CrisprClusterID, HostID, SpacerID, times.Timepoint, ContigID, Start, End 
        from times join (
            select crispr_host_connector.CrisprClusterID, 
                crispr_host_connector.HostID, 
                crispr_host_connector.SpacerID, 
                crispr_host_connector.SpacerTp, 
                spacer_matches.ContigID,
                spacer_matches.Start, 
                spacer_matches.End 
                FROM crispr_host_connector JOIN spacer_matches 
                ON crispr_host_connector.SpacerID = spacer_matches.SpacerID
            ) on times.Id = SpacerTp ''')

    cur.execute('''select CrisprClusterID, HostID, Start, End, Timepoint from
            spacer_host_match_times where ContigID = ?''', (name,))

    result = cur.fetchall()
    if not result:
        return dict()
        #raise RuntimeError('problem with returning spacers for %s' % name)
    ret = []
    for row in result:
        ret.append(Spacer(row[2], row[3], name, row[0], host=row[1], timepoint=row[4]))
    return ret

def get_relative_abundance(cur, name):
    cur.execute('''SELECT times.Timepoint, AVG(Relab) 
        FROM times 
        JOIN (
            SELECT Timepoint AS Tp, Relab 
            FROM contig_relabs WHERE ContigID IN (
                SELECT Id FROM contigs WHERE Cluster = (
                    SELECT Cluster FROM contigs WHERE Id = ?)
                )
            ) On times.Id=Tp 
            GROUP BY times.Timepoint''', (name,))

    result = cur.fetchall()
    if not result:
        raise RuntimeError('cannot get the relative abundance of %s' % name)
    
    return dict(result)

def get_contigs_for_scaffold(cur, scaff_id):
    fragments = {}
    cur.execute('''select contigs.Name, contigs.Id, scaffold_fragments.Start, contigs.Length, scaffold_fragments.Oreintation 
                    from scaffold_fragments 
                    join contigs 
                    on ContigID = Id 
                    where ScaffoldID = ?''', (scaff_id,))
    result = cur.fetchall()
    if result is None:
        raise RuntimeError('Cannot find contigs for scaffold %d' % scaff_id)
    else:
        for row in result:
            complement = False
            if row[4] == '-':
                complement = True
            fragments[row[1]] = ScaffoldFragment(row[1], row[0], row[2], row[3], complement)
    return fragments
            
        

def generate_contig_ids(cur, name=None, name_file=None):
    id_map = {}
    query = []
    if name is None and name_file is None:
        raise RuntimeError('must provide at least one of name or name_file')
        
    if name is not None:
        if isinstance(name, list):
            query.extend(name)
        if isinstance(name, str):
            query.append(name)
    
    if name_file is not None:
        for line in name_file:
            line = line.rstrip()
            query.append(line)
    try:
        cur.execute('''SELECT Name, Id, Flag, Length FROM contigs WHERE Name IN (?)''', tuple(query))
    except sqlite3.ProgrammingError:
        print query
        return {}
    result = cur.fetchall()
    if not result:
        raise RuntimeError('wierd')
    for row in result:
        if row[2] & 4:
            id_map[row[0]] = Scaffold(row[1], row[0], row[3])
            id_map[row[0]].fragments = get_contigs_for_scaffold(cur, row[1])
        else:
            id_map[row[0]] = Contig(row[1], row[0], row[3])
             
    return id_map

def get_all_timepoints(cur):
    cur.execute('''SELECT Timepoint FROM times''')
    result = cur.fetchall()
    if not result:
        raise RuntimeError(" can't get times, database corrupt")
    return sorted([x[0] for x in result])           

def calc_host_relab(cur, host):
    cur.execute('''SELECT Relab, times.Timepoint
        FROM times JOIN (
            SELECT Timepoint AS Tp, Relab 
            FROM host_relabs WHERE HostID = ?
        ) ON Tp=times.Id WHERE times.Reactor='SBR9' AND times.flag & 4''', (host,))

    result = cur.fetchall()
    if not result:
        raise RuntimeError('cannot get the relative abundance of %s' % cluster)
    return result

def unzip(data):
    v, h = zip(*data)
    v = [x if x is not None else 0 for x in v ]
    return v, h
###############################################################################
###############################################################################
###############################################################################
###############################################################################

mpl_line_styles = ['-' , '--' , '-.' , ':']
mpl_line_marker_styles = [ 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd']
mpl_line_color_styles = ['b', 'g', 'r', 'c', 'm', 'y']
mpl_combinations = []
for k in mpl_line_styles:
    for i in mpl_line_marker_styles:
        for j in mpl_line_color_styles:
            mpl_combinations.append((i, j, k))

def doWork( args ):
    """ Main wrapper"""
    conn=sqlite3.connect(args.database, detect_types=sqlite3.PARSE_DECLTYPES)
    cur=conn.cursor()
    sorted_dates = get_all_timepoints(cur)
    id_map = generate_contig_ids(cur, args.contigs, args.contigs_file)
    for contigobj in id_map.values():
        if isinstance(contigobj, Scaffold):
            for fragid, fragobj in contigobj.fragments.items():
                fragobj.spacers.extend(get_spacers_for_contig(cur, fragid))
                fragobj.rel_abs = get_relative_abundance(cur, fragid)
        else:
            contigobj.spacers.extend(get_spacers_for_contig(cur, contigobj.id))
            contigobj.rel_abs = get_relative_abundance(cur, contigobj.id)

    # open the vcf file containing SNPs get the timepoints that the SNP is
    # present
    vcf_reader = vcf.Reader(open(args.snps))
    for record in vcf_reader:
        if record.QUAL is not None and int(record.QUAL) < int(args.snp_quality):
            continue
        if record.CHROM not in id_map:
            continue
        for sample in record.samples:
            if sample['GT'] != '0/0' and sample['GT'] != './.' and sample['GT'] != '0':
                id_map[record.CHROM].snps.append([record.POS,
                    names_to_times[sample.sample]])

            

    #-----
    # make a 2d plot
    nullfmt = NullFormatter()
    for contig_name, data in id_map.items():

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        left_h = left+width+0.02

        rect_scatter = [left, bottom, width, height]
        rect_histy = [left_h, bottom, 0.2, height]
        
        fig = plt.figure()
        axScatter = fig.add_axes(rect_scatter)
        axHisty = fig.add_axes(rect_histy)

        axHisty.set_xlabel('relative abundance')
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        axHisty.set_ylim((mpl.dates.date2num(sorted_dates[0]),
            mpl.dates.date2num(sorted_dates[-1])))
            
        axHisty.xaxis.set_major_locator(MaxNLocator(4))

        ra_time, ra_point = data.relative_abundance()
        #ra_time = [x[1] for x in data['rel_abs']]
        ra_time = mpl.dates.date2num(ra_time)

        axHisty.plot_date(ra_point, ra_time, xdate=False,
                ydate=True, color='0.7', marker=' ', linestyle='--')

        #ax.xaxis.set_ticks_position('top')
        axScatter.set_xlim([0, len(data)])
        axScatter.set_ylim((mpl.dates.date2num(sorted_dates[0]),
            mpl.dates.date2num(sorted_dates[-1])))
        axScatter.set_xlabel('genome position (bp)')
        axScatter.set_title(contig_name)

        # extra axis for plotting host data
        axHisty2 = axHisty.twiny()
        axHisty2.set_ylim((mpl.dates.date2num(sorted_dates[0]),
            mpl.dates.date2num(sorted_dates[-1])))
        
        axHisty2.xaxis.set_major_locator(MaxNLocator(4))
        axHisty2.yaxis.set_major_formatter(nullfmt)

        labels = axHisty2.get_xticklabels()
        for label in labels:
            label.set_rotation(-30)
            label.set_horizontalalignment('right')
            label.set_size('small')
        plotted_host = False
        marker_style_index = 0
        
        spacer_data = data.get_spacers()
        
        for crispr_cluster, cluster_data in spacer_data.items():
            if cluster_data['host'] is not None:
                plotted_host = True
                host_relab_data = calc_host_relab(cur, cluster_data['host'])
                v, h = unzip(host_relab_data)
                axHisty2.plot_date( v, list(h), xdate=False, ydate=True,
                        label=cluster_data['host'], alpha=0.5,
                        ls=mpl_combinations[marker_style_index][2],
                        marker=mpl_combinations[marker_style_index][0],
                        color=mpl_combinations[marker_style_index][1])

            sp_points, sp_dates = unzip(cluster_data['data'])
            sp_dates = mpl.dates.date2num(sp_dates)
            axScatter.plot_date(sp_points, sp_dates, label=crispr_cluster, alpha=0.5,
                    xdate=False, ydate=True,
                    marker=mpl_combinations[marker_style_index][0],
                    color=mpl_combinations[marker_style_index][1])
            marker_style_index += 1
    
        if plotted_host is False:
            axHisty2.xaxis.set_major_formatter(nullfmt)

        sn_points, sn_dates = unzip(data.snps)
        sn_dates = mpl.dates.date2num(sn_dates)
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
        del fig

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('database',
            help="path to sqlite database containing information about phage")
    parser.add_argument('-c', '--contigs', dest='contigs_file', type=argparse.FileType('r'),
            help='A file containing contigs to consider.')
    parser.add_argument('-q', '--snp-quality', dest='snp_quality', default=20,
            help="Minimum quality for a SNP")
    parser.add_argument('-s', '--snp', dest='snps',
            help="vcf file containing SNPs")
    parser.add_argument('-o', '--output', dest='output', default='.',
            help="output directory for the image")
    parser.add_argument('contigs', nargs='*', help='Name of contigs to consider')

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
