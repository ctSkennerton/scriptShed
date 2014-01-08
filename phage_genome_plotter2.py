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
from __future__ import division, print_function
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


class PlotFormatter(object):
    """docstring for PlotFormatter"""
    def __init__(self):
        super(PlotFormatter, self).__init__()
        self.combinations = []
        self.index = 0
        mpl_line_styles = ['-' , '--' , '-.' , ':']
        mpl_line_marker_styles = [ 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd']
        mpl_line_color_styles = ['b', 'g', 'r', 'c', 'm', 'y']
        mpl_combinations = []
        for k in mpl_line_styles:
            for i in mpl_line_marker_styles:
                for j in mpl_line_color_styles:
                    self.combinations.append((i, j, k))

    def __len__(self):
        return len(self.combinations)

    def __call__(self):
        return self.format()

    def format(self):
        f = self.combinations[self.index]
        self.index += 1
        if self.index >= len(self):
            self.index = 0
        return f




class Crispr(object):
    """holds many spacers and some metadata about formatting for plots"""
    def __init__(self, format):
        super(Crispr, self).__init__()
        self.spacer_clusters = dict()
        self.host = None
        self.format = format

        self.marker = self.format[0]
        self.line = self.format[2]
        self.colour = self.format[1]
    
    def __len__(self):
        return len(self.spacer_clusters)
    
    def __getitem__(self, key):
        return self.spacer_clusters[key]
    
    def __setitem__(self, key, value):
        self.spacer_clusters[key] = value
    
    def __delitem__(self, key):
        del self.spacer_clusters[key]
    
    def __contains__(self, item):
        return item in self.spacer_clusters  
    
    def __iter__(self):
        return iter(self.spacer_clusters)  

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
    def __init__(self, spid=None, start=0, end=0, contig=None, crispr=None, cluster=None, host=None, timepoint=None, cov=None):
        if start > end:
            self.start = end
            self.end = start
        else:
            self.start = start
            self.end = end
        self.id = spid
        self.contig = contig
        self.timepoint = timepoint
        self.host = host
        self.crispr = crispr
        self.cluster = cluster
        self.cov = cov

class Contig(object):
    def __init__(self, Id, Name, length, timepoint=None):
        super(Contig, self).__init__()
        self.id = Id
        self.Name = Name
        self.length = length
        self.timepoint = timepoint
        # list of 2-tuple of the positions of and timepoints of snps and
        # spacers

        self.snps = []
        self.spacers = {}
        self._relabs = 0.0
        
    def __len__(self):
        return self.length
    
    @property
    def relative_abundance(self):
        return self._relabs
    
    @relative_abundance.setter
    def relative_abundance(self, value):
        self._relabs = value
        
    def get_spacers(self):
        ret = {} #defaultdict(dict)
        for sp in self.spacers.values():
            if sp.crispr in ret:
                ret[sp.crispr]['host'] = sp.host
                ret[sp.crispr]['data'].append([sp.start, sp.timepoint])
            else:
                ret[sp.crispr] = dict()
                ret[sp.crispr]['host'] = sp.host
                ret[sp.crispr]['data'] = []
                ret[sp.crispr]['data'].append([sp.start, sp.timepoint])
                
        return ret

    def get_spacer_ids(self):
        ret = set()
        for spid in self.spacers.keys():
            ret.add(spid)
        return ret

    def get_spacer_positions(self):
        ret = {}
        for sp in self.spacers.values():
            ret[sp.cluster] = sp.start
        return ret

        
    def has_spacer_cluster(self, spid):
        if spid in self.spacers:
            return True
        return False

    def get_spacer_clusters(self):
        ret = set()
        for sp in self.spacers.values():
            ret.add(sp.cluster)
        return ret

    def get_crisprs(self, formatter):
        """Return a mapping of crispr clusters and their spacers as Crispr objects"""
        ret = {}
        for sp in self.spacers.values():
            if sp.crispr not in ret:
                ret[sp.crispr] = Crispr(formatter.format())
            ret[sp.crispr][sp.id] = sp
            ret[sp.crispr].host = sp.host
        return ret

    def get_spacer_from_id(self, spid):
        try:
            return self.spacers[spid]
        except KeyError:
            return None

    def get_spacers_from_spacer_cluster(self, cid):
        for sp in self.spacers.values():
            if sp.cluster == cid:
                yield sp


class ScaffoldFragment(Contig):
    def __init__(self, Id, Name, start, Length, complement, timepoint=None):
        super(ScaffoldFragment, self).__init__(Id, Name, Length,
                timepoint=timepoint)
        self.start = start
        self.complement = complement
        
class Scaffold(object):
    def __init__(self, Id, Name, Length, timepoint=None):
        super(Scaffold, self).__init__()
        self.id = Id
        self.name = Name
        self.length = Length
        self.fragments = {}
        self.timepoint = timepoint
        self.snps = []

    def __len__(self):
        return self.length
    
    @property
    def relative_abundance(self):
        sum = 0.0
        for fragid, fragobj in self.fragments.items():
            sum += fragobj.relative_abundance
        frag_count = len(self.fragments)
        return sum / len(self.fragments)
    
    def get_spacers(self):
        ret = {}  #defaultdict(dict)
        for frag in self.fragments.values():
            for sp in frag.spacers.values():
                if frag.complement:
                    spacer_pos = sp.start + frag.start - frag.length
                else:
                    spacer_pos = sp.start + frag.start
                    
                if sp.crispr in ret:
                    ret[sp.crispr]['host'] = sp.host

                    ret[sp.crispr]['data'].append([spacer_pos, sp.timepoint])
                else:
                    ret[sp.crispr] = dict()
                    ret[sp.crispr]['host'] = sp.host
                    ret[sp.crispr]['data'] = []
                    ret[sp.crispr]['data'].append([spacer_pos, sp.timepoint])
 
        return ret   

    def get_spacer_positions(self):
        ret = {}
        for frag in self.fragments.values():
            for sp in frag.spacers.values():
                if frag.complement:
                    spacer_pos = sp.start + frag.start - frag.length
                else:
                    spacer_pos = sp.start + frag.start

                ret[sp.cluster] = spacer_pos
        return ret

    def get_spacer_ids(self):
        ret = set()
        for frag in self.fragments.values():
            for spid in frag.spacers.keys():
                ret.add(spid)
        return ret
              
    def has_spacer_cluster(self, spid):
        for frag in self.fragments.values():
            if spid in frag.spacers:
                return True
        return False

    def get_spacer_clusters(self):
        ret = set()
        for frag in self.fragments.values():
            for sp in frag.spacers.values():
                ret.add(sp.cluster)
        return ret

    def get_crisprs(self, formatter):
        """Return a mapping of crispr clusters and their spacers as Crispr objects"""
        ret = {}
        for frag in self.fragments.values():
            for sp in frag.spacers.values():
                if sp.crispr not in ret:
                    ret[sp.crispr] = Crispr(formatter.format())
                ret[sp.crispr][sp.id] = sp
                ret[sp.crispr].host = sp.host
        return ret

    def get_spacer_from_id(self, spid):
        for frag in self.fragments.values():
            try:
                return frag.spacers[spid]
            except KeyError:
                pass
        return None

    def get_spacers_from_spacer_cluster(self, cid):
        for frag in self.fragments.values():
            for sp in frag.spacers.values():
                if sp.cluster == cid:
                    yield sp



class ContigCluster(object):
    ''' Hold all the contigs for a cluster
        The main job of this class is to hold some metadata about contig
        clusters such as the 'reference' contig, which is the one that is
        called the reference in the VCF file for this cluster, from which the
        snp and spacer positions will be determined.
    '''
    def __init__(self):
        self.contigs = {}
        self.reference = None
        self.spacer_matrix = {}
        self.spacer_pos = {}
    
    def __len__(self):
        return len(self.contigs)
    
    def __getitem__(self, key):
        return self.contigs[key]
    
    def __setitem__(self, key, value):
        self.contigs[key] = value
    
    def __delitem__(self, key):
        del self.contigs[key]
    
    def __contains__(self, item):
        return item in self.contigs  
    
    def __iter__(self):
        return iter(self.contigs)   
    
    def values(self):
        return self.contigs.values()

    def keys(self):
        return self.contigs.keys()

    def get_contigs_times(self, timepoint):
        for contig in self.contigs.values():
            if contig.timepoint == timepoint:
                yield contig
        #return [contig if contig.timepoint == timepoint for contig in self.contigs.values()]
            
    
    def populate_spacer_matrix(self, sorted_times):
        reference_spacers = self.reference.get_spacer_ids()
        reference_spacer_clusters = self.reference.get_spacer_clusters()
        self.spacer_pos = self.reference.get_spacer_positions()

        for spid in reference_spacers:
            spacer = self.reference.get_spacer_from_id(spid)
            if spacer.crispr not in self.spacer_matrix:
                self.spacer_matrix[spacer.crispr] = dict()

            if self.reference.timepoint not in self.spacer_matrix[spacer.crispr]:
                self.spacer_matrix[spacer.crispr][self.reference.timepoint] = set()
            
            self.spacer_matrix[spacer.crispr][self.reference.timepoint].add(spacer.cluster)
            if spacer.host is not None:
                #print "Host for crispr %s is %s" % (str(spacer.crispr), str(spacer.host))
                self.spacer_matrix[spacer.crispr]['host'] = spacer.host
            #self.spacer_pos[spacer.cluster] = spacer.start

        for contig in self.contigs.values():
            if contig.timepoint == self.reference.timepoint:
                continue
            spacer_clusters_for_current = reference_spacer_clusters & contig.get_spacer_clusters()
            for spc in spacer_clusters_for_current:
                for spacer in contig.get_spacers_from_spacer_cluster(spc):
                    if spacer.crispr not in self.spacer_matrix:
                        #self.spacer_matrix[spacer.crispr] = dict()
                        #print "WARNING: CRISPR %s in timepoint %s not in reference %s" % (str(spacer.crispr), str(contig.timepoint), str(self.reference.name))  
                        continue                

                    if contig.timepoint not in self.spacer_matrix[spacer.crispr]:
                        self.spacer_matrix[spacer.crispr][contig.timepoint] = set()

                    self.spacer_matrix[spacer.crispr][contig.timepoint].add(spacer.cluster)

    @property
    def relative_abundance(self):
        times = defaultdict(list)
        for contig_data in self.contigs.values():
            times[contig_data.timepoint].append(contig_data.relative_abundance)
        
        for t, data in times.items():
            s = sum(data)
            times[t] = s / len(data)
        sorted_times = sorted(times.keys())
        sorted_points = [times[x] for x in sorted_times ]
        return sorted_times, sorted_points


def generate_placeholders(l):
    placeholder= '?' # For SQLite. See DBAPI paramstyle.
    placeholders= ', '.join(placeholder * len(l))
    return placeholders            

def get_contig_length(cur, name):
    cur.execute('''SELECT Length FROM contigs WHERE Id = ?''', (name,))
    result = cur.fetchone()
    if not result:
        raise RuntimeError('could not find contigs name %s in the database' % name)
    return result[0]


def get_spacers_for_contig(cur, name):
    cur.execute('''CREATE TEMP VIEW IF NOT EXISTS spacer_match_times AS 
                SELECT
                    ContigID, 
                    SpacerID, 
                    SpacerCluster, 
                    SpacerTp, 
                    contigs.Timepoint AS ContigTp, 
                    Start, 
                    End,
                    CrisprClusterID,
                    HostID
                FROM 
                    contigs 
                JOIN 
                    (
                    SELECT 
                        spacer_matches.ContigID, 
                        crispr_host_connector.SpacerID, 
                        SpacerCluster, 
                        SpacerTp, 
                        spacer_matches.Start, 
                        spacer_matches.End,
                        crispr_host_connector.CrisprClusterID,
                        crispr_host_connector.HostID 
                    FROM 
                        crispr_host_connector 
                    JOIN 
                        spacer_matches 
                    ON 
                        crispr_host_connector.SpacerID = spacer_matches.SpacerID
                    ) 
                ON 
                    ContigID = contigs.Id'''
            )


    cur.execute('''select CrisprClusterID, HostID, Start, End, ContigTp, SpacerID, SpacerCluster from
            spacer_match_times where ContigID = ?''', (name,))

    result = cur.fetchall()
    if not result:
        return dict()
        #raise RuntimeError('problem with returning spacers for %s' % name)
    ret = {}
    for row in result:
        ret[row[5]] = Spacer(row[5], row[2], row[3], name, row[0], host=row[1], timepoint=row[4], cluster=row[6])
    return ret

def get_relative_abundance(cur, name):
    cur.execute(''' SELECT 
                        Relab 
                    FROM 
                        times 
                    JOIN 
                        (
                        SELECT 
                            Timepoint AS Tp, 
                            Relab 
                        FROM
                            contig_relabs
                        WHERE 
                            ContigID = ?
                        ) On times.Id=Tp''', (name,))

    result = cur.fetchone()
    if not result:
        raise RuntimeError('cannot get the relative abundance of %s' % name)
    #print result
    return result[0]

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
            
        
def generate_contig_query_list(name, name_file):
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
    return query

def get_cluster_id_from_contigs(cur, name=None, name_file=None):
    cl = generate_contig_query_list(name, name_file)
    query = "SELECT\
                Cluster, Name\
            FROM\
                contigs\
            WHERE\
                Name IN (%s)" % generate_placeholders(cl)
    cur.execute(query, tuple(cl))
    result = cur.fetchall()
    if not result:
        raise RuntimeError("Cannot identify any clusters for given contigs")
    else:
        ret = {}
        for row in result:
            if row[0] not in ret:
                ret[row[0]] = ContigCluster()
                ret[row[0]].reference = row[1]
        #ret = [ x[0] if x[0] for x in result ]
        return ret

def get_contigs_for_clusters(cur, clusters):
    
    query ="SELECT \
                contigs.Name,\
                contigs.Id, \
                contigs.Flag,\
                contigs.Length,\
                contigs.Cluster,\
                times.Timepoint\
            FROM \
                contigs\
            JOIN\
                times\
            ON\
                contigs.Timepoint = times.Id\
            WHERE \
                Cluster IN (%s)" % generate_placeholders(clusters.keys())

    cur.execute(query, tuple(clusters))
    result = cur.fetchall()
    for row in result:
        if row[2] & 4:
            clusters[row[4]][row[0]] = Scaffold(row[1], row[0], row[3], timepoint=row[5])
            clusters[row[4]][row[0]].fragments = get_contigs_for_scaffold(cur, row[1])
        else:
            clusters[row[4]][row[0]] = Contig(row[1], row[0], row[3], timepoint=row[5])

        if row[0] == clusters[row[4]].reference:
            clusters[row[4]].reference = clusters[row[4]][row[0]]

    return clusters


def generate_contig_ids(cur, name=None, name_file=None):
    id_map = {}
    query = generate_contig_query_list(name, name_file)
    try:
        cur.execute('''SELECT 
                            contigs.Name, 
                            contigs.Id, 
                            contigs.Flag,
                            contigs.Length,
                            contigs.Cluster,
                            times.Timepoint
                        FROM 
                            contigs
                        JOIN
                            times
                        ON
                            contigs.Timepoint = times.Id
                        WHERE 
                            Name IN (?)''', tuple(query))
        #cur.execute('''SELECT Name, Id, Flag, Length FROM contigs WHERE Name IN (?)''', tuple(query))
    except sqlite3.ProgrammingError:
        print(query)
        return {}
    result = cur.fetchall()
    if not result:
        raise RuntimeError('wierd')
    for row in result:
        if row[2] & 4:
            id_map[row[0]] = Scaffold(row[1], row[0], row[3], timepoint=row[5])
            id_map[row[0]].fragments = get_contigs_for_scaffold(cur, row[1])
        else:
            id_map[row[0]] = Contig(row[1], row[0], row[3], timepoint=row[5])
             
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
        raise RuntimeError('cannot get the relative abundance of %s' % host)
    return result

def unzip(data):
    v, h = zip(*data)
    v = [x if x is not None else 0 for x in v ]
    return v, h
    
def set_plot(sorted_dates, title=None, yscale=None):
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
    axHisty.yaxis.set_major_formatter(NullFormatter())
    
    axHisty.set_ylim((mpl.dates.date2num(sorted_dates[0]),
        mpl.dates.date2num(sorted_dates[-1])))
        
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    
    axHisty2 = axHisty.twiny()
    axHisty2.set_ylim((mpl.dates.date2num(sorted_dates[0]),
        mpl.dates.date2num(sorted_dates[-1])))
    
    axHisty2.xaxis.set_major_locator(MaxNLocator(4))
    axHisty2.yaxis.set_major_formatter(NullFormatter())

    labels = axHisty2.get_xticklabels()
    for label in labels:
        label.set_rotation(-30)
        label.set_horizontalalignment('right')
        label.set_size('small')
    
    return fig, axScatter, axHisty, axHisty2
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################



def doWork( args ):
    """ Main wrapper"""
    conn=sqlite3.connect(args.database, detect_types=sqlite3.PARSE_DECLTYPES)
    cur=conn.cursor()
    sorted_dates = get_all_timepoints(cur)

    clusters = get_cluster_id_from_contigs(cur, args.contigs, args.contigs_file)
    clusters = get_contigs_for_clusters(cur, clusters)
    for cluster_name, contig_data in clusters.items():
        for contigobj in contig_data.values():
            if isinstance(contigobj, Scaffold):
                for fragid, fragobj in contigobj.fragments.items():
                    fragobj.spacers.update(get_spacers_for_contig(cur, fragid))
                    fragobj.relative_abundance = get_relative_abundance(cur, fragid)
            else:
                contigobj.spacers.update(get_spacers_for_contig(cur, contigobj.id))
                contigobj.relative_abundance = get_relative_abundance(cur, contigobj.id)

    cluster_reference = None
    if os.path.exists(args.snps) and os.stat(args.snps).st_size != 0:
        vcf_reader = vcf.Reader(open(args.snps))
        for record in vcf_reader:
            if cluster_reference is None:
                for clust_id, clust in clusters.items():
                    try:
                        cluster_reference = clust[record.CHROM]
                        clusters[clust_id].reference = clust[record.CHROM]
                        #print clusters[clust_id].reference
                        #print cluster_reference
                    except KeyError:
                        print("Cannot find %s for cluster" % (record.CHROM, ))
                
            if record.QUAL is not None and int(record.QUAL) < int(args.snp_quality):
                continue
            if cluster_reference is None:
                continue
            for sample in record.samples:
                if sample['GT'] != '0/0' and sample['GT'] != './.' and sample['GT'] != '0':
                    cluster_reference.snps.append([record.POS,
                        names_to_times[sample.sample]])

    #-----
    # make a 2d plot  

    for cluster_id, cluster_data in clusters.items():
        fig, axScatter, axHisty, axHisty2 = set_plot(sorted_dates)
        ra_time, ra_point = cluster_data.relative_abundance
        ra_time = mpl.dates.date2num(ra_time)
        #print(ra_time, ra_point, sep="\t")

        axHisty.plot_date(ra_point, ra_time, xdate=False,
                ydate=True, color='0.7', marker=' ', linestyle='--')


        axScatter.set_xlim([0, len(cluster_data.reference)])
        axScatter.set_ylim((mpl.dates.date2num(sorted_dates[0]),
            mpl.dates.date2num(sorted_dates[-1])))
        axScatter.set_xlabel('genome position (bp)')
        axScatter.set_title("Phage %s" % (str(cluster_id)))
                
        plotted_host = False
        marker_style_index = 0
    
        formatter = PlotFormatter()

        cluster_data.populate_spacer_matrix(sorted_dates)

        for crispr_cluster, crispr_data in cluster_data.spacer_matrix.items():
            #print(crispr_cluster, crispr_data, sep=' ')
            format = formatter()
            try:
                if crispr_data['host'] is not None:
                    #print "Crispr %s has host %s" % (str(crispr_cluster), str(crispr_data['host']))
                    host_relab_data = calc_host_relab(cur, crispr_data['host'])
                    v, h = unzip(host_relab_data)
                    h = list(h)
                    h = mpl.dates.date2num(h)

                    axHisty2.plot_date(v, h, xdate=False, ydate=True,
                            label=str(crispr_data['host']), alpha=0.5,
                            ls=format[2],
                            marker=format[0],
                            color=format[1])
                    plotted_host = True
            except KeyError:
                pass

            for timepoint, spacers in crispr_data.items():
                if timepoint == 'host':
                    continue
                sp_dates = []
                sp_points = []
                for spcid in spacers:
                    sp_dates.append(timepoint)
                    try:
                        sp_points.append(cluster_data.spacer_pos[spcid])
                    except TypeError, e:
                        print(e, file=sys.stderr)
                        print(cluster_data.spacer_pos, spcid, sep=' ')
                        print(cluster_id, cluster_data.reference, sep="\t")
                        for timepoint, spacers in crispr_data.items():
                            for spcid in spacers:
                                print(timepoint, spcid)
                        sys.exit(1)


                sp_dates = mpl.dates.date2num(sp_dates)
                axScatter.plot_date(sp_points, sp_dates, label=crispr_cluster, alpha=0.5,
                    xdate=False, ydate=True,
                    marker=format[0],
                    color=format[1],
                    ms=5)

    
        if plotted_host is False:
            axHisty2.xaxis.set_major_formatter(NullFormatter())
            
        if len(cluster_data.reference.snps):
            sn_points, sn_dates = unzip(cluster_data.reference.snps)
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
            
        plt.savefig(os.path.join(args.output, "cluster_%d.png" % cluster_id), dpi=300,
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
