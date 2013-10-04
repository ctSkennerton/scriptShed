#!/usr/bin/env python
###############################################################################
#
# plot_seq_record 
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
__copyright__ = "Copyright 2013"
__credits__ = ["uqcskenn"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "uqcskenn"
__email__ = "uqcskenn@hawke"
__status__ = "Development"

###############################################################################

import argparse
import sys

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')     
from biograpy import Panel, tracks, features
from Bio import SeqIO
from BCBio import GFF
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

def doWork( args ):
    panel=Panel(fig_width=900, padding = 25, grid=None, xmin=0)
    seq_length = 0
    for gff in args.gffs:
        seqrecord = GFF.parse(gff).next()
        if len(seqrecord) > seq_length:
            seq_length = len(seqrecord)
        #seqrecord = SeqIO.parse(args.infile, "genbank").next()
        cds_track = tracks.BaseTrack(sort_by = 'collapse')
        for feature in seqrecord.features:
            if feature.type == 'CDS':
                #print feature.qualifiers['product']
                if feature.qualifiers['product'][0] == 'hypothetical protein':
                    col = '#BDBDBD'
                else:
                    col = '#2B8CBE'
                feat = features.GenericSeqFeature(feature, color_by_cm=False,
                        fc=col )
                cds_track.append(feat)
            elif feature.type == 'source':
                cds_track.append(features.GenericSeqFeature(feature,
                    color_by_cm=False, alpha=0.0, fc='1.0', ec='1.0'))
            else:
                cds_track.append(features.GenericSeqFeature(feature,
                    color_by_cm=False, fc='0.0', ec='0.0'))
        panel.add_track(cds_track)
    panel.save(args.outfile, xmin=0,xmax=seq_length)
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile',  help="name of output image")
    parser.add_argument('gffs', nargs='+', type=argparse.FileType('r'),
            help="gff file")
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

