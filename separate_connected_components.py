#!/usr/bin/env python
###############################################################################
#
# separate_connected_components 
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
import os
import networkx as nx
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

def doWork( args ):
    if args.outdir is None:
        outdir = './'
    else:
        outdir = args.outdir

    reader = getattr(nx, "read_"+args.format)
    G = reader(args.infile)
    cc = nx.connected_components(G)
    for i, c in enumerate(cc):
        if len(c) == 1:
            break
        with open(outdir+'/component_'+str(i), 'w') as fp:
            for n in c:
                print(n, file=fp)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help="graph file in gexf format for"
            " partiitioning")
    parser.add_argument('outdir', nargs='?', 
            help="output directory name for separated components.  If not given"
            " the current directory will be used.")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    parser.add_argument('-f', '--format', default='gexf', help="format of the
            graph: gexf, gml")
    
    # parse the arguments
    args = parser.parse_args()        

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

