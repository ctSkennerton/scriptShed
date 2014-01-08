#!/usr/bin/env python
###############################################################################
#
# __script_name__.py - description!
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

import os
import errno

import numpy as np
np.seterr(all='raise')

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText



###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################
def make_stats_box(data, loc=2):
    avg = np.mean(data)
    count = len(data)
    stats = "A: %.1f\nC: %d" % (avg, count)
    fp = dict(size=6)
    return AnchoredText(stats, loc=loc, prop=fp)

def plot_pair(ax, data, title=None, xtitle=None, ytitle=None, stats=True):
    ax.hist(data, bins=100)#, histtype='step', alpha=1.0,
            #lw=0.3)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(6))
    if title:
        ax.set_title(title)
    if xtitle:
        ax.set_xlabel(xtitle)
    if ytitle:
        ax.set_ylabel(ytitle)
    if stats:
        ax.add_artist(make_stats_box(data))

def matrix_plots(data, nrow, ncol, genome_names, output, xaxismin=50,
        xaxismax=105):
    fig, subplots = plt.subplots(nrows=nrow, ncols=ncol, sharex=True,
            figsize=(10,8))
    row_counter = 0
    data_counter = len(genome_names) - 1
    for row in subplots:
        col_counter = 0
        first_col = True

        for ax in row:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize('xx-small')
                tick.label.set_rotation('vertical')
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize('xx-small')
            if first_col:
                ax.set_ylabel(genome_names[data_counter][0])
                first_col = False
            if col_counter <= row_counter:
                xtitle=None
                if row_counter == nrow - 1:
                    xtitle = genome_names[data_counter][1]
                ax.set_xlim(xaxismin, xaxismax)
                plot_pair(ax, data[data_counter], xtitle=xtitle)
                data_counter -= 1
            else:
                ax.get_yaxis().set_ticks([])
                ax.set_frame_on(False)
                for t in ax.xaxis.get_ticklines():
                    t.set_visible(False)
                #ax.set_axis_off()
            col_counter += 1
        row_counter += 1
    fig.subplots_adjust(hspace=0.1, wspace=0.35, left=0.08, bottom=0.08,
            top=0.9, right=0.9 )
    #plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches=0)

    #-----
    # clean up!
    plt.close(fig)
    del fig


def generate_matrix(args):
    genomes = set()
    aai_results = []
    names = []
    try:
        with open(args.filename, "r") as fh:
            for line in fh:
                line = line.rstrip()
                fields = line.split(",")
                name_pair = fields[0].split("_")
                names.append(name_pair)
                genomes.add(name_pair[1])
                genomes.add(name_pair[0])
                aai_results.append(map(float, fields[1:]))
    except ValueError, e:
        raise e

    ncol = len(genomes) - 1
    matrix_plots(aai_results, ncol, ncol, names, args.output,
            xaxismin=args.xmin, xaxismax=args.xmax)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-x' '--xmin', type=int, default=50, help='minimum '\
            'value for the x-axis', dest='xmin')
    parser.add_argument('-X' '--xmax', type=int, default=105, help='maximum '\
            'value for the x-axis', dest='xmax')
    parser.add_argument('filename', help="CSV file containing ANI results")
    parser.add_argument('output', help="Output image. Remember to include the "\
                        "file extension: valid options are: pdf, png, svg, ps")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    generate_matrix(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

