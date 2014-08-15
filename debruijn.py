#!/usr/bin/env python
from __future__ import print_function, division
import sys
import argparse
from Bio import SeqIO
import networkx as nx
from collections import deque


class CheckOdd(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            i = int(values)
            if i % 2 != 1:
                raise
            else:
                setattr(namespace, self.dest, i)
        except:
            print("The argument to %s must be an odd integer number" % option_string)
            sys.exit(1)

def collapse_linear_paths(graph):
    seen_nodes = set()
    nodes = graph.nodes()
    nodes_to_delete = set()
    for i in range(len(nodes)):
        path_end = None
        path_begin = None
        current_path = deque()
        node = nodes[i]
        if node in seen_nodes:
            continue

        if graph.in_degree(node) == 1 and graph.out_degree(node) == 1:
            # we have the start of a simple path
            current_path.append(node)
            root_node = node
            seen_nodes.add(node)
            while 1:
                sucs = graph.successors(node)
                preds = graph.predecessors(node)
                if len(sucs) == 1 and len(preds) == 1:
                    node = sucs[0]
                    current_path.append(node)
                    seen_nodes.add(node)
                else:
                    path_end = sucs
                    for i in sucs:
                        seen_nodes.add(i)
                    break

            node = root_node
            while 1:
                sucs = graph.successors(node)
                preds = graph.predecessors(node)
                if len(preds) == 1 and len(sucs) == 1:
                    node = preds[0]
                    current_path.appendleft(node)
                    seen_nodes.add(node)
                else:
                    path_begin = preds
                    for i in preds:
                        seen_nodes.add(i)
                    break

            # now join the sequence,
            # make a new node and link it to the end of the path
            seq = []
            if len(current_path) > 1:
                for i, s in enumerate(current_path):
                    if i == 0:
                        seq.append(s)
                    else:
                        seq.append(s[-1])
                    if s in nodes_to_delete:
                        print("Error node %s seen more than once" % s)
                    nodes_to_delete.add(s)

                seq = ''.join(seq)
                graph.add_node(seq)
                for i in path_end:
                    graph.add_edge(seq, i)

                for i in path_begin:
                    graph.add_edge(i, seq)

    for i in nodes_to_delete:
        graph.remove_node(i)


def consume_reads(graph, fastaFile, readFormat, k, countMax):
    for count, record in enumerate(SeqIO.parse(fastaFile, readFormat)):
        if countMax is not None and count > countMax:
            break

        prev_kmer=None
        rc = record.reverse_complement()
        if str(rc.seq) < str(record.seq):
            record = rc

        for i in range(len(record.seq) - k):
            kmer = str(record.seq[i:i+k])
            graph.add_node(kmer)
            if prev_kmer is not None:
                graph.add_edge(prev_kmer, kmer)

            prev_kmer = kmer



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', type=argparse.FileType(),
            help="input fasta file to generate debruijn graph from")
    parser.add_argument('-o', '--outfile',
            help='Output file for generating gexf file of graph')
    parser.add_argument('-f', '--format',
            help='input read format', default='fasta')
    parser.add_argument('-k', '--kmer', action=CheckOdd, default=63,
            help='kmer size for the graph')
    parser.add_argument('-m', '--max', default=None, type=int, help='read this meny records from the file')
    parser.add_argument('--collapse', default=False, action='store_true',
            help='collapse linear paths into a single node')
    args = parser.parse_args()

    G = nx.DiGraph()
    consume_reads(G, args.infile, args.format, args.kmer, args.max)
    if args.collapse:
        collapse_linear_paths(G)
    nx.write_gexf(G, args.outfile)
