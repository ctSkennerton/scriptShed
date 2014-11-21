#!/usr/bin/env python
from __future__ import print_function
import requests
import argparse
import itertools
import StringIO

def grouper(n, iterable):
    it = iter(iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           return
       yield chunk

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given a list of kegg orthology identifiers, '\
            'will query the kegg server and download all genes in kegg that have that KO')
    parser.add_argument('-n', '--ntseq', help='download as nucleotide sequences, the default is protein sequences')
    parser.add_argument('ko', nargs='+', help='List of kegg orthology identifiers to download')
    args = parser.parse_args()

    ko_link = 'http://rest.kegg.jp/link/genes/'
    ko_get = 'http://rest.kegg.jp/get/'
    seqtype = '/aaseq'
    if args.ntseq is not None:
        seqtype = '/ntseq'
    for i in grouper(10, args.ko):
        ko_request='+'.join(i)
        r = requests.get(ko_link+ko_request)
        str_file = StringIO.StringIO(r.text)
        genes = []
        for line in str_file:
            ko_id, gene = line.rstrip().split()
            genes.append(gene)

        for j in grouper(10, genes):
            gene_request = '+'.join(j)
            r = requests.get(ko_get+gene_request+seqtype)
            print(r.text, end='')

