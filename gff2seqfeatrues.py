#!/usr/bin/env python
from __future__ import print_function
import sys
from Bio import SeqIO
from BCBio import GFF
import argparse

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("gff")
    p.add_argument('infasta')
    #p.add_argument('outffn')
    args = p.parse_args()

    seq_dict = SeqIO.to_dict(SeqIO.parse(args.infasta, "fasta"))
    with open(args.gff) as gffp:
        for rec in GFF.parse(gffp, base_dict=seq_dict):
            for feature in rec.features:
                feat = feature.extract(rec.seq)
                try:
                    print(">{}".format(feature.qualifiers['locus_tag'][0]))
                    print(feat)
                except KeyError as e:
                    print("skipping feature as it doesn't have a locus tag", file=sys.stderr)
                    print(feature, file=sys.stderr)
