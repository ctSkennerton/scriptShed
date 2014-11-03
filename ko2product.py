#!/usr/bin/env python
import requests
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=argparse.FileType(),
            help="Input file containing two columns: the first containing"
            " the name of the gene and the second column containing the KO"
            " identifier")
    args = parser.parse_args()
    for line in args.infile:
        fields = line.rstrip().split()
        if len(fields) == 2:
            r = requests.get('http://rest.kegg.jp/find/ko/'+fields[1])
            print fields[0]+"\t"+r.text.rstrip()
        else:
            print fields[0]+"\t\thypothetical protein"
