import argparse
from ftplib import FTP
import re


def acc_to_ftp_path(acc):
    match = re.search('(\w+)_(\d+)\.(\d)', acc)
    if match:
        prefix = match.group(1)
        accn = match.group(2)
        version = match.group(3)
        accp = re.findall('.{3}', accn)
        accp = '/'.join(accp)
        return (prefix, accp, version)
    else:
        raise ValueError("could not get FTP path from ".format(acc))


def get_files_from_types(types, base_name, ftp, path_only=False):
    for g in ftp.mlsd():
        for t in types:
            if g[0] == base_name + t:
            #if g[0].endswith(t):
                with open(g[0], 'wb') as ofp:
                    if path_only:
                        print("ftp://ftp.ncbi.nlm.nih.gov/{}/{}".format(ftp.pwd(), g[0]))
                    else:
                        res = ftp.retrbinary('RETR {}'.format(g[0]), ofp.write)
                        print("Downloaded {} to current directory".format(g[0]))


def map_types_to_file_suffix(types):
    mapping = {'fna': '_genomic.fna.gz',
               'faa': '_protein.faa.gz',
               'ffn': '_cds_from_genomic.fna.gz',
               'gb': '_genomic.gbff.gz',
               'gff': '_genomic.gff.gz'
              }
    return [mapping[i] for i in types]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='file containing accession numbers from NCBI assembly in the form of GCA_001871445.1, one per line')
    parser.add_argument('-f', '--outfmt', required=True, help="select the type of data to download. Can be specified multiple times if more than one file type is required. fna: genomic fasta; faa: translated protein sequences; gb: genbank; gff: gff file of annotations; ffn: coding regions, not translated", choices=['fna', 'faa', 'gb', 'gff', 'ffn'], action='append')
    parser.add_argument('-p', '--path-only', dest='path_only', default=False, action='store_true', help='only print the file paths, don\'t actually download anything')
    args = parser.parse_args()

    types = map_types_to_file_suffix(args.outfmt)
    with open(args.infile) as fp, FTP("ftp.ncbi.nlm.nih.gov") as ftp:
        ftp.login()
        for line in fp:
            accession = line.rstrip()
            prefix, accp, version = acc_to_ftp_path(accession)
            base_path = "genomes/all/{}/{}".format(prefix,accp)
            ftp.cwd(base_path)
            for f in ftp.mlsd():
                if accession in f[0]:
                    ftp.cwd(f[0])
                    get_files_from_types(types, f[0], ftp, args.path_only)
            ftp.cwd("/")
