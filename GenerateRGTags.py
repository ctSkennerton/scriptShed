import os
import sys
import glob
import pysam
import argparse
import multiprocessing

def get_args():
    '''Parse sys.argv'''
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True,
                        help='The input bam file.')
    parser.add_argument('-CN',type=str,
                        help="Name of sequencing center producing the read. \
                        GATK not required.")
    parser.add_argument('-DS',type=str,
                        help="Description. GATK Not Required. ")
    parser.add_argument('-DT',type=str,
                        help="Date the run was produced (ISO8601 date or date/time). \
                        GATK Not Required. ")
    parser.add_argument('-PI',type=int,
                        help="Predicted median insert size. GATK Not Required.")
    parser.add_argument('-PL',type=str, required=True,
                        choices = ['CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 
                                    'HELICOS', 'IONTORRENT', 'PACBIO'],
                        help="Platform/technology used to produce the reads.")
    parser.add_argument('file', metavar='FILE', type=str, nargs='+', 
                        help="files containing read headers for each of the \
                                conditions")

    args = parser.parse_args()
    return args

def read2RGId(filename, readMap):
    """Given a list of files containing read headers. Generate a map of those
    headers to a corresponding RG identifier"""
    f = open(filename, 'r')
    for line in f:
        readMap[line.rstrip()] = filename
    

def addRG2Header(readMap, files, args):
    """Add read group info to a header."""

    samfile = pysam.Samfile(args.input, 'r')
    new_header = samfile.header.copy()
    samfile.close()
    new_header['RG'] = []
    for fn in files:
        # process the reads to get their RG
        read2RGId(fn, readMap)

        # Create the RG tags in the header
        # CREATE TEMPLATE
        # Read group. Unordered multiple @RG lines are allowed.
        RG_template = { 'ID': '',           # Read group identifier. e.g., Illumina flowcell + lane name and number
                        'CN': '',           # GATK Not Required. Name of sequencing center producing the read.
                        'DS': '',           # GATK Not Required. Description
                        'DT': '',           # GATK Not Required. Date the run was produced (ISO8601 date YYYY-MM-DD or YYYYMMDD)
                        'PI': '',           # GATK Not Required. Predicted median insert size.
                        'PU': '',           # GATK Not Required. Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD).
                        'SM': '',           # Sample. Use pool name where a pool is being sequenced.
                        'PL': 'ILLUMINA'}   # Platform/technology used to produce the reads.

        # ADD INFO TO TEMPLATE
        RG_template['ID'] = fn
        if args.CN: RG_template['CN'] = args.CN.upper()
        if args.DS: RG_template['DS'] = args.DS
        if args.DT: RG_template['DT'] = args.DT
        if args.PI: RG_template['PI'] = args.PI
        if args.PL: RG_template['PL'] = args.PL

        #RG_template['LB'] = sam_info["sample_name"]
        #RG_template['SM'] = sam_info["sample_name"]
        #RG_template['DS'] = "{0}.{1}".format(sam_info['sample_name'], sam_info['locality'])
        #RG_template['PU'] = '{0}.{1}'.format(sam_info['flowcell_id'], sam_info['lane'])
        print RG_template
        new_header['RG'].append( RG_template )
    print new_header['RG']
    return new_header


def add_RGs_2_BAMs_runner(newHeader, readMap, args):
    """Generates the correct @RG header and adds a RG field to a bam file."""
    # Make Bam of Sam if it doesn't exist
    filename = args.input
    if filename.endswith('sam'):
        sam_handle = pysam.Samfile(filename)
        bam_name = os.path.splitext(filename)[0]+".bam"
        bam_handle = pysam.Samfile( bam_name, "wb", template = sam_handle )
        for s in sam_handle:
            bam_handle.write(s)
        filename = bam_name

    # Massage paths and make outputfiles
    pysam.sort(filename, os.path.splitext(filename)[0]+".sorted")
    pysam.index(os.path.splitext(filename)[0]+".sorted.bam")
    filename = os.path.splitext(filename)[0]+".sorted.bam"
    path, filename = os.path.split(filename)
    name, ext = os.path.splitext(filename)
    new_name = name + '.wRG.' + 'bam'
    outfile_name =  os.path.join(path,new_name)
    outfile = pysam.Samfile( outfile_name, 'wb', header = newHeader )
    # Step 2: Process Samfile adding Read Group to Each Read
    samfile = pysam.Samfile(os.path.join(path, filename), 'rb')
    samfile.fetch()
    for count, read in enumerate(samfile.fetch()):
        name = read.qname
        read_group = readMap[name]
        if read.tags is None:
            read.tags = [("RG", read_group)]
        else:
            read.tags = read.tags + [('RG',read_group)]
        outfile.write(read)
    outfile.close()

    # Step 3: Make index of read group enabled samfile
    #pysam.index(outfile_name)
    #sys.stdout.write(".")
    #sys.stdout.flush()
    return

#def parseFileName(files):

#    return info

def add_RGs_2_BAM(args):
    print 'Making RG headers.'
        #if filename.endswith('sam') or filename.endswith('bam'):
            #sam_info = parseFileName(filename)
    read_map = {}
    new_RG_header = addRG2Header(read_map, args.file, args)
    add_RGs_2_BAMs_runner(new_RG_header, read_map, args)
    return

def main():
    args = get_args()
    add_RGs_2_BAM(args)

if __name__ == '__main__':
    main()
