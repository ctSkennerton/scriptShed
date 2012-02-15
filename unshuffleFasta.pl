#!/usr/bin/env perl

if (!@ARGV) {
    print "Usage: $0 infile.fa forward_reads.fa reverse_reaads.fa\n";
    print "\tforward_reads.fa / reverse_reads.fa : reads to be unshuffled\n";
    print "\tinfile.fa :file containing shuffled fasta reads\n";
    exit;
}

$filename = $ARGV[0];
$filenameOutA = $ARGV[1];
$filenameOutB = $ARGV[2];

die "Could not open $filename" unless (-e $filename);
#die "Could not open $filenameB" unless (-e $filenameB);

open FILEA, "> $filenameOutA";
open FILEB, "> $filenameOutB";

open INFILE, "< $filename";

my $lineA;

$lineA = <INFILE>;
#$lineB = <FILEB>;

while(defined $lineA) {
    print FILEA $lineA;
    $lineA = <INFILE>;
    while (defined $lineA && (substr($lineA, 0,1) != '>')) {
        print FILEA $lineA;
        $lineA = <INFILE>;
    }

    print FILEB $lineA;
    $lineA = <INFILE>;
    while (defined $lineA && (substr($lineA, 0,1) != '>')) {
        print FILEB $lineA;
        $lineA = <INFILE>;
    }
}
