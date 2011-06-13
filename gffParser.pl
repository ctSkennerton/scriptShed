#!/usr/bin/env perl

use warnings;
use strict;

my $gff = $ARGV[1];
my $query = $ARGV[0];

open (GFF, '<', $gff) or die $!;
open (QUERY, '<', $query) or die $!;

my %gff_hash;
my %query_hash;
while(<QUERY>)
{
	chomp $_;
	$query_hash{$_} = 1;
}
close QUERY;

while (<GFF>)
{
	my ($contig, $annotation) = split(/\t/, $_, 2);
	if (exists $query_hash{$contig})
	{
	print $_;
	}
}

close GFF;
exit;