#!/usr/bin/perl
###############################################################################
#
#    this script extracts sequences from a multiple fasta file based on the 
#	 matches obtained from a blast file in m8 format
#
#    Copyright (C) 2010, 2011, 2012 Connor Skennerton
#
#    Support added for the Manotator by Mike Feb 2012
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Std;
use Pod::Usage;
#CPAN modules

#locally-written modules

BEGIN 
	{
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
	}

# get input params and print copyright

my $options = checkParams();


my $query = \*STDIN;
if (defined $options->{'i'}) {
    open($query, $options->{'i'}) or die;
}
my $database = $options->{'d'};
my $outfile = \*STDOUT; 
if (defined $options->{'o'}) {
    open($outfile, ">", $options->{'o'}) or die;
}

my %seqs;
#printAtStart();

if($options->{'f'})
{
    my @aux = undef;
    my ($name, $seq, $qual);
    while (($name, $seq, $qual) = readfq($query, \@aux)) 
    {
        $seqs{$name} = 1;
    }
}
elsif (! defined $options->{'c'})
{
    while (my $line = <$query>) 
    {
        chomp $line;
        if($options->{'l'})
        {
            list($line);
        }
        elsif ($options->{'b'})
        {
            blast($line);
        }
        elsif ($options->{'S'})
        {
            sam($line);
        }
        else #($options->{'U'})
        {
            mannotator("UniRef90_".$line);
        }
    }
} else {
    my @s = split(/,/, $options->{'c'});
    foreach my $e (@s) {
        $seqs{$e} = 1;
    }
}
close $query;

my @aux = undef;
my ($name, $seq, $qual);
open(DB,$database) or die;
while (($name, $seq, $qual) = readfq(\*DB, \@aux)) 
{
	if (exists $seqs{$name})
	{
     	unless($options->{'v'})
     	{
           	print_seq(\$name,\$seq,\$qual, $outfile);
     	}
	}
	elsif ($options->{'v'})
	{
           	print_seq(\$name,\$seq,\$qual, $outfile);
	}
}


sub print_seq{
    my ($name_ref, $seq_ref, $qual_ref, $fh) = @_;
    if (defined $$qual_ref)
    {
        # fastq file
        print $fh "@".$$name_ref."\n".$$seq_ref."\n+".$$name_ref."\n".$$qual_ref."\n";
    }
    else
    {
        print $fh ">".$$name_ref."\n".$$seq_ref."\n";
    }
}

sub readfq {
	my ($fh, $aux) = @_;
	@$aux = [undef, 0] if (!defined(@$aux));
	return if ($aux->[1]);
	if (!defined($aux->[0])) {
		while (<$fh>) {
			chomp;
			if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
				$aux->[0] = $_;
				last;
			}
		}
		if (!defined($aux->[0])) {
			$aux->[1] = 1;
			return;
		}
	}
	my $name = /^.(\S+)/? $1 : '';
	my $seq = '';
	my $c;
	$aux->[0] = undef;
	while (<$fh>) {
		chomp;
		$c = substr($_, 0, 1);
		last if ($c eq '>' || $c eq '@' || $c eq '+');
		$seq .= $_;
	}
	$aux->[0] = $_;
	$aux->[1] = 1 if (!defined($aux->[0]));
	return ($name, $seq) if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
		chomp;
		$qual .= $_;
		if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
			return ($name, $seq, $qual);
		}
	}
	$aux->[1] = 1;
	return ($name, $seq);
}

sub list{
  	my ($line) = shift;
	$seqs{$line} = 1;
}

sub blast{
  my ($line) = shift;
  my @columns = split(/\t/, $line);
  	if (exists $options->{'s'})
	{
		$seqs{$columns[1]} = $columns[0];
	}
	else
	{
		$seqs{$columns[0]} = $columns[1];
	}
}

sub sam{ 
    my ($line) = shift;
    my @c = split(/\t/,$line);
    # test whether the third bit is set - query unmapped
    unless($c[1] & 4)
    #if ($line !~ /\*\t0\t0\t\*\t\*\t0\t0/) 
    {
        $seqs{$c[0]} = 1;
    }
}

sub mannotator{
    my ($line) = shift;
    my @columns = split /\^/, $line;
    $seqs{$columns[0]} = 1;
}

sub checkParams 
{
    my %options;

    # Add any other command line options, and the code to handle them
    getopts( "i:d:so:c:lbShvfU",\%options );

    # if no arguments supplied print the usage and exit
    #
   pod2usage if (0 == (keys (%options) ));

    # If the -h option is set, print the usage and exit
    #
    pod2usage if ($options{'h'});
    unless ($options{'c'}) {
        unless ($options{'S'} || $options{'U'} || $options{'b'} || $options{'l'} || $options{'f'} )
        {
            pod2usage('-msg' => "Please specify one of  -S -b -l -f -U");
        }
    }
    unless (($options{'i'} || $options{'c'}) && $options{'d'})
    {
        pod2usage('-msg' => "You must specify -d and either -i or -c");
    }
    
    if (defined $options{'s'} && !(defined $options{'b'}))
    {
        pod2usage('-msg' => "The subject flag can only be specified with the blast flag\n");
    }

    
    return \%options;
}


sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2010, 2011, 2012 Connor Skennerton
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME
	
 contig_extractor

=head1 COPYRIGHT

 copyright (C) 2010, 2011, 2012 Connor Skennerton

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   Used for extracting whole contigs / sequences from a multiple fasta file that contain 
   significant matches to reads/sequences/contigs from a variety of list formats

=head1 SYNOPSIS

 contig_extractor { [-c CONTIG_NAMES] | [-i FILE] -l|b|S|f|U } -d SEQUENCE_FILE [-o FILE] [-s] [-h] [-v] 

      [-h]              Displays basic usage information
      -d                Name of the subject file containing the contigs
      [-i]              Name of the file containing the matches to the contigs. Incompatible with -c
      [-o]              Name of the output file
      [-l]              Use a list of identifiers (one per line) to populate the list
      [-b]              Use a m8 blast file as the input
      [-s]              Generate headers based on the subject of a blast file.  default is to use the query
      [-S]              Use a sam file as the input
      [-f]              Use a fasta/fastq file as the input
      [-v]              Invert the match
      [-c]              list the names of sequences to extract as a comma separated list. Incompatible with -i
      [-U]              Input file is a mannotator annotation mappings file (default is the uniref file)


=cut

