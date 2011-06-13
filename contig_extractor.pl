#!/usr/bin/perl
###############################################################################
#
#    this script extracts sequences from a multiple fasta file based on the 
#	 matches obtained from a blast file in m8 format
#
#    Copyright (C) 2010 Connor Skennerton
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
use Getopt::Long;
use Bio::SeqIO;
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

my $query = $options->{"i"};
my $subject = $options->{"s"};
my $outfile = ">".$options->{"o"};
my $sf = 'fasta';
if (exists ($options->{"sf"}))
	{
	$sf = $options->{"sf"};
	}

open(QUERY, $query) or die;
open(OUT, $outfile) or die;

my %seqs;
my %contigs;

my $seq_in = Bio::SeqIO->new('-file' => $subject,
                             	 '-format' => $sf);
while (my $seqobj = $seq_in->next_seq())
		{
		$contigs{$seqobj->primary_id} = $seqobj->seq();
		}

while (my $line = <QUERY>) 
{
	chomp $line;
	my @columns = split(/\s+/, $line);

	if (exists $options->{'l'})
	{
		$seqs{$line} = 1;
		next;
	}

	elsif (exists $options->{'h'})
	{
		$seqs{$columns[1]} = $columns[0];
	}
	else
	{
		$seqs{$columns[0]} = $columns[1];
	}
	
}
close QUERY;
foreach my $match (sort keys %seqs)
	{                                                                                 
    print OUT ">$match\n";                                                                                                         
    print OUT "$contigs{$match}\n";
	}	
		
close OUT;
	

printAtStart();
sub checkParams {
    my @standard_options = ( "help+", "i:s", "s:s", "o:s", "sf:s", "l|list:+", "h|hit:+" );
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
   if (0 == (keys (%options) ))
   		{
   		usage() and die;
		}
    # If the -help option is set, print the usage and exit
    #
    if ($options{'help'})
    	{
    	help() and die;
    	}

    return \%options;
}


sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2010 Connor Skennerton
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

sub usage {

print "contig_extractor -i QUERY_FILE -s SUBJECT_FILE -o FILE_NAME [-sf] FORMAT [-help] \n";
}
sub help{
print "

   $0


   copyright (C) 2010 Connor Skennerton

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

DESCRIPTION
   
   used for extracting whole contigs from a multiple fasta file that contain 
   significant matches to reads/sequences/contigs from an m8 or m9 blast output file


SYNOPSIS\n\n";

   usage();
print"
      [-help]           Displays basic usage information
      [-sf]				the format of the file containing the contigs,the default is fasta    						
      -s				name of the subject file containing the contigs
      -q				name of the m8 blast file containing the matches to the contigs
      -o				name of the output file
      
";
}

exit;

