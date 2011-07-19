#!/usr/bin/perl
###############################################################################
#
#    this script extracts sequences from a multiple fasta file based on the 
#	 matches obtained from a blast file in m8 format
#
#    Copyright (C) 2010, 2011 Connor Skennerton
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

my $query = $options->{"input"};
my $database = $options->{"database"};
my $outfile = $options->{"output"};

my $database_format = 'fasta';
if (exists ($options->{"database_format"}))
{
  $database_format = $options->{"database_format"};
}
my $outfile_format = $database_format;
if (exists ($options->{'output_format'}))
{
  $outfile_format = $options->{'output_format'};
}

open(QUERY, $query) or die;

my %seqs;

while (my $line = <QUERY>) 
{
	chomp $line;
	if($options->{'list'})
	{
        list($line);
     }
     elsif ($options->{'blast'})
     {
        blast($line);
     }
     else
     {
        sam($line);
     }
}
close QUERY;

my $seq_in = Bio::SeqIO->new('-file' => $database,
                             	 '-format' => $database_format);

my $seq_out = Bio::SeqIO->new('-file' => '>'.$outfile, '-format' => $outfile_format);
while (my $seqobj = $seq_in->next_seq())
{
	if (exists $seqs{$seqobj->primary_id})
	{
     	unless($options->{'inverse'})
     	{
           	$seq_out->write($seqobj);
     	}
	}
	elsif ($options->{'inverse'})
	{
     	$seq_out->write($seqobj);
	}
}

printAtStart();

sub list{
  	my ($line) = shift;
	$seqs{$line} = 1;
}

sub blast{
  my ($line) = shift;
  my @columns = split(/\t/, $line);
  	if (exists $options->{'subject'})
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
  if ($line !~ /\*\t0\t0\t\*\t\*\t0\t0/) 
    {
      my @columns = split(/\t/, $line);
      $seqs{$columns[0]} = 1;
    }
}


sub checkParams 
{
    my @standard_options = ( "i|input:s", "if|input_format:s", "d|database:s", "s|subject:s", "o|output:s", "of|output_format:s", "df|database_format:s", "l|list:+", "b|blast:+", "S|sam:+", "h|help:+", "v|inverse:+" );
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
   pod2usage if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    pod2usage if ($options{'help'});
    
    pod2usage('-msg' => "please select one of the list, blast or sam options\n") unless ($options{'sam'} && $options{'blast'} && $options{'list'});
    
    pod2usage('-msg' => "The subject flag can only be specified with the blast flag\n") if ($options->{'subject'} && !$options->{'blast'});
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

__DATA__

=head1 NAME
	
 contig_extractor

=head1 COPYRIGHT
 
 copyright (C) 2010, 2011 Connor Skennerton
 
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
   
   used for extracting whole contigs from a multiple fasta file that contain 
   significant matches to reads/sequences/contigs from an m8 or m9 blast output file

=head1 SYNOPSIS
 
 contig_extractor -i|input FILE  -d|database SEQUENCE_FILE -o|output FILE {-l|list || -b|blast || -S|sam }  
                  [-s|subject] [-if|input_format FORMAT] [-of|output_format FORMAT] [-df|database_format FORMAT] 
                  [-h|help] [-v|inverse] 

      [-help]           Displays basic usage information
      [-sf]             The format of the file containing the contigs,the default is fasta    						
      -s                Name of the subject file containing the contigs
      -i                Name of the m8 blast file containing the matches to the contigs
      -o                Name of the output file
      [-h]              Use the hit (second column of blast table) to populate the list [default: use query]
      [-l]              Use a list of identifiers (one per line) to populate the list
      


=cut

