#!/usr/bin/perl
###############################################################################
#
#    parses a blast output file based on user defined characteristics
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
use Bio::SearchIO; 
#CPAN modules

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright

my $options = checkParams();

my $out= $options->{'in'}."-blast_parser";
if(exists ($options->{"out"})) 
	{ 
	$out = $options->{"out"}; 
	}
my $p = 1; 
if(exists ($options->{"p"})) 
	{
	$p = $options->{"p"}; 
	}
my $l = 1; 
if(exists ($options->{"l"})) 
	{
	$l = $options->{"l"}; 
	}
my $e = 10; 
if(exists ($options->{"e"}))
	{
	$e = $options->{"e"};
	}
my $b = 1; 
if(exists ($options->{"b"}))
	{
	$b = $options->{"b"};
	}
my $f = 'blasttable';
if (exists($options->{"f"}))
	{
	$f = $options->{"f"};
	}
my $m = 10000000;
if (exists($options->{"m"}))
	{
	$m = $options->{"m"};
	}
my $in = new Bio::SearchIO(-format => $f, 
                           -file   => ($options->{"in"}));
open(OUT, ">".$out) or die;	

while( my $result = $in->next_result ) 
{
	while( my $hit = $result->next_hit ) 
    {
		while( my $hsp = $hit->next_hsp ) 
        {
			my $mismatchcount = $hsp->length('total') - 
			($hsp->num_conserved + $hsp->gaps('total'));
			if( $hsp->length('total') >= $l) 
            {
        		if ( $hsp->percent_identity >= $p) 
                {
         			if ( $hsp->significance <= $e) 
                    {
         				if ( $hsp->bits >= $b)
         				{	
                            if ( $mismatchcount <= $m)
                            {
                                print OUT join("\t", ( 	$result->query_name,
                                $hit->name,
                                sprintf("%.2f",$hsp->percent_identity),
                                $hsp->length('total'),
                                $mismatchcount,
                                $hsp->gaps('total'),
                                $hsp->query->strand < 0 ?
                                ( $hsp->query->end,
                                $hsp->query->start ) :
                                ( $hsp->query->start,
                                $hsp->query->end ),
                                $hsp->hit->strand < 0 ?
                                ( $hsp->hit->end,
                                $hsp->hit->start ) :
                                ( $hsp->hit->start,
                                $hsp->hit->end ),
                                $hsp->evalue,
                                $hsp->bits)), $hit->description,"\n";
                            }
         				}
                    }
                } 
            }
        }  
    }
}
close OUT;
printAtStart();
exit;

sub checkParams {
    my @standard_options = ( "help+", "in:s", "out:s", "l:s", "p:s", "e:s", "f:s", "b:s", "m:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    #exec("pod2usage $0") if $options{'help'};
	
	# if there is no input file, print the usage and exit
	#exec("pod2usage $0") if (! defined ($options->{"in"}));
	

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

    blast_parser.pl

=head1 COPYRIGHT

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

=head1 DESCRIPTION

   parses a blast file on using various characteristics - uses bioperl

=head1 SYNOPSIS

   blast_parser.pl -in FILENAME [-out FILENAME] [-help] [-q INTEGER] [-p INTEGER] [-e INTEGER] [-b INTEGER] [-m INTEGER] [-f FORMAT]
   				   
     
      -in FILENAME                blast file to be parsed
      [-out FILENAME]           Output file name [default: infile-blast_parsed]
      [-help]                        Displays basic usage information
      [-l INTEGER]                the length of the match that the HSP must be higher than or 
                                        equal to, default is 1
      [-p INTEGER]              the required percent identity that the HSP must 
                                       be higher than or equal to, default is 1
      [-e INTEGER]              the required e-value that the HSP must be lower than or equal to,
                                       default is 10
      [-b INTEGER]              the required bits score that the HSP must be higher than or equal to,
                                       the default is 1  
      [-f FORMAT]               the file format to be used, default is 'blasttable'		 
      [-m INTEGER]            the number of mismatches that the HSP must be lower than or
                                      equal to, the default is 10000000
=cut

