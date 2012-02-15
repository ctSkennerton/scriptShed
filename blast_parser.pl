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
my ($input, $output);

if(defined $options->{'i'}) {
    open($input, '<', $options->{'i'}) or die $!;
} else {
    $input = \*STDIN;
}
if(exists $options->{"o"}) { 
	  open($output, '>', $options->{"o"}) or die $!; 
} else {
    $output = \*STDOUT;
}



my $in = new Bio::SearchIO(-format => $options->{'f'} -fh => $input);

while( my $result = $in->next_result ) 
{
    while( my $hit = $result->next_hit ) 
    {
        while( my $hsp = $hit->next_hsp ) 
        {
            my $mismatchcount = $hsp->length('total') - ($hsp->num_conserved + $hsp->gaps('total'));
            if(defined $options->{'l'}) {
                if( $hsp->length('total') < $options->{'l'}) {
                next;
                }
            }
            if(defined $options->{'p'}) {
                if ( $hsp->percent_identity < $options->{'p'}) {
                    next;
                }
            }
            if(defined $options->{'e'}) {
                if ( $hsp->significance > $options->{'e'}) {
                    next;
                }
            }
            if(defined $options->{'b'}) {	
                 if( $hsp->bits < $options->{'b'}) {
                    next;
                }
            }
            if(defined $options->{'m'}) {
                if ( $mismatchcount > $options->{'m'}) {
                    next;
                }
            }
            print $output join("\t", ( 	$result->query_name, $hit->name, sprintf("%.2f",$hsp->percent_identity), $hsp->length('total'),
            $mismatchcount, $hsp->gaps('total'), 
            $hsp->query->strand < 0 ?  ( $hsp->query->end, $hsp->query->start ) : ( $hsp->query->start, $hsp->query->end ),
            $hsp->hit->strand < 0 ?  ( $hsp->hit->end, $hsp->hit->start ) : ( $hsp->hit->start, $hsp->hit->end ),
            $hsp->evalue, $hsp->bits)), $hit->description,"\n";
        }  
    }
}
exit;

sub checkParams {
    my @standard_options = ( "h+", "i:s", "o:s", "l:s", "p:s", "e:s", "f:s", "b:s", "m:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    #exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'h'};
	
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

   copyright (C) 2010 2011 Connor Skennerton

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

   blast_parser.pl [-i FILENAME] [-o FILENAME] [-h] [-q INTEGER] [-p INTEGER] [-e INTEGER] [-b INTEGER] [-m INTEGER] [-f FORMAT]
   				   
     
      -i FILENAME               Blast file to be parsed [default: STDIN]
      [-o FILENAME]             Output file name [default: STDOUT]
      [-h]                      Displays basic usage information
      [-l INTEGER]              The length of the match that the HSP must be higher than
      [-p INTEGER]              The required percent identity that the HSP must be higher than
      [-e INTEGER]              The required e-value that the HSP must be lower than
      [-b INTEGER]              The required bits score that the HSP must be higher than
      [-f FORMAT]               The file format to be used, default is 'blasttable'
      [-m INTEGER]              The number of mismatches that the HSP must be lower than
=cut

