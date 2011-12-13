#!/usr/bin/env perl
###############################################################################
#
#    genbank2mfasta.pl -- extract all nucleotide cds to a multiple fasta file
#    Copyright (C) 2010 Michael Imelfort
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

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();
my $options = checkParams();

# get outptut file name
my $out_filename = $options->{'in'}."-features";
if(exists ($options->{"out"})) { $out_filename = $options->{"out"}; }

my $target_feature = $options->{"f"};
#opent the files
open (OUT_FILE, ">".$out_filename) or die $!;
my $seqio_object = Bio::SeqIO->new(-file => $options->{"in"});
my $seq_object = $seqio_object->next_seq;
 
for my $feat_object ($seq_object->get_SeqFeatures) 
{
    if ($feat_object->primary_tag eq "$target_feature") 
    {
        print OUT_FILE ">", $target_feature, "\t",$feat_object->start, "_", $feat_object->end,"\t";
        for my $tag ($feat_object->get_all_tags) 
        {  
            if( $tag =~ /note/)
            { 
                for my $value ($feat_object->get_tag_values($tag))
                {
                    print OUT_FILE "'$value'";
                    last;
                }
             }
        }
        print OUT_FILE "\n", $feat_object->spliced_seq->seq,"\n";
    }
}

# close the files
close OUT_FILE;

sub checkParams {
    my @standard_options = ( "help+", "in:s", "out:s", "f:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    return \%options;
}


sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2010 Michael Imelfort, Connor Skennerton
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    genbank2mfasta.pl

=head1 COPYRIGHT

   copyright (C) 2010 Michael Imelfort, Connor Skennerton

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

   Extract features from a genbank file to a multiple fasta file - uses bioperl

=head1 SYNOPSIS

   genbank2mfasta.pl -in FILENAME -f FEATURE [-out FILENAME] [-help]
     
      -in FILENAME       Genbank file to be parsed
      -f FEATURE         feature type to be extracted
      [-out FILENAME]    Output file name [default: infile-parsed.fasta]
      [-help]            Displays basic usage information
         
=cut

