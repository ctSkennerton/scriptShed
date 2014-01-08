#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;


BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

my $options = checkParams();

# get outptut file name
my $out_file; # = $options->{'in'}.".outf";
my $in_file;
if(exists ($options->{"out"})) {  open($out_file, '>', $options->{"out"}) || die $!; }
else{ $out_file = \*STDOUT;}

if(exists ($options->{"in"})) {  open($in_file, '<', $options->{"in"}) || die $!; }
else{ $in_file = \*STDIN;}



# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-fh' => $in_file,
                             '-format' => $options->{"inf"});
my $seq_out = Bio::SeqIO->new('-fh' => $out_file,
                             '-format' => $options->{"outf"});

# write each entry in the input file to the output file
while (my $seq_in = $seq_in->next_seq) {
	$seq_out->write_seq($seq_in);
}

sub checkParams {
    my @standard_options = ( "help+", "in:s", "out:s", "inf:s", "outf:s" );
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
 Copyright (C) 2010 Connor Skennerton
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    fileconverter

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

   this script reads in any given sequence format and converts it to another
   format

=head1 SYNOPSIS

   script_name [-help] [-in FILENAME] [-inf FORMAT] [-out FILENAME] [-outf FORMAT]

      [-help]            Displays basic usage information
      [-in]				specify the input filename
      [-out]			specify the output filename - if none is specified
      					the input file name will be appended with the output file format
      [-inf]			the file format of the input file
      [-outf]			the desired output format
         
=cut
