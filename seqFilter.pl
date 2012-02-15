#!/usr/bin/env perl
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

use warnings;
use strict;

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


#globals
my (@gc, @seq_length, @coverage);

# genes {contig_name => [length, gc]}
my %genes;
my ($total_base, $seq_count, $total_gc, $seq_lt, $cov_lt, $gc_lt, $n_stat) = 0;
my $len_range = 0;
my $cov_range = 0;
my $gc_range = 1;
my $count = 0;
my ($outfh, $infh);
generate_parse_params();
if(defined $options->{'o'}) {
    open($outfh, '>', $options->{'o'}) or die $!;
} else {
    open($outfh, '>', \*STDOUT);
}

if(defined $options->{'i'}) {
    open($infh, '<', $options->{'i'}) or die $!;
} else {
    open($infh, '<', \*STDIN);
}

# set up the values needed for reading the fastx file
my @aux = undef;
my ($name, $seq, $qual);
while (($name, $seq, $qual) = readfq($infh, \@aux)) {
    $seq_counter++;
    my $length =  length($seq);
    $total_base += $length; 
    my $gc = calcgc($seq,$name);
    $total_gc += $gc;
    $genes{$name} = [$seq,$length,$gc];

    if(defined $options->{'l'}) {
        unless(filter_length($length)) {
            next;
        }
    }

    if(defined $options->{'g'}) {
        unless(filter_gc($gc)) {
            next;
        }
    }
    
    if(defined $options->{'c'}) {
        
      unless(filter_cov($cov)) {
            next;
        }
    }
    #calculate the specified n-statistic and print all sequences that are greater then
    #or equal to that value
    if (exists $options->{"n"}) {
        $genes{$name} = [$length, $gc, $seq];
    } else {
        print_out( $name, $seq, $length, $gc, $outfh );
    }
    $seq_count++;
}

if(exists $options->{'n'}) {
    $n_stat = $total_base / ($seq_count / (100 / $options->{"n"}));
    foreach my $genename (sort { $genes{$b}->[0] <=> $genes{$a}->[0] } keys %genes) {
        if ($genes{$genename}->[0] >= $n_stat) {
            print_out($genename,$genes{$genename}->[2],$genes{$genename}->[0],$genes{$genename}->[1], $outfh );
        }
    }
}
# print an aggregate report of all the sequences
if( exists ($options->{"a"})) {
    if ($seq_count > 1) {
        if (exists $options->{'n'}) {
            printf "\n\nAverage GC content is %.4f out of %d bases with an average length of %d and an n%s of %d\n\n", $total_gc / $total_base, $total_base, $total_base / $seq_count, $options->{'n'}, $n_stat;
        } else {
            my $n50 = $total_base / ($seq_count / 2);
            printf "\n\nAverage GC content is %.4f out of %d bases with an average length of %d and an n50 of %d\n\n", $total_gc / $total_base, $total_base, $total_base / $seq_count, $n50;
        }
		} else {
        printf "\n\nAverage GC content is %.4f out of %d bases\n\n", $total_gc / $total_base, $total_base;
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
sub checkParams 
{
    my %options;

    # Add any other command line options, and the code to handle them
    getopts( "i:g:o:hral:n:",\%options );

    # if no arguments supplied print the usage and exit
    #
   pod2usage if (0 == (keys (%options) ));

    # If the -h option is set, print the usage and exit
    #
    pod2usage if ($options{'h'});
        
    unless ($options{'g'} || $options{'n'} || $options{'l'} )
    {
        pod2usage('-msg' => "Please specify one of  -g -n -l ");
    }

    return \%options;
}

sub calcgc {
    my ($seq, $length) = @_;
    my $count = 0;
    $count++ while $seq =~ /[CG]/gi;
    return $count / $length;
}

sub generate_parse_params {
    if (exists ($options->{"g"})) {
          if ($options->{"g"} =~ /:/) {
              $gc_range = 1;
              @gc = split (/:/, $options->{"g"});
          } elsif ($options->{"g"} =~ /\,(0.\d{1,})/) {
          # if the user wants anything less than a specified gc percent
              $gc_lt = 1;
              push(@gc, $1);
          } elsif ($options->{"g"} =~ /(0.\d{1,})\,/) {
              # if the user wants anything greater than a specified gc percent
              $gc_lt = 0;
              push(@gc, $1);
          } else {
              pod2usage (-msg=> "\nERROR: input parameter error for -g\n");
          }
    }
    if (exists ($options->{"l"})) {
        if ($options->{"l"} =~ /:/) {
            $len_range = 1;
            @seq_length = split (/:/, $options->{"l"});
            } elsif ($options->{"l"} =~ /\,(\d+)/) {
                # if the user wants anything less than a specified length
                $seq_lt = 1;
                push(@seq_length, $1);
            } elsif ($options->{"l"} =~ /(\d+)\,/) {
                # if the user wants anything greater than a specified length
                $seq_lt = 0;
                push(@seq_length, $1);
            } else {
                pod2usage (-msg=> "\nERROR: input parameter error for -l\n");
            }
      }
}

sub filter_length {
    my ($length) = @_;
    if($len_range == 1) {
        if($length <= $seq_length[1] && $length >= $seq_length[0]) {
            return 1;
        } else {
            return 0;
        }
    } elsif ($seq_lt == 0) {
        if ($length >= $seq_length[0]) {
            return 1;
        } else {
            return 0;
        }
    } elsif ($seq_lt == 1) {
        if ($length <= $seq_length[0]) {
            return 1;
        } else {
            return 0;
        }
    } else {
      warn "wierd stuff goning on in length filter\n";
      return 0;
    }
}

sub filter_gc {
    my ($gc) = @_;
    if($gc_range == 1) {
        if($gc <= $gc[1] && $gc >= $gc[0]) {
            return 1;
        } else {
            return 0;
        }
    } elsif ($seq_lt == 0) {
        if ($gc >= $gc[0]) {
            return 1;
        } else {
            return 0;
        }
    } elsif ($seq_lt == 1) {
        if ($gc <= $gc[0]) {
            return 1;
        } else {
            return 0;
        }
    } else {
      warn "wierd stuff goning on in gc filter\n";
      return 0;
    }
}

sub print_out {
    my($name, $seq, $length, $gc, $outfh) = @_;
    if(exists $options->{'r'}) {
        print $outfh "$name\t$length\t$gc\n";
    } else {
        print $outfh ">$name\n$seq\n";
    }
}

sub printAtEnd 
{
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2010, 2011 Connor Skennerton
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME
 
 seqFilter.pl
 
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
 
 parses a fasta file based on length, gc, n-statistic
 Can also produce a report of sequences containing
 just the name, length, gc content for all sequences. 
 the aggregate option is provied for an at-a-glance look that the 
 overall sequence statistics of the file.
 
=head1 SYNOPSIS
 
 seqstats2.1.pl -i INPUT [-l INT:INT or ,INT or INT,] [-g DECIMAL:DECIMAL or ,DECIMAL or DECIMAL,] [-n INT] [-r] [-a] [-help]
 
    [-i]                      the name of the input fasta file [default: STDIN]
    [-o]                      name of the output file [default: STDOUT]
    [-l]                      parse based on specified length, can be a range or greater or less than a value 
    [-g]                      parse the sequences based on gc, can be a range or greater or less than a value
    [-n]                      perform a scecified n-statistic calculation and print and parse
    [-a]                      print aggregate stats for all the sequences
    [-r]                      Report format: print only the header information without the sequence 
    [-h]                   Displays detailed usage information
 
=head1 OPTIONS
 
 -i input file        the name of the fasta file to be parsed
 
 -l length            parse sequences based on the length value, specified as an integer.
                      there are two forms of this option. first is range where two values
                      are separated by a colon.  this will be interperated as the lower and
                      upper limits for the sequences; any sequences that fall between these 
                      two values (inclusive) will be printed
                      
                      in the second form a single comma is placed before or after the specified
                      value.  if the comma is before then it is interperated that the user would
                      like all values less than the specified length (inclusive); or greater than
                      if the comma is placed after the value.
 
 -g GC-content        parse sequences based on the GC value, specified as a decimal.
                      there are two forms of this option. first is range where two values
                      are separated by a colon.  this will be interperated as the lower and
                      upper limits for the sequences; any sequences that fall between these 
                      two values (inclusive) will be printed
                      
                      in the second form a single comma is placed before or after the specified
                      value.  if the comma is before then it is interperated that the user would
                      like all values less than the specified GC (inclusive); or greater than
                      if the comma is placed after the value.

 -c coverage          parse sequences based on the coverage value, specified as a integer.
                      there are two forms of this option. first is range where two values
                      are separated by a colon.  this will be interperated as the lower and
                      upper limits for the sequences; any sequences that fall between these 
                      two values (inclusive) will be printed
                      
                      in the second form a single comma is placed before or after the specified
                      value.  if the comma is before then it is interperated that the user would
                      like all values less than the specified coverage (inclusive); or greater than
                      if the comma is placed after the value.
                      
                      NOTE:  coverage information cannot be calculated from the sequence data.  
                      however popular short read assemblers (eg. Velvet, SaSSY) output coverage
                      information in the header of their contigs.  seqstats2.1 parses this header
                      for the coverage information using a regular expression.  COVERAGE CALCULATIONS
                      MAY FAIL IF YOU ARE PARSING CONTIGS FROM OTHER ASSEMBLERS THAT DO NOT USE THE 
                      SAME HEADER STYLE!
                      
                      Velvet header style = NODE_xx_len_xx_cov_xx
                      SaSSY header style = Contigxx_l_xx_c_xx
 
 -n n-statistic       the n-statistic is typically in the form of n50, which can be defined as the 
                      length of the shortest contig that is greater than 50% of the cumulative length
                      divided by the total number of contigs of the entire assembly.  seqstats2.1 can
                      parse based on any specified percentage given by the user.
                     
 -r report            modifies the print function so that the name, length, GC and coverage are 
                      printed in tab delimated format. default is to print the name and the 
                      sequence of the contig

 -a aggregate stats   prints a summary of the entire assembly including the average GC content, 
                      average length, total number of bases, and n-statistic.  if -n is specified
                      the n-statistic will be calculated using that variable else the default is
                      n50.
                      
                      

=head1 EXAMPLE USAGE

 seqstats2.1.pl -i seqs.fa -l 100:12000 >out_file
  
  <> open the file 'seqs.fa' and print to an ouput file all sequences that are between the 
  range of 100bp to 12kb
 
 seqstats2.1.pl -i seqs.fa -l 10000, -a -r
  
  <> print in report format, the largest 20 sequences larger than 10kb and print the aggregate 
  statistics for the entire assembly
  
 seqstats2.1.pl -i seqs.fa -g ,0.69 >out_file
  
  <> print to an output file all sequences that have a GC-content of less than 69%
  
 seqtats2.1.pl -i seqs.fa -n 70 -a -r
  
  <> print in report format all sequences over the n70 value and print the aggregate statistics
  for the entire assembly
  
  
=head1 TODO

  <> combine parsing based on coverage, GC and length, n-statistic
  <> implement a binning system for multiple output files based on parsing

=cut

