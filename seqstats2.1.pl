#!/usr/bin/env perl
###############################################################################
#
#    seqstats2.1.pl - parses sequences based on various parameters
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
use Pod::Usage;
use Bio::SeqIO;
use Bio::Tools::Seqstats;
use Data::Dumper;
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

#globals
my (@gc, @seq_length, @coverage);
my %genes;
my ($total_base, $seq_count, $total_gc, $seq_lt, $cov_lt, $gc_lt, $n_stat) = 0;
my $len_range = 0;
my $cov_range = 0;
my $gc_range = 0;
my $max;
my $count = 0;
generate_parse_params();

my $seqin = new Bio::SeqIO(-format => 'fasta',
                           -file   => $options->{"i"});


#create the hash
while( my $seq = $seqin->next_seq() ) 
    {
    #if there is no sequence in the fasta file or it appears that the file
    #is amino acids skip the sequence
    next if( $seq->length == 0 );
    if( $seq->alphabet eq 'protein' ) 
        {
        print $seq->display_id, "\t";
        warn("does not work on amino acid sequences ...skipping this seq");
        next;
        }
    $seq_count++;
    $total_base += $seq->length;
    
    #//CTS// placeholder for now - extracts coverage info from a velvet or sassy header
    # in the future should be replaced with something more flexable
    my @header_columns = split(/_/, $seq->display_id());
    

    my $gc_calc = calcgc($seq->seq(), $seq->length());
	
    $genes{$seq->primary_id} = { "seq" => $seq->seq, "length" => $seq->length, "GC" => $gc_calc, "cov" => $header_columns[-1] }
	}

    
if (exists $options->{"m"})
    {
    $max = $options->{"m"};
    }
else
    {
    $max = $seq_count;
    }
#print Dumper(%genes);
# sort the hash based on length and print
if (exists $options->{"l"})
{
sort_length();
}
#sort the hash based on gc and print
if (exists $options->{"g"})
{
sort_gc();
}
#sort based on coverage and print
if (exists $options->{"c"})
{
sort_cov();
}

#calculate the specified n-statistic and print all sequences that are greater then
#or equal to that value

if (exists $options->{"n"})
    {    
    $n_stat = $total_base / ($seq_count / (100 / ($options->{"n"})));
    foreach my $genename (sort { $genes{$b}->{"length"} <=> $genes{$a}->{"length"} } keys %genes)
        {
        if ($genes{$genename}->{'length'} >= $n_stat)
            {
            print_out($genename);# ">", "$genename\n$genes{$genename}->{'seq'}\n";
            $count++;
            last if ($count == $max);
            }
        }
    }

# print an aggregate report of all the sequences
if( exists ($options->{"a"})) 
    {
    if ($seq_count > 1)
        {
        if (exists $options->{'n'})
        	{
        	printf "\n\nAverage GC content is %.4f out of %d bases with an average length of %d and an n%s of %d\n\n", $total_gc / $total_base, $total_base, $total_base / $seq_count, $options->{'n'}, $n_stat;
        	}
        else
        	{
        	my $n50 = $total_base / ($seq_count / 2);
        	printf "\n\nAverage GC content is %.4f out of %d bases with an average length of %d and an n50 of %d\n\n", $total_gc / $total_base, $total_base, $total_base / $seq_count, $n50;
            }
		}
    else
        {
        printf "\n\nAverage GC content is %.4f out of %d bases\n\n", $total_gc / $total_base, $total_base;
        }
    }
    

#close OUT;
printAtEnd();
exit;

sub calcgc {
    my ($seq, $length) = @_;
    my @seqarray = split('',$seq);
    my $count = 0;
    foreach my $base (@seqarray) 
    	{
		if ($base =~ /[G|C]/i)
        	{
            $count++ 
            }
    	}
    return $count / $length;
}

sub generate_parse_params
{
if (exists ($options->{"g"}))
	{
    if ($options->{"g"} =~ /:/) 
    	{
        $gc_range = 1;
        @gc = split (/:/, $options->{"g"});
        }
    # if the user wants anything less than a specified gc percent
    #//CTS// check that regex is correct
    elsif ($options->{"g"} =~ /\,(0.\d{1,})/)
    	{
        $gc_lt = 1;
        push(@gc, $1);
        }
    # if the user wants anything greater than a specified gc percent
    #//CTS// check that regex is correct    
    elsif ($options->{"g"} =~ /(0.\d{1,})\,/)
    	{
        $gc_lt = 0;
        push(@gc, $1);
        }
    else
    	{
        pod2usage (-msg=> "\nERROR: input parameter error for -g\n");
        }
	}
if (exists ($options->{"l"}))
	{
    if ($options->{"l"} =~ /:/) 
    	{
        $len_range = 1;
        @seq_length = split (/:/, $options->{"l"});
        }
    # if the user wants anything less than a specified length
    elsif ($options->{"l"} =~ /\,(\d+)/)
    	{
		$seq_lt = 1;
        push(@seq_length, $1);
        }
    # if the user wants anything greater than a specified length
    elsif ($options->{"l"} =~ /(\d+)\,/)
    	{
        $seq_lt = 0;
        push(@seq_length, $1);
        } 
    else
    	{
        pod2usage (-msg=> "\nERROR: input parameter error for -l\n");
        }
	}
if (exists ($options->{"c"}))
	{
    if ($options->{"c"} =~ /:/) 
    	{
        $cov_range = 1;
        @coverage = split (/:/, $options->{"c"});
        }
    # if the user wants anything less than a specified length
    elsif ($options->{"c"} =~ /\,(\d+)/)
    	{
        $cov_lt = 1;
        push(@coverage, $1);
        }
    # if the user wants anything greater than a specified length
    elsif ($options->{"c"} =~ /(\d+)\,/)
    	{
        $cov_lt = 0;
        
        push(@coverage, $1);

        } 
    else
    	{
        pod2usage (-msg=> "\nERROR: input parameter error for -c\n");
        }
	}
}

sub sort_length
{
foreach my $genename (sort { $genes{$b}->{"length"} <=> $genes{$a}->{"length"} } keys %genes)  
	{
        if ($len_range == 1)
    	{
            if (($genes{$genename}->{'length'} <= $seq_length[1]) and ($genes{$genename}->{'length'} >= $seq_length[0]))
        	{
                print_out($genename); #">", "$genename\n$genes{$genename}->{'seq'}\n";
				$count++;
                last if ($count == $max);            
            }
        }
        elsif ($seq_lt == 0) 
        {
            if ($genes{$genename}->{'length'} >= $seq_length[0])
            {
                print_out($genename); #">", "$genename\n$genes{$genename}->{'seq'}\n";
				$count++;
                last if ($count == $max);
            }
        }
        elsif ($seq_lt == 1)
        {
            if($genes{$genename}->{'length'} <= $seq_length[0])
            
            {
                print_out($genename); #">", "$genename\n$genes{$genename}->{'seq'}\n";
				$count++;
                last if ($count == $max);
            }
        }    	
    else
    	{
    	print STDERR "Oops... no sequences match your criteria\n"; exit;
    	}
	}
}

sub sort_gc
{
foreach my $genename (sort { $genes{$a}->{"GC"} <=> $genes{$b}->{"GC"} } keys %genes)  
	{
        if ($gc_range == 1)
    	{
            if (($genes{$genename}->{'GC'} <= $gc[1]) and ($genes{$genename}->{'GC'} >= $gc[0]))
        	{
                print_out($genename); #">", "$genename\n$genes{$genename}->{'seq'}\n";
				$count++;
                last if ($count == $max);

            }
        }
        elsif ($gc_lt == 0) 
        {
            if ($genes{$genename}->{'GC'} >= $gc[0])
            {
                print_out($genename); #">", "$genename\n$genes{$genename}->{'seq'}\n";
				$count++;
                last if ($count == $max);

            }
        }
        elsif ($gc_lt == 1)
        {
            if($genes{$genename}->{'GC'} <= $gc[0])
            
            {
                print_out($genename); 
				$count++;
                last if ($count == $max);

            }
        }    	
        else
    	{
    	print STDERR "Oops... no sequences match your criteria\n"; exit;
    	}
    }
}

sub sort_cov
{
foreach my $genename (sort { $genes{$b}->{"cov"} <=> $genes{$a}->{"cov"} } keys %genes)  
	{
        if ($cov_range ==1) #//CTS use of uninitialized value//
    	{
            if (($genes{$genename}->{'cov'} <= $coverage[1]) and ($genes{$genename}->{'cov'} >= $coverage[0]))
        	{
                print_out($genename); 
				$count++;
                last if ($count == $max);
            }
        }
        elsif ($cov_lt == 0) 
        {
            if ($genes{$genename}->{'cov'} >= $coverage[0]) #//CTS use of uninitialized value
            {
                print_out($genename); 
				$count++;
                last if ($count == $max);
            }
        }
        elsif ($cov_lt == 1)
        {
            if($genes{$genename}->{'cov'} <= $coverage[0])
            
            {
                print_out($genename); 
				$count++;
                last if ($count == $max);
            }
        }    	
    else
    	{
    	print STDERR "Oops... no sequences match your criteria\n"; exit;
    	}
	}
}

sub print_out
{
my ($genename) = @_;
if (exists $options->{'r'})
	{
    print "$genename\t$genes{$genename}->{'length'}\t$genes{$genename}->{'GC'}\t$genes{$genename}->{'cov'}\n";
    }
else
	{
    print ">", "$genename\n$genes{$genename}->{'seq'}\n"
    }
}
sub checkParams 
{
    my @standard_options = ( "help+", "i:s", "l:s", "g:s", "c:s", "m:s", "n:s", "r+", "a+" );
    my %options;
    
    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );
    
    # if no arguments supplied print the usage and exit
    #
    pod2usage if (0 == (keys (%options) ));
    
    # If the -help option is set, print the usage and exit
    #
    pod2usage (-verbose=>2) if $options{'help'};
    
    #if there is no input file print error msg and usage
    pod2usage (-msg=> "\nERROR: no input file provided\n") if (!exists ($options{"i"}));
    
    #if there are no options print error msg and usage
    if ((!exists $options{'g'}) and(!exists $options{'l'}) and (!exists $options{'c'}) and (!exists $options{'n'}))
    	{
        pod2usage(-msg=> "\nERROR: you haven't given me any options to parse on!\n\tplease specify either -l or -g or -c or -n \n");
        }
    
    return \%options;
    
}


sub printAtEnd 
{
print STDERR <<"EOF";
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
 
 seqstats2.1.pl
 
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
 
 parses a fasta file based on length, gc, n-statistic, coverage or number  
 of sequences to keep. Can also produce a report of sequences containing
 just the name, length, gc content and coverage for all sequences. The 
 report can be sorted based on  gc content, length or coverage if no option
 is provided then the sequences will be printed in a random order.
 the aggregate option is provied for an at-a-glance look that the 
 overall sequence statistics of the file.
 
=head1 SYNOPSIS
 
 seqstats2.1.pl -i INPUT [-l INTEGER:INTEGER or ,INTEGER or INTEGER,] [-g DECIMAL:DECIMAL or ,DECIMAL or DECIMAL,] 
                       [-c INTEGER:INTEGER or ,INTEGER or INTEGER,] [-n INTEGER] [-m INTEGER] [-r] [-a] [-help]
 
     -i                       the name of the input fasta file
    [-l]                      parse based on specified length, can be a range or greater or less than a value 
    [-g]                      parse the sequences based on gc, can be a range or greater or less than a value
    [-c]                      parse the sequences based on coverage, can be a range or greater or less than a value
    [-n]                      perform a scecified n-statistic calculation and print and parse
    [-m]                      maximum number of records to print out
    [-a]                      print aggregate stats for all the sequences
    [-r]                      Report format: print only the header information without the sequence 
    [-help]                   Displays detailed usage information
 
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
                     
 -m maximum           modifier of -c -g -l -n prints only the specified maximum number of sequences
 
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
 
 seqstats2.1.pl -i seqs.fa -l 10000, -a -r -m 20 
  
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



