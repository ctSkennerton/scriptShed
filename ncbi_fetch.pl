#!/usr/bin/env perl
#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long qw(:config no_ignore_case bundling);
use Carp;

#CPAN modules
use Bio::DB::EUtilities;
use Bio::SeqIO;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands 
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

my $global_options = checkParams();
my $format = overrideDefault('gb', 'outformat');
my $database = overrideDefault('protein', 'db');
my $file = overrideDefault('file.'.$format, 'outfile');

# set optional history queue
if(defined $global_options->{'query'}) {
    my $factory = Bio::DB::EUtilities->new(-eutil      => 'esearch',
                                           -email      => 'mymail@foo.bar',
                                           -db         => $database,
                                           -term       => $global_options->{'query'},
                                           -usehistory => 'y');
     
    my $count = $factory->get_count;
    # get history from queue
    my $hist  = $factory->next_History || die 'No history data returned';
    print "History returned\n";
    # note db carries over from above
    $factory->set_parameters(-eutil   => 'efetch',
                             -rettype => $format,
                             -history => $hist);
     
    my $retry = 0;
    my ($retmax, $retstart) = (500,0);
     
    open (my $out, '>', $file) || die "Can't open file:$!";
     
    RETRIEVE_SEQS:
    while ($retstart < $count) {
        $factory->set_parameters(-retmax   => $retmax,
                                 -retstart => $retstart);
        eval{
            $factory->get_Response(-cb => sub {my ($data) = @_; print $out $data} );
        };
        if ($@) {
            die "Server error: $@.  Try again later" if $retry == 5;
            print STDERR "Server error, redo #$retry\n";
            $retry++ && redo RETRIEVE_SEQS;
        }
        print "Retrieved $retstart\n";
        $retstart += $retmax;
    }
     
    close $out;
} else {

    my @ids;
    if(defined $global_options->{'patterns'}) {
        open(PAT, '<', $global_options->{'patterns'}) || die $!;
        while(<PAT>) {
            chomp;
            push @ids, $_;
        }
        close PAT;
    } else {
        @ids = @ARGV;
    }

    my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
        -db      => $database,
        -rettype => $format,
        -email   => 'mymail@foo.bar',
        -id      => \@ids);

    # dump HTTP::Response content to a file (not retained in memory)
    $factory->get_Response(-file => $file);
}

# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "db|d=s", "outfile|o=s", "outformat|F=s", "patterns|f=s", "query|q=s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    #exec("pod2usage $0") if (0 == (scalar (@ARGV) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    #if(!exists $options{''} ) { printParamError (""); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #  
    my ($error) = @_;  
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

    ncbi_fetch.pl

=head1 COPYRIGHT

   copyright (C) Connor Skennerton

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

   Reusable script for downloading information from NCBI.using bioperl's eutils.
   See http://www.ncbi.nlm.nih.gov/books/NBK25499/ for a good resource on eutils
   and http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.chapter4_table1/?report=objectonly
   for a look at the database names and output formats that can be used with efetch

=head1 SYNOPSIS

    ncbi_fetch.pl  [-h] [-d DATABASE] [-o FILE] [-F FORMAT] [-f FILE] [-q QUERY] [ID [ID ...]]

      [-h]               Displays basic usage information
      [-f FILE]          Specify a file that contains IDs, one per line to fetch
      [-d DATABASE]      Specify NCBI database to fetch from. Common options would
                         be 'protein' or 'nuccore'
      [-F FORMAT]        Specify the output format. Default: genbank (gb)
                         The available formats are based on the queried database,
                         however common options would be 'gb' for genbank or 'fasta'
                         for fasta.  Other handy formats are 'fasta_cds_na' which will
                         give you a multiple fasta file containing all the CDS regions
                         if say you specify a genome as the ID!
      [-o FILE]          Specify a file to print results to. Default: file.FORMAT
      [-q QUERY]         Specify a query used to search NCBI, matches to the query will
                         be downloaded. e.g. 'Viruses[Organism] AND "Complete Genome"[All Fields]'
=cut
