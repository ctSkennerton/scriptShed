#!/usr/bin/env perl
###############################################################################
#
#    geneiousTableToSequin.pl
#    Created by Connor Skennerton on 1/02/11.
#
#    Copyright 2011 Connor Skennerton. All rights reserved.
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

open (IN, $options->{'i'}) or die;
my $gene_count = 1;
print ">Feature ",$options->{'n'}," Table1\n";
while (my $line = <IN>)
	{
    chomp $line;

    my @columns = split(/,/, $line);
    if ($columns[0] =~ /CDS/)
    	{
        #print "made it to here\n";
        CDS_parsing(\@columns);
    	}
    elsif ($columns[0] =~ /term/)
    	{
        terminator_parsing(\@columns);
        }
    elsif ($columns[0] =~ /prom/)
    	{
        promotor_parsing(\@columns);
        }
    elsif ($columns[0] =~ /tRNA/)
    	{
        tRNA_parsing(\@columns);
        }
    else
        {
        misc_feature(\@columns);
        }
    }
    
     
       
printAtEnd();
close IN;
#close OUT;
exit;

sub CDS_parsing
{
#print "inside CDS sub\n";
my ($columns_ref) = @_;
my @columns = @{$columns_ref};
    my $start;
    my $end;
    my $feature_key;
    my $qualifier_key;
    my $qualifier_value;

if ($columns[-1] =~ /rev/)
    	{
        $end = $columns[2];
        $start = $columns[3];
        }
    else
    	{
        $end = $columns[3];
        $start = $columns[2];
        }
    $feature_key = $columns[0];
	$qualifier_key = 'gene';
    if ($columns[1] =~ /gp(\d+)/)
        {
        $qualifier_value = $options->{'n'}.'_'.$1;
        }
    elsif ($columns[1] =~ /Node/i)
    {
    	$qualifier_value = $options->{'n'}.'_'.$gene_count;
    }
    elsif ($columns[1] =~ /\d{1,2}/)
    	{
    	$qualifier_value = $options->{'n'}.'_'.$columns[1];
    	}
    else
    	{
    	$qualifier_value = $columns[1];
        }
    print "$start\t$end\tCDS\n\t\t\tproduct\t$qualifier_value\n$start\t$end\tgene\n\t\t\tgene\tgp";print"$gene_count\n";
    $gene_count++;

}

sub terminator_parsing
{
my ($columns_ref) = @_;
my @columns = @{$columns_ref};
    my $start;
    my $end;
    my $feature_key;
    my $qualifier_key;
    my $qualifier_value;

if ($columns[-1] =~ /rev/)
    	{
        $end = $columns[2];
        $start = $columns[3];
        }
    else
    	{
        $end = $columns[3];
        $start = $columns[2];
        }
    $feature_key = $columns[0];
    $qualifier_value = $columns[1];
    
    print "$start\t$end\tterminator\n\t\t\tgene\t$qualifier_value\n";

}

sub promotor_parsing
{
my ($columns_ref) = @_;
my @columns = @{$columns_ref};
    my $start;
    my $end;
    my $feature_key;
    my $qualifier_key;
    my $qualifier_value;

if ($columns[-1] =~ /rev/)
    	{
        $end = $columns[2];
        $start = $columns[3];
        }
    else
    	{
        $end = $columns[3];
        $start = $columns[2];
        }
    $feature_key = $columns[0];
    $qualifier_value = $columns[1];
    
    print "$start\t$end\tpromoter\n\t\t\tgene\t$qualifier_value\n";

}

sub tRNA_parsing
{
my ($columns_ref) = @_;
my @columns = @{$columns_ref};
    my $start;
    my $end;
    my $feature_key;
    my $qualifier_key;
    my $qualifier_value;

if ($columns[-1] =~ /rev/)
    	{
        $end = $columns[2];
        $start = $columns[3];
        }
    else
    	{
        $end = $columns[3];
        $start = $columns[2];
        }
    $feature_key = $columns[0];
    $qualifier_value = $columns[1];
    
    print "$start\t$end\ttRNA\n\t\t\tproduct\t$qualifier_value\n";

}

sub misc_feature
{
    my ($columns_ref) = @_;
    my @columns = @{$columns_ref};
    my $start;
    my $end;
    my $feature_key;
    my $qualifier_key;
    my $qualifier_value;
    
    if ($columns[-1] =~ /rev/)
    {
        $end = $columns[2];
        $start = $columns[3];
    }
    else
    {
        $end = $columns[3];
        $start = $columns[2];
    }
    $feature_key = $columns[0];
    $qualifier_value = $columns[1];
    
    print "$start\t$end\tmisc_feature\n\t\t\tnote\t$qualifier_value\n";
    
}


   sub checkParams 
{
    my @standard_options = ( "help+", "i:s", "n:s" );
    my %options;
    
    # Add any other command line options, and the code to handle them
    GetOptions( \%options, @standard_options );
    
    # if no arguments supplied print the usage and exit
    #
    pod2usage if (0 == (keys (%options) ));
    
    # If the -help option is set, print the usage and exit
    #
    pod2usage (-verbose=>2) if $options{'help'};
    
    return \%options;
    
}


sub printAtEnd {
print STDERR <<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2011 Connor Skennerton
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}


__DATA__

=head1 NAME

    geneiousTableToSequin.pl

=head1 COPYRIGHT

   copyright 2011 Connor Skennerton

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

  write detailed description here

=head1 SYNOPSIS

  geneiousTableToSequin.pl <options>
  
     -i                       the name of the input .csv file
     -n                       the name of the organism eg. EPV1
    [-help]                   Displays basic usage information

=cut



