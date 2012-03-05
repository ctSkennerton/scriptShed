#!/usr/bin/env perl
###############################################################################
#
#    __Script__Name__
#    
#    <one line to give the program's name and a brief idea of what it does.>
#
#    Copyright (C) 2012 Michael Imelfort
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
my $global_options = checkParams();

my $in_fh = \*STDIN;
if(defined $global_options->{'in'}) {
    open $in_fh, '<', $global_options->{'in'} or die $!;
}

my $pair_fh = \*STDOUT;
if(defined $global_options->{'paired'}) {
    open $pair_fh, '>', $global_options->{'paired'} or die $!;
}

my $single_fh = \*STDERR;
if(defined $global_options->{'single'}) {
    open $single_fh, '>', $global_options->{'single'} or die $!;
}

my @aux = undef;
my ($name, $seq, $qual);
my %pairs_hash;
while (($name, $seq, $qual) = readfq($in_fh, \@aux)) {
    if($name =~ /(.*)\/(\d).*/) {
        push @{$pairs_hash{$1}}, [$name, $seq, $qual, $2];
    }
}
close $in_fh;
while(my($k,$v) = each %pairs_hash) {
    if(scalar @{$v} > 1) {
        foreach my $e (sort {$a->[3] <=> $b->[3]} @{$v}) {
            if(defined $e->[2]){
                print $pair_fh "\@$e->[0]\n$e->[1]\n\+$e->[0]\n$e->[2]\n";
            } else {
                print $pair_fh ">$e->[0]\n$e->[1]\n";
            }
        }
    } else {
        if(defined $v->[0]->[2]){
            print $single_fh "\@$v->[0]->[0]\n$v->[0]->[1]\n\+$v->[0]->[0]\n$v->[0]->[2]\n";
        } else {
            print $single_fh ">$v->[0]->[0]\n$v->[0]->[1]\n";
        }
    }

}
close $single_fh;
close $pair_fh;
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

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "paired|p:s","in|i:s", "single|s:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    #pod2usage() if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    pod2usage() if $options{'help'};

    # Compulsosy items
    #if(!exists $options{''} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2012 Connor Skennerton 
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or die "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or die "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
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

    __Script__Name__

=head1 COPYRIGHT

   copyright (C) 2012 Michael Imelfort

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

   Insert detailed description here

=head1 SYNOPSIS

    __Script__Name__  [-help|h] [-in|i] [-paired|p] [-single|s]

      [-help -h]                   Displays basic usage information
      [-in|i]                      Input file [stdin]
      [-paired|p]                  Output for reads with pairs [stdout]
      [-single|s]                  Output for reads without pairs [stderr]
=cut

