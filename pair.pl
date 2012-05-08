#!/usr/bin/env perl
###############################################################################
#
#    pair
#    
#    lots of different operations that work with paired reads
#
#    Copyright (C) 2012 Connor Skennerton
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
checkParams();

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

sub seg_help {
print "pair segregate [-help|h] [-in|i] [-paired|p] [-single|s]

      [-help -h]                   Displays basic usage information
      [-in|i]                      Input file [stdin]
      [-paired|p]                  Output for reads with pairs [stdout]
      [-single|s]                  Output for reads without pairs [stderr]\n"

}

sub seg_main {
    #my $ARGV = shift;
    my @seg_options = ("help|h+", "paired|p:s","in|i:s", "single|s:s");
    my %options;

    # since I've shifted ARGV earlier the global will
    # be right for what I want to do
    GetOptions( \%options, @seg_options );

    if($options{'help'} || scalar keys %options == 0) { &seg_help; exit;}

    my $in_fh = \*STDIN;
    if(defined $options{'in'}) {
        open $in_fh, '<', $options{'in'} or die $!;
    }

    my $pair_fh = \*STDOUT;
    if(defined $options{'paired'}) {
        open $pair_fh, '>', $options{'paired'} or die $!;
    }

    my $single_fh = \*STDERR;
    if(defined $options{'single'}) {
        open $single_fh, '>', $options{'single'} or die $!;
    }

    my @aux = undef;
    my ($name, $seq, $qual);
    my %pairs_hash;
    # go through the input file and take note of which is the first and second read
    while (($name, $seq, $qual) = readfq($in_fh, \@aux)) {
        if($name =~ /(.*)\/(\d).*/) {
            push @{$pairs_hash{$1}}, [$name, $seq, $qual, $2];
        }
    }
    close $in_fh;

    # go through all the reads and determine which are paired and which are single
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


}

sub match_help {
    print "pair match [-help|h] [-in|i] -1 FILE -2 FILE [-out|o]

      [-help|h]                    Displays basic usage information
      [-in|i]                      Input file [stdin]
       -1 FILE                     Database for first member of pair
       -2 FILE                     Database for second member of pair
      [-out|o]                     Output file [stdout]\n";
}

sub match_print_seq {
    # pass in the seq and the mate
    # the mate is considered the second read
    my ($seq,$mate,$fh) = @_;
    
    # I'm assuming here that both the pairs are in the same format
    if(defined $seq->[1]) {
        printf($fh, "@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s\n", 
            $seq->[2], 
            $seq->[0], 
            $seq->[1],
            $mate->[2], 
            $mate->[0], 
            $mate->[1]);
    } else {
        printf($fh, ">%s\n%s\n>%s\n%s\n", 
            $seq->[2], 
            $seq->[0], 
            $mate->[2], 
            $mate->[0]);
    }
}
sub match_main {

    my @match_options = ( "help|h+", "in|i:s", "1:s", "2:s", "3:s", "out|o:s" );
    my %options;

    # since I've shifted ARGV earlier the global will
    # be right for what I want to do
    GetOptions( \%options, @match_options );

    if($options{'help'} || scalar keys %options == 0) { &match_help; exit;}

    my $in_fh; # = \*STDIN;
    if(defined $options{'in'}) {
        open $in_fh, '<', $options{'in'} or die $!;
    } else {
        $in_fh = \*STDIN;
    }
    my $out_fh;
    if(defined $options{'out'}) {
        open $out_fh, '>', $options{'out'} or die $!;
    } else {
        $out_fh = \*STDOUT;
    }

    my ($one_fh, $two_fh, $three_fh);
    if(defined $options{'3'}) {
        $three_fh = &openRead($options{'3'});
    } else {
        $one_fh = &openRead($options{'1'});
        $two_fh = &openRead($options{'2'});
    }

    #create two separate lists one for the first and second pair members
    my @aux = undef;
    my ($name, $seq, $qual);
    
    my %one_hash;
    my %two_hash;
    # go through the input file and take note of which is the first and second read
    while (($name, $seq, $qual) = readfq($in_fh, \@aux)) {
        # remove the trailing segment ID
        if($name =~ /(.*)\/(\d).*/) {
            if($2 == 1) {
                push @{$one_hash{$1.'/2'}}, [$seq, $qual, $name];
            } else {
                push @{$two_hash{$1.'/1'}}, [$seq, $qual, $name];
            }
        }
    }
    close $in_fh;

    # now go through each of the database files looking for
    # the corresponding mates
    while (($name, $seq, $qual) = readfq($one_fh, \@aux)) {
        if(defined $two_hash{$name}) {
            &match_print_seq([$seq,$qual,$name], $two_hash{$name}, $out_fh );
        }
    }
    # and now for the other file
    while (($name, $seq, $qual) = readfq($two_fh, \@aux)) {
        if(defined $one_hash{$name}) {
            &match_print_seq( $one_hash{$name}, [$seq,$qual,$name], $out_fh );
        }
    }




}
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    #
    
    # figure out subcommands
    my $arg = shift @ARGV;
    if ($arg eq "segregate") {
        &seg_main;   
    } elsif ($arg eq "match") {
        &match_main; 
    } else {
        pod2usage();
    }


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
    my ($default_value, $option_name, $options) = @_;
    if(exists $options->{$option_name})
    {
        return $options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

    pair

=head1 COPYRIGHT

   copyright (C) 2012 connor skennerton

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

    pair <command> [options]
    
    Commands:
        segregate       Given a single stream of reads, segregate into two streams of paired and unpaired reads
        match           Given an input stream of reads and two database files, get the other member of the pair

=cut

