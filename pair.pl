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
use IO::File;
use IO::Zlib;
use IO::Uncompress::Bunzip2;
use Class::Struct Seq => {name => '$', seq => '$', comment => '$', qual => '$', direction => '$' };
#CPAN modules
use Bio::Tools::CodonTable;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
checkParams();

# define a class for out sequences

# constants
use constant {FORWARD => 0, REVERSE => 1};


sub format_seq {
    my ($seq,$no_comments) = @_;
    if (defined ${$seq}->qual) {
        return sprintf "@%s%s\n%s+\n%s\n", ${$seq}->name, (defined ${$seq}->comment ^ defined $no_comments) ? ${$seq}->comment : '', ${$seq}->seq, ${$seq}->qual;
    } else {
        return sprintf ">%s%s\n%s", ${$seq}->name, (defined ${$seq}->comment ^ defined $no_comments) ? ${$seq}->comment : '', ${$seq}->seq;
    }
}

sub fastaCut {
    #-----
    # Cut up a fasta sequence
    #
    my ($string, $prot, $line_wrap) = @_;
    
    # translate if need be
    if(0 != $prot)
    {
        my $codon_table = Bio::Tools::CodonTable -> new ( -id => $prot );
        $string = $codon_table->translate($string);
    }
    
    # wrap the line if need be
    if(0 != $line_wrap)
    {
        my $return_str = "";
        my $len = length $string;
        my $start = 0;
        while($start < $len)
        {
            $return_str .= substr $string, $start, $line_wrap;
            $return_str .="\n";
            $start += $line_wrap;
        }
        return $return_str;
    }
    return "$string\n";
}

sub print_seq{
    my ($seq_ref, $fh, $prot, $line_wrap, $no_comments) = @_;
    
    if(defined $line_wrap)
    { 
        if(defined $prot)
        {
            ${$seq_ref}->seq( fastaCut(${$seq_ref}->seq, $prot, $line_wrap) );
        }
        else
        {
            ${$seq_ref}->seq( fastaCut(${$seq_ref}->seq, 0, $line_wrap) );
        }
    }
    elsif(defined $prot)
    {
         ${$seq_ref}->seq( fastaCut(${$seq_ref}->seq, $prot, 0) );
    }
    else
    {
        ${$seq_ref}->seq( ${$seq_ref}->seq."\n");
    }

    print $fh format_seq($seq_ref, $no_comments);

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
    my $current_seq = Seq->new();
	/^.(\S+)(.*)/;
    $current_seq->name($1);
    $current_seq->comment($2);
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
    $current_seq->seq($seq);
	return $current_seq if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
		chomp;
		$qual .= $_;
		if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
            $current_seq->qual($qual);
			return $current_seq;
		}
	}
	$aux->[1] = 1;
	return $current_seq;
}

sub guess_read_direction {
    my $seq_ref = shift;
    if(${$seq_ref}->name =~ /(.*)\.(f|r)$/) {
        ${$seq_ref}->direction(($2 eq 'f') ? FORWARD : REVERSE);
    } elsif (${$seq_ref}->name =~ /(.*)\/(\d).*/) {
        ${$seq_ref}->direction(($2 == 1) ? FORWARD : REVERSE);
    } elsif ( ${$seq_ref}->comment =~ /([12]):\w:\d:\w+/) {
        ${$seq_ref}->direction( ($1 == 1) ? FORWARD : REVERSE);
    } 
}

# newbler (.f or .r), illumina 1.3 (\1 or \2), illumina 1.8 - space separated
sub determine_pairing_convention {
    my $seq_ref = shift;
    if(${$seq_ref}->name =~ /(.*)\.([fr])$/) {
        ${$seq_ref}->direction(($2 eq 'f') ? FORWARD : REVERSE);
        return ($1, 'newbler', ($2 eq 'f') ? FORWARD : REVERSE);
    } elsif (${$seq_ref}->name =~ /(.*)\/(\d)$/) {
        ${$seq_ref}->direction(($2 == 1) ? FORWARD : REVERSE);
        return ($1, 'ill13', ($2 eq '1') ? FORWARD : REVERSE);
    } elsif ( ${$seq_ref}->comment =~ /([12]):\w:\d:\w+$/) {
        ${$seq_ref}->direction( ($1 == 1) ? FORWARD : REVERSE);
        return (${$seq_ref}->name, 'ill18', ($1 == 1) ? FORWARD : REVERSE);
    } else {
        ${$seq_ref}->direction(FORWARD);
        return (${$seq_ref}->name, 'unk', FORWARD);
    }

}

sub seg_help {
    print "pair segregate [-help|h] [-in|i FILE] { [-paired|p FILE] [-1 FILE] [-2 FILE] } [-single|s FILE] 

      [-help -h]                   Displays basic usage information
      [-in|i]                      Input file [stdin]
      [-paired|p]                  Output for reads with pairs [stdout]
      [-1]                         Output file for the first read in the pair [stdout]
      [-2]                         Output file for the second read in the pair [stdout]
      [-single|s]                  Output for reads without pairs [stderr] \n"

}

sub seg_main {
    #my $ARGV = shift;
    my @seg_options = ("help|h+", "1:s", "2:s", "paired|p:s","in|i:s", "single|s:s" );
    my %options;

    # since I've shifted ARGV earlier the global will
    # be right for what I want to do
    GetOptions( \%options, @seg_options );

    if($options{'help'}) { &seg_help; exit;}


    if (defined $options{'paired'} & (defined $options{'1'} | defined $options{'2'})) {
        &seg_help;
        exit;
    }

    my $one_fh = \*STDOUT;
    my $two_fh = \*STDOUT;
    if (defined $options{'1'} and defined $options{'2'}) {
        if ($options{'1'} eq $options{'2'}) {
            $options{'paired'} = $options{'1'};
            $options{'1'} = undef;
            $options{'2'} = undef;
        }
    }

    if (defined $options{'paired'} ) {
        $one_fh = IO::File->new($options{'paired'},'w') || die "cannot write to ".$options{"paired"}.": $!";
        $two_fh = $one_fh;
    } 
    
    if (defined $options{'1'}) {
        $one_fh = IO::File->new($options{'1'},'w') || die $!;
    } 

    if (defined $options{'2'}) {
        $two_fh = IO::File->new($options{'2'},'w') || die $!;
    } 


    my $in_fh = \*STDIN;
    if(defined $options{'in'}) {
        $in_fh = IO::File->new($options{'in'}, 'r') or die $!;
    }

    my $single_fh = \*STDERR;
    if(defined $options{'single'}) {
        $single_fh = IO::File->new($options{'single'}, 'w') or die $!;
    }

    my @aux = undef;
    my %pairs_hash;
    # go through the input file and take note of which is the first and second read
    while (my $current_seq  = readfq($in_fh, \@aux)) {
        my ($name, $type, $direction) = determine_pairing_convention(\$current_seq);
        if ($type eq 'unk') {
            warn "Cannot determine pairing convention for:".$current_seq->name."\n";
        }
        push @{$pairs_hash{$name}}, \$current_seq;
    }
    close $in_fh;

    # go through all the reads and determine which are paired and which are single
    print "\n";
    while(my($k,$v) = each %pairs_hash) {
        my $number_of_reads = scalar @{$v};
        if($number_of_reads > 1) {
            if($number_of_reads > 2) { 
                warn "more than one mapping for each mate for $k ($number_of_reads)\n";
                my $done_forward = 0;
                my $done_reverse = 0;
                foreach my $e (sort {${$a}->direction <=> ${$b}->direction} @{$v}) {
                    if (${e}->direction == FORWARD) {
                        unless($done_forward) {
                            print_seq($e, $one_fh, undef,undef,undef);
                            $done_forward = 1;
                        }
                    } else {
                        unless($done_reverse) {
                            print_seq($e, $two_fh, undef,undef,undef);
                            $done_reverse = 1;
                        }
                    }
                }
            } else {
                foreach my $e (sort {${$a}->direction <=> ${$b}->direction} @{$v}) {
                    (${$e}->direction == FORWARD) ? print_seq($e, $one_fh, undef,undef,undef) : print_seq($e, $two_fh,undef,undef,undef);
                }
            }
        } else {
            print_seq($v->[0], $single_fh,undef,undef,undef);
        }
    }
}

sub match_help {
    print "pair match [-help|h] [-in|i FILE] [-1 FILE] [-2 FILE] [-gzip|z] [-bzip2|j] [-append|a] -d1 FILE -d2 FILE

      [-help|h]                    Displays basic usage information
      [-in|i]                      Input file [stdin]
      [-1 FILE]                    Output for the first member of the pair [stdout]
      [-2 FILE]                    Output for the second member of the pair [stdout]
      [-gzip|z]                    The database files are gzipped
      [-bzip2|j]                   The database files are bzipped
      [-append|a]                  Append onto, rather than overwrite the files given with -1 -2
      -d1 FILE                     Database for the first member of the pair
      -d2 FILE                     Database for the second member of the pair\n";
}

sub match_main {

    my @match_options = ( "help|h+", "in|i:s", "gzip|z+","bzip2|j+","d1:s", "d2:s", "1:s", "2:s", "append|a+" );
    my %options;

    # since I've shifted ARGV earlier the global will
    # be right for what I want to do
    GetOptions( \%options, @match_options );

    if($options{'help'} || scalar keys %options == 0) { &match_help; exit;}

    my $in_fh = (defined $options{'in'}) ? openRead($options{'in'}) : \*STDIN; # = \*STDIN;
    #my $out_fh = (defined $options{'out'}) ? openWrite($options{'out'}) : \*STDOUT;


    my $one_fh = \*STDOUT;
    my $two_fh = \*STDOUT;
    my $d_one_fh;
    my $d_two_fh;
    if (defined $options{'d1'}) {
        $d_one_fh = openRead($options{'d1'}, $options{'gzip'}, $options{'bzip2'})# IO::File->new($options{'1'}, (defined $options{'append'}) ? 'a' : 'w') || die $!;
    } 

    if (defined $options{'d2'}) {
        $d_two_fh = openRead($options{'d2'}, $options{'gzip'}, $options{'bzip2'}) #IO::File->new($options{'2'},(defined $options{'append'}) ? 'a' : 'w') || die $!;
    } 
    if (defined $options{'1'}) {
        $one_fh = openWrite($options{'1'}, $options{'append'});
    } 
    if (defined $options{'2'}) {
        if($options{'2'} eq $options{'1'}) {
            $two_fh = $one_fh;
        } else  {
            $two_fh = openWrite($options{'2'}, $options{'append'});
        }
    } 

    #create two separate lists one for the first and second pair members
    my @aux = undef;
    
    my %one_hash;
    my %two_hash;
    # go through the input file and take note of which is the first and second read
    while (my $current_seq = readfq($in_fh, \@aux)) {
        # remove the trailing segment ID
        my ($name, $type, $direction) = determine_pairing_convention(\$current_seq);
        if($direction == FORWARD) {
            if($type eq 'ill13') {
                $name .= '/2';
            } elsif($type eq 'newbler') {
                $name .= '.r';
            }
            $one_hash{$name} = \$current_seq;
        } else {
            if($type eq 'ill13') {
                $name .= '/1';
            } elsif($type eq 'newbler') {
                $name .= '.f';
            }
            $two_hash{$name} = \$current_seq;
        }
    }
    close $in_fh;
    # now go through each of the database files looking for
    # the corresponding mates
    @aux = undef;
    while (my $current_seq = readfq($d_one_fh, \@aux)) {
        if(defined $two_hash{$current_seq->name}) {
            print_seq(\$current_seq, $one_fh, undef, undef,undef);
            print_seq($two_hash{$current_seq->name}, $two_fh, undef, undef,undef);
        }
    }
    # and now for the other file
    @aux = undef;
    while (my$current_seq = readfq($d_two_fh, \@aux)) {
        if(defined $one_hash{$current_seq->name}) {
            print_seq($one_hash{$current_seq->name}, $one_fh, undef, undef,undef);
            print_seq(\$current_seq, $two_fh, undef, undef,undef);
        }
    }
}

sub unshuffle_help {
    print "pair unshuffle [-help|h] [-in|i] [-1 FILE] [-2 FILE] 

      [-help|h]                    Displays basic usage information
      [-in|i FILE]                 Input file [stdin]
      [-1 FILE]                    Output for first member of pair [stdout]
      [-2 FILE]                    Output for second member of pair [stderr]\n";

}

sub unshuffle_main {

    my @unshuffle_options = ( "help|h+", "in|i:s", "1:s", "2:s" );
    my %options;

    # since I've shifted ARGV earlier the global will
    # be right for what I want to do
    GetOptions( \%options, @unshuffle_options );

    if($options{'help'} || scalar keys %options == 0) { &unshuffle_help; exit;}

    my $in_fh = (defined $options{'in'}) ? openRead($options{'in'}, undef, undef) : \*STDIN; # = \*STDIN;


    my $one_fh = \*STDOUT;
    my $two_fh = \*STDERR;

    if (defined $options{'1'}) {
        $one_fh = openWrite($options{'1'}, $options{'append'});
    } 

    if (defined $options{'2'}) {
        $two_fh = openWrite($options{'2'}, $options{'append'});
    } 

    my @aux = undef;
    
    my $read_counter = 1;
    while (my $current_seq = readfq($in_fh, \@aux)) {
        if($read_counter & 1) {
            print_seq(\$current_seq, $one_fh, undef, undef,undef)
        } else {
            print_seq(\$current_seq, $two_fh, undef, undef,undef)
        }
        $read_counter++;
    }
}


sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    #
    
    # figure out subcommands
    my $arg = shift @ARGV;
    if(! defined $arg) {
        pod2usage();
        exit;
    } elsif ($arg eq "segregate") {
        &seg_main;   
    } elsif ($arg eq "match") {
        &match_main; 
    } elsif ($arg eq "unshuffle") {
        &unshuffle_main;
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
    my ($fn,$append) = @_;
    my $fh = IO::File->new($fn,(defined $append) ? 'a' : 'w') || die "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn, $gzip, $bzip) = @_;
    my $fh;
    if(defined $gzip) {
        $fh = IO::Zlib->new($fn,"rb") || die $!;
    } elsif(defined $bzip){
        $fh = IO::Uncompress::Bunzip2->new($fn) || die $!;
    } else {
       $fh = IO::File->new($fn, 'r') || die $!;
    }
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
        unshuffle       Given an input file of paired reads that are interleaved, separate each member of the pair into a separate file

=cut

