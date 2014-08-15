#!/usr/bin/env perl
###############################################################################
#
#    fxrename
#    
#    Copyright (C) 2013 uqcskenn
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
use Carp;
use IO::Zlib;
use IO::File;
use IO::Uncompress::Bunzip2;
use Data::Dumper;

#CPAN modules

#locally-written modules
use Class::Struct Seq => {name => '$', seq => '$', comment => '$', qual => '$' };

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################

my $dfh;
if(scalar @ARGV >= 1) {
    if($global_options->{'gzip'}) {
        $dfh = IO::Zlib->new($ARGV[0],"rb") || die "$!: $ARGV[0]";
    } elsif($global_options->{'bzip2'}){
        $dfh = IO::Uncompress::Bunzip2->new($ARGV[0]) || die "$!: $ARGV[0]";
    }else {
       $dfh = IO::File->new($ARGV[0], 'r') || die "$!: $ARGV[0]";
    }
} else {
    $dfh = \*STDIN;
}

my $ofp;
if (scalar @ARGV >= 2) {
    $ofp = openWrite($ARGV[1]); #open($outfile, ">", $ARGV{'-o'}) or die;
} else {
    $ofp = \*STDOUT;
}

my @aux = undef;
while (my $seq = readfq($dfh, \@aux)) 
{
    my $modified_name = $seq->name();
    if(defined $global_options->{regex}) {
        my $eval_string = "\$modified_name=~ ".$global_options->{regex};
        eval $eval_string;
    }
    if(defined $global_options->{suffix}) {
        $modified_name .= $global_options->{suffix};
    }
    if(defined $global_options->{prefix}) {
        $modified_name = $global_options->{prefix}.$modified_name;
    }
    $seq->name($modified_name);
    print_seq(\$seq, $ofp);
}
######################################################################
# CUSTOM SUBS
######################################################################
sub format_seq {
    my $seq = shift;
    if (defined ${$seq}->qual) {
        return sprintf "@%s%s\n%s+\n%s\n", ${$seq}->name, (defined ${$seq}->comment ^ defined $global_options->{'comment'}) ? ${$seq}->comment : '', ${$seq}->seq, ${$seq}->qual;
    } else {
        if (defined ${$seq}->comment ^ defined $global_options->{'comment'}) {
            return sprintf ">%s%s\n%s", ${$seq}->name, ${$seq}->comment, ${$seq}->seq;
        }
        return sprintf ">%s\n%s", ${$seq}->name, ${$seq}->seq;
    }
}

sub fastaCut {
    #-----
    # Cut up a fasta sequence
    #
    my ($string, $line_wrap) = @_;
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
    my ($seq_ref, $fh) = @_;
    
    if(defined $global_options->{'wrap'})
    { 
        ${$seq_ref}->seq( fastaCut(${$seq_ref}->seq, $global_options->{'wrap'}) );
    }
    else
    {
        ${$seq_ref}->seq( ${$seq_ref}->seq."\n");
    }

    print $fh format_seq($seq_ref);

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
######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "prefix=s", "suffix=s", "regex|r=s", "comment|c+", "wrap|w:i", "gzip|z+", "bzip2|j+");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

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

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}


######################################################################
# MISC

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 fxrename
 Copyright (C) 2013 uqcskenn
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    fxrename

=head1 COPYRIGHT

   copyright (C) 2013 uqcskenn

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

    fxrename  [options] -regex|r REGEX <file.fx>...
    fxrename [options] [-prefix STR] [-suffix STR] <file.fx>...

      [-help -h]                   Displays basic usage information
      -r REGEX                     perl regular expression to be run over the sequence identifiers
      -suffix STR                  A string to be appended to the identifier
      -prefix STR                  A string to be prepended to the identifier
=cut

