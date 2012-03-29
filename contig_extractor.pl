#!/usr/bin/perl
###############################################################################
#
#    this script extracts sequences from a multiple fasta file based on the 
#	 matches obtained from a blast file in m8 format
#
#    Copyright (C) 2010, 2011, 2012 Connor Skennerton
#
#    Support added for the Manotator by Mike Feb 2012
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
#use Getopt::Std;
use Pod::Usage;
#CPAN modules
use Getopt::Euclid;
#locally-written modules

BEGIN 
	{
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
	}

# get input params and print copyright

#my $options = checkParams();


my $query = \*STDIN;
if (defined $ARGV{'-i'}) {
    open($query, $ARGV{'-i'}) or die;
}
my $outfile = \*STDOUT; 
if (defined $ARGV{'-o'}) {
    open($outfile, ">", $ARGV{'-o'}) or die;
}

my %seqs;
#printAtStart();

if($ARGV{'-f'}) {
    my @aux = undef;
    my ($name, $seq, $qual);
    while (($name, $seq, $qual) = &readfq($query, \@aux)) {
        if($ARGV{'-Ri'}) {
            # perform the split according to the user
            my @p = split(/$ARGV{'-Ri'}->{separator}/, $name);
            $name = $p[$ARGV{'-Ri'}->{field_num}];
        }
        $seqs{$name} = 1;
    }
} elsif (! defined $ARGV{'-n'}) {
    while (my $line = <$query>) 
    {
        my $name = undef;
        chomp $line;
        if($ARGV{'-l'}) {
            $name = &list($line);
        } elsif ($ARGV{'-b'}) {
           $name = &blast($line);
        } elsif ($ARGV{'-s'}) {
            # skip header lines
            next if $line =~ /^@/;
            $name =&sam($line);
        } elsif ($ARGV{'-m'}) {
            $name =&mannotator("UniRef90_".$line);
        } else {
            next if ($line =~ /^\#/);
            last if ($line =~ /\#+FASTA/i);
            $name =&gff($line);
        }
        unless(defined $name ) {
            die "Crazy error!\n";
        } else {
            if($ARGV{'-Ri'}) {
                # perform the split according to the user
                my @p = split(/$ARGV{'-Ri'}->{separator}/, $name);
                $name = $p[$ARGV{'-Ri'}->{field_num}];
            }
            $seqs{$name} = 1;
        }
    }
} else {
    foreach my $e (@{$ARGV{'-n'}}) {
        if($ARGV{'-Ri'}) {
            # perform the split according to the user
            my @p = split(/$ARGV{'-Ri'}->[0]/, $e);
            $e = $p[$ARGV{'-Ri'}->[1]];
        }
        $seqs{$e} = 1;
    }
}
close $query;

my @aux = undef;
my ($name,$name2, $seq, $qual);
foreach my $database (@{$ARGV{'-d'}}) {
    open(DB,'<',$database) or die $!;
    while (($name, $seq, $qual) = readfq(\*DB, \@aux)) 
    {
        $name2 = $name;
        if($ARGV{'-Rd'}) {
            # perform the split according to the user
            my @p = split(/$ARGV{'-Rd'}->{separator_d}/, $name);
            $name = $p[$ARGV{'-Rd'}->{field_num_d}];
        }
        if (exists $seqs{$name})
        {
            unless($ARGV{'-v'})
            {
                print_seq(\$name2,\$seq,\$qual, $outfile);
            }
        }
        elsif ($ARGV{'-v'})
        {
                print_seq(\$name2,\$seq,\$qual, $outfile);
        }
    }
    close DB;
}

sub print_seq{
    my ($name_ref, $seq_ref, $qual_ref, $fh) = @_;
    if (defined $$qual_ref)
    {
        # fastq file
        print $fh "@".$$name_ref."\n".$$seq_ref."\n+".$$name_ref."\n".$$qual_ref."\n";
    }
    else
    {
        print $fh ">".$$name_ref."\n".$$seq_ref."\n";
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

sub list{
  	my ($line) = shift;
	return $line;
}

sub blast{
  my ($line) = shift;
  my @columns = split(/\t/, $line);
  	if (exists $ARGV{'-S'})
	{
		return $columns[1];
	}
	else
	{
		return $columns[0];
	}
}

sub sam{
    my ($line) = shift;
    my @c = split(/\t/,$line);
    # test whether the third bit is set - query unmapped
    unless($c[1] & 4)
    {
        # test whether the read is paired
        if (($c[1] & 1 ) && ($c[1] & 128)) {
            # test whether the read is the second pair
            return $c[0]."/2";
        } else {
            return $c[0]."/1";
        }
    }
}
sub gff {
    my ($line) = shift;
    my @c = split(/\t/, $line);
    return $c[0];
}
sub mannotator{
    my ($line) = shift;
    my @columns = split /\^/, $line;
    return $columns[0];
}

#sub checkParams 
#{
#    my %options;
#
#    # Add any other command line options, and the code to handle them
#    getopts( "i:d:so:c:lbShvfUg",\%options );
#
#    # if no arguments supplied print the usage and exit
#    #
#   pod2usage if (0 == (keys (%options) ));
#
#    # If the -h option is set, print the usage and exit
#    #
#    pod2usage if ($options{'h'});
#    unless ($options{'c'}) {
#        unless ($options{'S'} || $options{'g'} || $options{'U'} || $options{'b'} || $options{'l'} || $options{'f'} )
#        {
#            pod2usage('-msg' => "Please specify one of  -g -S -b -l -f -U");
#        }
#    }
#    unless ($options{'d'}) {
#        pod2usage('-msg' => "You must specify -d");
#    }
#    
#    if (defined $options{'s'} && !(defined $options{'b'}))
#    {
#        pod2usage('-msg' => "The subject flag can only be specified with the blast flag\n");
#    }
#
#    
#    return \%options;
#}


#sub printAtStart {
#print<<"EOF";
#---------------------------------------------------------------- 
# $0
# Copyright (C) 2010, 2011, 2012 Connor Skennerton
#    
# This program comes with ABSOLUTELY NO WARRANTY;
# This is free software, and you are welcome to redistribute it
# under certain conditions: See the source for more details.
#---------------------------------------------------------------- 
#EOF
#}

__END__

=head1 NAME
	
 contig_extractor

=head1 COPYRIGHT

 copyright (C) 2010, 2011, 2012 Connor Skennerton

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


=head1 REQUIRED ARGUMENTS

=over

=item -d <file>

file of reads where a subset needs to be extracted (can be FASTA or FASTQ, automatically detected).  Option can be specified multiple times.

=for Euclid:
    repeatable
    file.type: readable

=back

=head1 OPTIONS

=over

=item -i <input_file>

File containing sequences or identifiers to be extracted from the sequence database

=for Euclid:
    input_file.type: readable
    input_file.excludes: name
    input_file.excludes.error: You cannot specify both names on the command line and give an input file

=item -o <output_file>

Output file name

=for Euclid
    output_file.type: writable

=item -l

Input file is a list of identifiers, one per line

=item -b

Input file is in tabular blast format

=item -s

Input is in Sam format

=item -f

Input is in fasta or fastq format

=item -g

Input is in gff3 format

=item -m

Input is a mannotator formated annotations file

=item -S

Used only when the input is in blast format; sets the subject as the list of identifiers. Default: query

=item -n <name>...

A list of sequence names to extract in the form of a space separated list

=item -v

Invert the match. ie extract non-matching reads

=item -Ri <separator> <field_num>

The header in the input file contains additional information that must be removed.
This option takes two arguements, a separator, that can be any valid perl regular expression
to split the input identifier and the field to be used as the key for extraction (zero indexed).

=for Euclid
    field_num.type: i
    field_num.type.error: "Please specify an integer for the field_num"

=item -Rd <separator_d> <field_num_d>

The header in the database file contains additional information that can be removed.
Specify the options the same as option -Ri above

=for Euclid
    field_num_d.type: i
    field_num_d.type.error: "Please specify an integer for the field number"

=back

=head1 VERSION

 0.5.1

=head1 DESCRIPTION

   Used for extracting whole contigs / sequences from a multiple fasta file that contain 
   significant matches to reads/sequences/contigs from a variety of list formats

=cut

