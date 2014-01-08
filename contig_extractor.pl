#!/usr/bin/perl
###############################################################################
#
#
#    Copyright (C) 2010-2013 Connor Skennerton
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
use IO::Zlib;
use IO::File;
use IO::Uncompress::Bunzip2;
use Data::Dumper;
#CPAN modules
use Getopt::Euclid;
use Bio::Tools::CodonTable;
#locally-written modules
use Class::Struct Seq => {name => '$', seq => '$', comment => '$', qual => '$' };


BEGIN 
	{
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
	}

my $query = \*STDIN;
if (defined $ARGV{'-i'}) {
    open($query, $ARGV{'-i'}) or die;
}
my $outfile = \*STDOUT; 
if (defined $ARGV{'-o'}) {
    open($outfile, ">", $ARGV{'-o'}) or die;
}

my %seqs;
my %cluster_map;
#printAtStart();
if(! defined $ARGV{'-r'}) {
    if($ARGV{'-f'}) {
        my @aux = undef;
        while (my $seq = &readfq($query, \@aux)) {
            my $name2 = $seq->name;
            if($ARGV{'-Ri'}) {
                $name2 =~ s/$ARGV{'-Ri'}/$1/;
            }
            $seqs{$name2} = 1;
        }
    } elsif (! defined $ARGV{'-n'}) {
        if (defined $ARGV{'-c'}){
            &mapping($ARGV{'-c'}, \%seqs);
        } else {

            while (my $line = <$query>) 
            {
                chomp $line;
                if ($line =~ /^$/) {
                    next;
                }
                my $name = undef;
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
                if(! defined $name ) {
                    die "Crazy error!\n";
                } else {
                    if($ARGV{'-Ri'}) {
                        $name =~ s/$ARGV{'-Ri'}/$1/;
                    }
                    $seqs{$name} = 1;
                }
            }
        }
    } elsif(defined $ARGV{'-n'}) {
        foreach my $e (@{$ARGV{'-n'}}) {
            if($ARGV{'-Ri'}) {
                $e =~ s/$ARGV{'-Ri'}/$1/;
            }
            $seqs{$e} = 1;
        }
    }
}
close $query;
my $keys_to_find = scalar keys %seqs;
if(defined $ARGV{'-r'}) {
    $keys_to_find = 1;
}

foreach my $database (@{$ARGV{'-d'}}) {
    # check to see if there are any keys left
    last unless($keys_to_find > 0);
    
    my $dfh;
    if($ARGV{'-z'}) {
        $dfh = IO::Zlib->new($database,"rb") || die $!;
    } elsif($ARGV{'-j'}){
        $dfh = IO::Uncompress::Bunzip2->new($database) || die $!;
    }else {
       $dfh = IO::File->new($database, 'r') || die $!;
    }

    my @aux = undef;
    while (my $seq = readfq($dfh, \@aux)) 
    {
        if(defined $ARGV{'-r'}) {
            if($seq->name =~ /$ARGV{'-r'}/ || $seq->comment =~ /$ARGV{'-r'}/) {
                unless($ARGV{'-v'})
                {
                    print_seq(\$seq, $outfile);
                }
            }
            elsif ($ARGV{'-v'})
            {
                print_seq(\$seq, $outfile);
            }
            next;
        }
        
        my $name2 = $seq->name;
        if($ARGV{'-Rd'}) {
            # perform the split according to the user
            $name2 =~ s/$ARGV{'-Rd'}/$1/;
        }
        if (exists $seqs{$name2})
        {
            unless($ARGV{'-v'})
            {
                if (defined $ARGV{'-c'}) {
                    $outfile = $seqs{$name2};
                }
                print_seq(\$seq, $outfile);
                $keys_to_find-- unless ($ARGV{'-Force'});
            }
        }
        elsif ($ARGV{'-v'})
        {
            print_seq(\$seq, $outfile);
        }
        # check to see if there are any keys left
        last unless($keys_to_find > 0);
    }
    $dfh->close();
}
if (defined $ARGV{'-c'}) {
    while (my ($seq,$fh) = each %seqs) {
        close $fh;
    }
}

sub format_seq {
    my $seq = shift;
    if (defined ${$seq}->qual) {
        return sprintf "@%s%s\n%s+\n%s\n", ${$seq}->name, (defined ${$seq}->comment ^ defined $ARGV{'-C'}) ? ${$seq}->comment : '', ${$seq}->seq, ${$seq}->qual;
    } else {
        if (defined ${$seq}->comment ^ defined $ARGV{'-C'}) {
            return sprintf ">%s%s\n%s", ${$seq}->name, ${$seq}->comment, ${$seq}->seq;
        }
        return sprintf ">%s\n%s", ${$seq}->name, ${$seq}->seq;
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
    my ($seq_ref, $fh) = @_;
    
    if(defined $ARGV{'-w'})
    { 
        if(defined $ARGV{'-p'})
        {
            ${$seq_ref}->seq( fastaCut(${$seq_ref}->seq, $ARGV{'-p'}, $ARGV{'-w'}) );
        }
        else
        {
            ${$seq_ref}->seq( fastaCut(${$seq_ref}->seq, 0, $ARGV{'-w'}) );
        }
    }
    elsif(defined $ARGV{'-p'})
    {
         ${$seq_ref}->seq( fastaCut(${$seq_ref}->seq, $ARGV{'-p'}, 0) );
    }
    else
    {
        ${$seq_ref}->seq( ${$seq_ref}->seq."\n");
    }

    print $fh format_seq($seq_ref);

    #if (defined ${$seq_ref}->qual ^ defined $ARGV{'-F'})
    #{
        ## fastq file
        #print $fh "@".$$name_ref."\n".$seq."+".$$name_ref."\n".$$qual_ref."\n";
    #}
    #else
    #{
        #print $fh ">".$$name_ref."\n".$seq;
    #}
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

sub sam {
    my ($line) = shift;
    my @c = split(/\t/,$line);
    # test whether the third bit is set - query unmapped
    unless($c[1] & 4)
    {
        # test whether the read is paired
        if (($c[1] & 1 ) && ($c[1] & 128)) {
            # test whether the read is the second pair
            if($ARGV{'-I'}) { 
            	return $c[0]."/2";
            } else {
            	return $c[0];
            }
        } else {
            if ($ARGV{'-I'}) {
            	return $c[0]."/1";
            } else {
            	return $c[0];
            }
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

sub mapping {
    my ($mapfile, $seqs) = @_;
    my $fh = IO::File->new($mapfile, 'r') || die $!;
    my %seen_file;
    my $fh_count = 0;
    while(<$fh>) {
        chomp;
        my @mapping = split /\t/;
        if (exists $seen_file{$mapping[1]}) {
            $seqs->{$mapping[0]} = $seen_file{$mapping[1]}
        } else {
            $fh_count++;
            if ($fh_count > 255) {
                warn "there are more than 255 different groupings specified in the mappings file\n";
                warn "please split the mapping file such that there are less groupings then run\n";
                warn "contig_extractor using the split mapping files\n";
                exit 1;
            }
            my $new_fh = IO::File->new($mapping[1],'w') || die $!;
            $seen_file{$mapping[1]} = $new_fh;
            $seqs->{$mapping[0]} = $new_fh;
        }
    }
}


__END__

=head1 NAME
	
 contig_extractor

=head1 COPYRIGHT

 copyright (C) 2010-2013 Connor Skennerton

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

=item -h

print a short help msg

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

=item -I

If using read names from a sam input file, append Illumina style (/1 or /2). [Default: false]

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

=item -r <regex>

Specify a perl regular expression that will be matched against all sequences in the database.
This option is not compatible with other types of input specified with -n or -i.  Unlike
other methods this regex will match against both the name and the comment sections of the 
fasta record

=for Euclid
    regex.excludes: name, input_file

=item -S

Used only when the input is in blast format; sets the subject as the list of identifiers. Default: query

=item -n <name>...

A list of sequence names to extract in the form of a space separated list

=item -c <mapping_file>

    File containing a mapping of contigs and their groupings.  Each line must contain a single contig name
    and the name of an output file for that contig that is TAB sepatated. 
    Example:
    contig0001	cluster_1.fa
    contig0002	cluster_2.fa
    NOTE: You cannot specify more than 255 different output file names. If you have more groupings than this
    you will need to split your job into parts 

=for Euclid
    mapping_file.type: readable

=item -v

Invert the match. ie extract non-matching reads

=item -z

The database file is gziped

=item -j

The database file is bziped

=item -F

Force output to be in fasta format

=item -Ri <input_regex>

The header in the input file contains additional information that must be removed.
Specify a perl regular expression such that $1 (first capture group) contains the required information.

=item -Rd <database_regex>

The header in the database file contains additional information that can be removed.
Specify a perl regular expression such that $1 (first capture group) contains the required information.

=item -w [<wrap_length>]

Wrap the lines at wrap-length chars

=for Euclid
    wrap_length.type: i

=item -p [<protein_code>]

    Enter a code if you wish to convert to protein
    Specify a number from the following list (Uses: Bio::Tools::CodonTable)
    1 Standard
    2 Vertebrate Mitochondrial
    3 Yeast Mitochondrial
    4 Mold, Protozoan,_and_CoelenterateMitochondrial_and_Mycoplasma/Spiroplasma
    5 Invertebrate Mitochondrial
    6 Ciliate, Dasycladacean_and_Hexamita_Nuclear
    9 Echinoderm Mitochondrial
    10 Euplotid Nuclear
    11 Bacterial
    12 Alternative Yeast_Nuclear
    13 Ascidian Mitochondrial
    14 Flatworm Mitochondrial
    15 Blepharisma Nuclear
    16 Chlorophycean Mitochondrial
    21 Trematode Mitochondrial
    22 Scenedesmus obliquus_Mitochondrial
    23 Thraustochytrium Mitochondrial

Default: no conversion

=for Euclid
    protein_code.type: i

=item -C 

Do not print comments in fasta files

=item -Force

Force scanning the entire database file.  contig_etractor tries to speed up processing by keeping track
of the patterns that have been found and exiting when they all have been, even if the end of the database
file has not been reached.  However the combination of other options such as -Ri or -Rd means that
sometimes the patterns are not unique and so scanning should be forced till the end of the file to 
ensure that all matches are found
be unique, however if that is not the case 

=back

=head1 VERSION

 0.5.9

=head1 DESCRIPTION

   Used for extracting whole contigs / sequences from a multiple fasta file that contain 
   significant matches to reads/sequences/contigs from a variety of list formats

=cut

