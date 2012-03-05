#!/usr/bin/env perl
use warnings;
use strict;

if ($ARGV[0] =~ /[\-]\-h[elp]/) {
    print "Usage: unshufle.pl infile.fa forward_reads.fa reverse_reaads.fa\n";
    print "\tforward_reads.fa / reverse_reads.fa : reads to be unshuffled\n";
    print "\tinfile.fa :file containing shuffled fasta reads\n";
    print "\nBy default reads from STDIN and prints to STDOUT and STDERR\n";
    exit;
}

my ($in_fh, $out_1_fh, $out_2_fh);
if($ARGV[0]) {
    open $in_fh, '<', $ARGV[0] or die $!;
} else {
    $in_fh = \*STDIN;
}
if( $ARGV[1]) {
    open $out_1_fh, '>', $ARGV[1] or die $!;
} else {
    $out_1_fh = \*STDOUT;
}
if( $ARGV[2]) {
    open $out_2_fh, '>', $ARGV[2] or die $!;
} else {
    $out_2_fh = \*STDERR;
}

my @aux;
my $read_counter = 0;
while (my($name,$seq,$qual) = &readfq($in_fh,\@aux)) {
    $read_counter++;
    if($read_counter % 2) {
        # odd numbered read
        &print_seq(\$name,\$seq,\$qual,$out_1_fh);
    } else {
        #even numbered read
        &print_seq(\$name,\$seq,\$qual,$out_2_fh);
    }
}


#while(defined $lineA) {
#    print FILEA $lineA;
#    $lineA = <INFILE>;
#    while (defined $lineA && (substr($lineA, 0,1) != '>')) {
#        print FILEA $lineA;
#        $lineA = <INFILE>;
#    }
#
#    print FILEB $lineA;
#    $lineA = <INFILE>;
#    while (defined $lineA && (substr($lineA, 0,1) != '>')) {
#        print FILEB $lineA;
#        $lineA = <INFILE>;
#    }
#}

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
