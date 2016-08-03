#!/usr/bin/env perl
use warnings;
use strict;
use Class::Struct Seq => {name => '$', seq => '$', comment => '$', qual => '$' };


BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

my $query = \*STDIN;
if (defined $ARGV[0]) {
    if($ARGV[0] eq "-h" or $ARGV[0] eq "--help") {
        HELP_MESSAGE();
    } else {
        open($query, $ARGV[0]) or die;
    }
}

my $outfile = \*STDOUT; 

my @aux = undef;
while (my $seq = &readfq($query, \@aux)) {
    my $count = () = $seq->seq =~ /C..CH/gi;
    if($count > 0) {
        print $outfile $seq->name,"\t$count\n";
    }
}

sub HELP_MESSAGE {
    print "Usage:\nsearchMHC.pl <file.faa>\n\n";
    exit(1);
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
