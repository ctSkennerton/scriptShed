#!/usr/bin/perl -w

use strict;
use IO::File;
use Getopt::Std;
use vars qw/ %opt /;

sub usage()
{
       print "Usage: $0 -d [parent directory] -f [fasta file] \n";
       print "-d dir: parent dir\n";
       print "-f file: fasta file\n";
       exit;
}

sub init()
{
        my $opt_string = 'd:f:h';
        (@ARGV and getopts( "$opt_string", \%opt )) or usage();
        usage() if ($opt{'h'} or not $opt{'d'} or not $opt{'f'});
}

MAIN:
{
	init();	
	my $fasta_file = $opt{'f'};
	if($fasta_file =~ /\.gz$/gi)
	{
		die "Unzip file first\n";
	}
	else
	{
		open (FASTA, $fasta_file) or die "Error: $!\n";
	}	

	chdir $opt{'d'} || die "Could not chdir to ".($opt{'d'})."\n";
	
	my $file_ref = undef;
	my $id = undef;
		

	
	#writes one sequence to the file
	while(my $line = <FASTA>)
	{
		if($line =~ /^>([^\s]+)/gi)
		{
			$id = $1;
			if(defined $file_ref)
			{
				$file_ref->close();
			}

			$file_ref = IO::File->new($id.".fa", "w");
			print $file_ref $line;
		}
		elsif(defined $file_ref)
		{
			print $file_ref $line;
		}
	}
	#close the last file	
	$file_ref->close() if(defined $file_ref);
}
