#!/usr/bin/env perl
#
use warnings;
use strict;

my %a;
my %n;
my $i = 0;
my $j = 0;
open(IN,'<',$ARGV[0]) or die "$!\n\nUsage: $0 INPUT_FILE >file.lrn 2>file.names\n";
my $first_line = <IN>;
while(<IN>) {
    $i++;
    chomp;
    my @b = split(/\t/);
    $j = scalar @b;
    $n{$i} = shift @b;
    $a{$i} = join("\t",@b);
}

my @s = sort {$a <=> $b} keys %a;
print "%",scalar @s, "\n";
print "\%$j\n\%9";
for(my $x = 1;$x <= $j; $x++) {
    print "\t1";
}
print "\n\%$first_line\n";
foreach (@s) {
    print "$_\t$a{$_}\n";
    warn "$_\t$n{$_}\n";
}
