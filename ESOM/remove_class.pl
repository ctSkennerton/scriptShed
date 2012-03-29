#!/sur/bin/env perl
#
use warnings;
use strict;
my %h;

open(CLS, '<',$ARGV[0] ) or die $!;
while (<CLS>){
    chomp;
    my @f = split(/\t/);
    $h{$f[0]} = 1;
}
close CLS;

my @a;
my $b;
open(IN, '<', $ARGV[1]) or die $!;
while(my $line = <IN>) {
    chomp $line;
    next if $line =~ /^%/;
    my @c = split(/\t/, $line);
    $b = scalar @c;
    if (not exists $h{$c[0]}) {
       push @a, $line;
   }
}
close IN;
print "%",scalar @a, "\n";
print "\%$b\n\%9";
for(my $i = 1; $i <= $b; $i++){
    print "\t1";
}
print "\n";
foreach (@a) {
    print "$_\n";
}
exit;
