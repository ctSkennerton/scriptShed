#!/bin/bash
# blast with the subject and query transposed

function print_help {
    echo "rbhb.sh [-n] -t <threads> -o <outfile> <prot1> <prot2> "
    echo "  prot1             Protein ORFs of genome 1"
    echo "  prot2             Protein ORFs of genome 2"
    echo "  -n                run ANI, default is AAI"
    echo "  -t threads        Number of threads to run blastp with."
    echo "                    WARNING: two blast jobs will be run simultaneously"
    echo "  -o outfile        Final output file"
    exit
}

threads=1
outfile="rbhb.out"
ani=0
blast=blastp
while getopts "nt:o:h" o; do       
    case "$o" in
        h)  print_help;;
        t) threads=$OPTARG;;
        o) outfile=$OPTARG;;
        n) ani=1; blast=blastn;;
    esac
done
shift $((OPTIND-1))
array=( $@ )
if [ ${#array[@]} != 2 ]; then
    print_help
fi

if (($ani == 1)); then
    if [ ! -e ${1}.nin ]; then
        makeblastdb -in $1 -dbtype nucl
    fi

    if [ ! -e ${2}.nin ]; then
        makeblastdb -in $2 -dbtype nucl
    fi
else
    if [ ! -e ${1}.psq ]; then
        makeblastdb -in $1 -dbtype prot
    fi

    if [ ! -e ${2}.psq ]; then
        makeblastdb -in $2 -dbtype prot
    fi
fi

$blast -evalue 1e-10 -query ${array[0]} -db ${array[1]} -outfmt 6 -out rbhb_A_B -num_threads $threads &

$blast -evalue 1e-10 -db ${array[0]} -query ${array[1]} -outfmt 6 -out rbhb_B_A -num_threads $threads &

wait
# the following perl command I downloaded from
# http://sysbio.harvard.edu/csb/resources/computational/scriptome/UNIX/Protocols/Sequences.html
# not sure exactly what they do but I'm assuming that they make sure that 
# the top blast hit is independant on whether the files are as the query or subject
perl -e '$name_col=0;
$score_col=11;
while(<>) {
    chomp;
    @F=split /\t/, $_;
    ($n, $s) = @F[$name_col, $score_col];
    if (! exists($max{$n})) {push @names, $n};
    if (! exists($max{$n}) || $s > $max{$n}) {$max{$n} = $s; $best{$n} = ()};
    if ($s == $max{$n}) {$best{$n} .= "$_\n"};
}
for $n (@names) {print $best{$n}}' rbhb_A_B >A_B.best

perl -e '$name_col=0; $score_col=11; while(<>) {
chomp; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col];
if (! exists($max{$n})) {push @names, $n};
if (! exists($max{$n}) || $s > $max{$n}) {$max{$n} = $s; $best{$n} = ()};
if ($s == $max{$n}) {$best{$n} .= "$_\n"};} for $n (@names) {print $best{$n}}' rbhb_B_A >B_A.best

perl -e '$col1=1; $col2=0;' -e '($f1,$f2)=@ARGV; open(F1,$f1);
while (<F1>) {s/\r?\n//; @F=split /\t/, $_; $line1{$F[$col1]} .= "$_\n"}
open(F2,$f2); while (<F2>) {s/\r?\n//;@F=split /\t/, $_;
if ($x = $line1{$F[$col2]}) {$x =~ s/\n/\t$_\n/g; print $x}}' A_B.best B_A.best > A_B_A

perl -e '$colm=0; $coln=13; $count=0;
while(<>) {s/\r?\n//; @F=split /\t/, $_;
if ($F[$colm] eq $F[$coln]) {print "$_\n"; $count++}}
warn "\nChose $count lines out of $. where column $colm had same text as column $coln\n\n";' A_B_A > $outfile

rm A_B.best B_A.best A_B_A rbhb_B_A rbhb_A_B

#perl -e '@cols=(0, 1, 3, 11); while(<>) {s/\r?\n//; @F=split /\t/, $_;
#print join("\t", @F[@cols]), "\n"} warn "\nJoined columns ",
#join(", ", @cols), " for $. lines\n\n"' A_B_A.recip > $5
