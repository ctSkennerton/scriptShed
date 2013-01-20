#!/usr/bin/env ruby
#
require 'bio'
require 'bio-samtools'
require 'csv'
# array of contig names
contigs = {}
# load in sequences
contigs_file = Bio::FlatFile.auto(ARGV[1])
contigs_file.each_entry do |record|
  contigs[record.definition] = record.nalen
end
# only use contigs over this length
contigs.select! {|k,v| v >= 1000}

# create a csv file to hold the data
out_file = CSV.open(ARGV[2], "wb")

bam = Bio::DB::Sam.new(:bam => ARGV[0], :fasta => ARGV[1])
bam.open
contigs.each do |seq_name, seq_len|
  coverage = bam.chromosome_coverage(seq_name, 0, seq_len)
  coverage.unshift(seq_name)
  out_file << coverage
end

