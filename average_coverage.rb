#!/usr/bin/env ruby
#
require 'bio'
require 'bio-samtools'
require 'optparse'

options = { #TODO set defaults here
      :length => 1000,
      :outfile => $stdout
}
o = OptionParser.new do |opts|
      #TODO Fill in usage, description and option parsing below
    opts.banner = " Usage: average_coverage.rb [options] <bam file> <fasta file>\n"
    # Example option
    opts.on("-l", "--length INT", "minimum length for contig [default: #{options[:length]}]") do |f|
        options[:length] = f.to_i
    end
    opts.on("-o","--outfile FILE", "output file name [default: stdout]") do |f|
        options[:outfile] = File.open(f,"w")
    end
end
o.parse!
if ARGV.length != 2 #TODO require a set number of arguments?
    $stderr.puts o
    exit 1
end
# load in sequences
bam = Bio::DB::Sam.new(:bam => ARGV[0], :fasta => ARGV[1])
bam.open
contigs_file = Bio::FlatFile.auto(ARGV[1])
contigs_file.each_entry do |record|
  if record.nalen >= options[:length]
    coverage = bam.average_coverage(record.definition, 0, record.nalen)
    na_seq = Bio::Sequence::NA.new(record.seq)
    options[:outfile].puts "#{record.definition},#{record.nalen},#{na_seq.gc_content.to_f},#{coverage}"
  end
end
