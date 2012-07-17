#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'enumerator'
require 'bio'
# small class to hold the information about how each 
# spacer matches a particular genome
class ProtoSpacer
  attr_accessor :right_flank, :left_flank, :sequence, :genome, :start, :end
  def initialize(g, s, e)
    @genome = g
    @start = s
    @end = e
    @right_flank = nil
    @left_flank = nil
  end # initialize

end # Protospacer


# small class to hold all the hits for our spacer
class Spacer
  include Enumerable 
  attr_accessor :protospacers, :gid, :spid 
  
  def initialize(*args)
    if args.size == 1
      # split name here
      if args[0] =~ /G(\d+)(SP|FL)(\d+).*/
        @gid = $1
        @spid = $3
      else
        raise ArgumentError, "argument #{args[0]} does not fit the form: /G\d+SP\d+/"
      end
    elsif args.size == 2
      @gid = args[0]
      @spid = args[1]
    else
      raise ArgumentError, "wrong number of arguments (#{args.size} for 2"
    end
    @protospacers = Array.new
  end # initialize
  
  def each
    @protospacers.each{|i| yield i}
  end

  def <<(elem)
      @protospacers << elem
  end

end # Spacer


SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')

# Parse command line options into the options hash
options = { #TODO set defaults here
  :logger => 'stderr',
  :blastfile => nil,
  :sequences => nil,
  :length => 15
}
o = OptionParser.new do |opts|
  #TODO Fill in usage, description and option parsing below
  opts.banner = "
    Usage: #{SCRIPT_NAME} <arguments>
    
    Description of what this program does...\n"
  # Example option
  opts.on("-e", "--eg", "description [default: #{options[:eg]}]") do |f|
    options[:operation] = OVERALL
  end
  
  # logger options
  opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") do |q|
    Bio::Log::CLI.trace('error')
  end
  opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") do | name |
    options[:logger] = name
  end
  opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG") do | s |
    Bio::Log::CLI.trace(s)
  end
  opts.on("-b", "--blastfile FILE", "Input blast outfile in (6/8 format) containing spacer matches to phage contigs") do |b|
    options[:blastfile] = b
  end
  opts.on("-c", "--contigs FILE", "Input fasta file containing phage contigs") do |c|
    options[:sequences] = c
  end
  opts.on("-l","--length INT", "length of the flanking region of the protospacer [default: #{options[:length]}") do |l|
    options[:length] = l.to_i
  end
end
o.parse!
unless options[:blastfile] and options[:sequences] #TODO require a set number of arguments?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
log = Bio::Log::LoggerPlus.new(LOG_NAME)
Bio::Log::CLI.configure(LOG_NAME)


#TODO what are you looking at me for? This is your script. Do something.
#log.info 'e.g. logging'

spacers = Hash.new
targets = Hash.new{|h,k| h[k] = Array.new}
# parse the blast file and determine the locations of the spacers
File.open(options[:blastfile],"r") do |report|
  report.each do |hit|
    fields = hit.split(/\t/)
    unless spacers.has_key?(fields[0])
      sp = Spacer.new(fields[0])
      spacers[fields[0]] = sp
    end
    targets[fields[1]] << spacers[fields[0]]
    proto = ProtoSpacer.new(fields[1], fields[8].to_i, fields[9].to_i)
    spacers[fields[0]] << proto
  end
end

# Now go though the sequences file and cut ~15bp from either side of the spacers for each CRISPR
Bio::FlatFile.open(Bio::FastaFormat, options[:sequences]) do |ff|
  ff.each do |record|
    if targets.has_key?(record.definition)
      targets[record.definition].each do |spacer|
        spacer.find_all{|i| i.genome == record.definition}.each do |protospacer|
          rstart = nil
          rend = nil
          lstart = nil
          lend = nil
          sequence = Bio::Sequence::NA.new(record.seq)
          # fix things up if the spacer is on the negative strand
          if protospacer.start > protospacer.end
            lstart = protospacer.end - options[:length]
            lend = protospacer.end

            rstart = protospacer.start + 1
            rend = protospacer.start + options[:length]
            # reverse complement the sequence
            #sequence = Bio::Sequence::NA.new(record.seq)
            sequence.reverse_complement!# = Bio::Sequence::NA.new(record.seq).reverse_complement
          else
            lstart = protospacer.start - options[:length]
            lend = protospacer.start

            rstart = protospacer.end + 1
            rend = protospacer.end + options[:length]
          end
          if lstart > 0 
            protospacer.left_flank = sequence.subseq(lstart, lend)
          end
          if rend <= record.seq.length
            protospacer.right_flank = sequence.subseq(rstart, rend)
          end
          protospacer.sequence = sequence.subseq(lend + 1,rstart - 1)
          #puts "#{record.seq.subseq(lstart,rend)}"
          #puts "#{protospacer.left_flank}#{protospacer.sequence}#{protospacer.right_flank}"
          #puts
        end
      end
    end
  end
end

crisprs = Hash.new{|h,k| h[k] = Array.new}
spacers.each do |spid,sp|
  crisprs[sp.gid] << sp
end

crisprs.each do |group, spacers|
  rfile = File.open("G#{group}.protospacers.r.fa", "w")
  lfile = File.open("G#{group}.protospacers.l.fa", "w")
  spacers.each do |spacer|
    spacer.each do |proto|
      header = ">G#{group}SP#{spacer.spid} Target=#{proto.genome} #{proto.start} #{proto.end};"
      if proto.right_flank
        rfile.puts "#{header} right flank (#{options[:length]}bp)"
        rfile.puts proto.right_flank
      end
      if proto.left_flank
        lfile.puts "#{header} left flank (#{options[:length]}bp)"
        lfile.puts proto.left_flank
      end
    end
  end
end
# Feed everything into weblogo to make the images

