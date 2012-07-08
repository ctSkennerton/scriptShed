#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'enumerator'
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
      if args[0] =~ /G(\d+)SP(\d+)/
        @gid = $1
        @spid = $2
      else
        raise ArgumentError, "argument does not fit the form: /G\d+SP\d+/"
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
end
o.parse!
if ARGV.length != 0 #TODO require a set number of arguments?
  $stderr.puts o
  exit 1
end
# Setup logging
Bio::Log::CLI.logger(options[:logger]) #bio-logger defaults to STDERR not STDOUT, I disagree
log = Bio::Log::LoggerPlus.new(LOG_NAME)
Bio::Log::CLI.configure(LOG_NAME)


#TODO what are you looking at me for? This is your script. Do something.
log.info 'e.g. logging'

spacers = Hash.new
targets = Hash.new{|h,k| h[k] = Array.new}
# parse the blast file and determine the locations of the spacers
Bio::Blast::Report.new(options[:blastfile]) do |report|
  report.each do |hit|
    unless spacers.has_key?(hit.query_id)
      sp = Spacer.new(hit.query_id)
      spacers[hit.query_id] = sp
    end
    targets[hit.target_id] << sp
    proto = ProtoSpacer.new(hit.target_id, hit.target_start, hit.target_end)
    spacers[hit.query_id] << proto
  end
end

# Now go though the sequences file and cut ~15bp from either side of the spacers for each CRISPR
Bio::FlatFile.open(Bio::FastaFormat, options[:sequences]) do |ff|
  ff.each do |record|
    targets[record.definition].each do |spacer|
      spacer.find_all{|i| i.genome == record.definition}.each do |protospacer|
        rstart = nil
        rend = nil
        lstart = nil
        lend = nil
        # fix things up if the spacer is on the negative strand
        if protospacer.start > protospacer.end
          lstart = protospacer.end - options[:length]
          lend = protospacer.end
          
          rstart = protospacer.start
          rend = protospacer.start + options[:length]
        else
          lstart = protospacer.start - options[:length]
          lend = protospacer.start
          
          rstart = protospacer.end
          rend = protospacer.end + options[:length]
        end
        if lstart > 0
          protospacer.left_flank = record.seq.subseq(lstart, lend)
        end
        if rend <= record.seq.length
          protospacer.right_flank = record.seq.subseq(rstart, rend)
        end
        protospacer.sequence = record.seq.subseq(lend,rstart)
      end
    end
  end
end

crisprs = Hash.new{|h,k| h[k] = Array.new}
spacers.each do |sp|
  crisprs[sp.gid] << sp
end

crisprs.each do |group, spacers|
  rfile = File.open("#{group}.protospacers.r.fa", "w")
  lfile = File.open("#{group}.protospacers.l.fa", "w")
  spacers.each do |spacer|
    spacer.each do |proto|
      header = ">G#{group}SP#{spacer.spid} Target=#{proto.genome} #{proto.start} #{proto.end};"
      if proto.right_flank
        rfile.puts "#{header} right flank (#{options[:length]}bp)"
        rfile.puts proto.right_flank
      end
      if proto.left_flank
        lfile.puts "#{header} left flank (#{options[:length]}bp)"
        lfile.puts left_flank
      end
    end
  end
end
# Feed everything into weblogo to make the images

