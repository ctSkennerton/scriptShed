#!/usr/bin/env ruby

require 'bio'
require 'optparse'

format_identifiers = {
  'l' => 'locus_tag',
  'g' => 'gene',
  'p' => 'product',
  'n' => 'note',
  'f' => 'function',
  'P' => 'protein_id',
  'A' => 'accession'
}

options = {
  :format => 'l:g:p:n:f',
  :help => false,
  :split_records => false
}
o = OptionParser.new do |opts|
  opts.on('-f', '--format STR', "format string of identifiers separated by colons.\n#{format_identifiers.to_s}") do |arg|
    options[:format] = arg
  end
  opts.on('-h', '--help') do |opt|
    options[:help] = true
  end
  opts.on('-s', '--split-entries', "Output translations for each genbank entry to a different file, named after the accession of the genbank record") do |arg|
    options[:split_records] = true
  end
end
o.parse!

if options[:help]
  $stderr.puts o
  exit 1
end

format = options[:format].split(':')

ff = Bio::FlatFile.new(Bio::GenBank, ARGF)

# iterates over each GenBank entry
ff.each_entry do |gb|
  output = $stdout
  if options[:split_records]
    filename = gb.accession
    if filename.empty?
      filename = gb.entry_id
    end
    output = File.open("#{filename}.faa", 'w')
  end
  # iterates over each element in 'features'
  gb.features.each do |feature|
    position = feature.position
    hash = feature.assoc            # put into Hash

    # skips the entry if "/translation=" is not found
    next unless hash['translation']

    # collects gene name and so on and joins it into a string
    gene_info_part = []
    prepend_accession = false
    format.each do |i|
      if i != 'A'
        gene_info_part.push(hash[format_identifiers[i]])
      else
        prepend_accession = true
      end
    end
    gene_info = gene_info_part.compact.join(' ')

    # shows amino acid sequence in the database entry (/translation=)
    if prepend_accession
      output.puts ">#{gb.accession}_#{gene_info}"
    else
      output.puts ">#{gene_info}"
    end
    output.puts hash['translation']
  end
end
