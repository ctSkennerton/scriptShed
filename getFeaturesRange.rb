#!/usr/bin/env ruby
require	'rubygems'
require 'bio'


ff = Bio::FlatFile.new(Bio::GenBank, ARGF)

# iterates over each GenBank entry
ff.each_entry do |gb|
    
    # shows accession and organism
	#puts " #{gb.accession} - #{gb.organism}"
    
    # iterates over each element in 'features'
  	gb.features.each do |annotation|
  		if annotation.feature =='CDS'
            hash = annotation.assoc
            locations = annotation.locations
            locus =  gb.entry_id
            #gene_id = locus + "-" + hash['note']
            name = hash['locus_tag'].to_s
            locations.each do |location|
                puts name << "\t" << location.from.to_s << "\t" << location.to.to_s << "\t" << hash['product'].to_s
#                puts gb.naseq.splicing(annotation.position) 
            end
        end
    end	
end






