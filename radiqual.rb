#!/usr/bin/ruby

# This program has been tested with the following Ruby, Bowtie, and Samtools versions

# Ruby
# ruby 1.8.7 (2009-06-12 patchlevel 174) [i686-darwin9.7.0]
# ruby 2.0.0p648 (2015-12-16 revision 53162) [universal.x86_64-darwin15]

# Bowtie
# bowtie version 0.12.8

# Samtools
# Version: 0.1.18 (r982:295)

require 'time'
require 'optparse'
require "#{File.dirname(__FILE__)}/src/lib/radpack"

BOWTIE_TEST = %x(bowtie 2>&1)
SAM_TEST = %x(samtools 2>&1)
if BOWTIE_TEST.include?('command not found') || SAM_TEST.include?('command not found')
  puts 'Both bowtie and samtools must be installed and in the $PATH for this program to run.'
  puts 'Aborting...'
  exit
end

options ={}
optparse = OptionParser.new { |opts|
  opts.banner = <<-EOS
Usage: ruby radiqual.rb -c <cut site cohesive end sequence> -s <cut site sticky end sequence> -t </path/to/fasta/file/with/rad/tags> -o </path/to/output/dir> </path/to/fasta/file/with/contigs/1> [ ... </path/to/fasta/file/with/contigs/n>]

Example:     $ ./radiqual.rb -c C -s CGTAG -t data/input/mock_rad_tags.fasta -o data/output/ data/input/mock_contigs.fasta data/input/mock_contigs2.fasta
  EOS
  opts.on('-h', '--help', 'Display this screen') {
    puts opts
    exit
  }
  options[:rec_seq] = nil
  opts.on('-r', '--recognition_sequence SEQUENCE','Restriction enzyme recognition sequence') { |seq|
    options[:rec_seq] = seq}
  options[:rad_tag_file] = nil
  opts.on('-t', '--tag_file FILE', 'File Containing List of RAD Tags FILE') { |file|
    options[:rad_tag_file] = file
  }
  options[:out_dir] = 'radiqual_out'
  opts.on('-o', '--out_dir DIR', 'Directory to write output files to DIR') { |dir|
    if File.exists? dir
      unless File.directory? dir
        dir = "#{dir}_dir"
      end
    else
      %x(mkdir -p #{dir})
    end
    options[:out_dir] = dir
  }
}

optparse.parse!
if options[:rec_seq].nil?
  raise OptionParser::MissingArgument, "Recognition Sequence = \'#{options[:rec_seq]}\'"
elsif options[:rad_tag_file].nil?
  raise OptionParser::MissingArgument, "RAD Tag File = \'#{options[:rad_tag_file]}\'"
end

##########
# DRIVER #
##########

exec_id = Integer(Time.new)
rec_seq = options[:rec_seq]
rad_fasta_file = options[:rad_tag_file]
out_dir = options[:out_dir]
assem_score_file = "assembly_scores.#{exec_id}.txt"
assem_align_file = "assembly_aligns.#{exec_id}.txt"

puts "EXECUTION ID: #{exec_id}"
puts "RECOGNITION SEQ: #{rec_seq}"
puts "RAD FASTA FILE: #{rad_fasta_file}"
puts "OUTPUT DIR: #{out_dir}"
puts 'ASSEMBLY FILES: '

assembly_scores = Array.new
ARGV.each { |assem_file|
  puts assem_file
  assembly_scores << RadPack::AssemblyScore.new(assem_file)
}
puts

max_tag_length = -1

CONTIG_CORE_CAPTURE_RGX = /[atcgnATCGN]{100}([atcgnATCGN]+)[atcgnATCGN]{100}/
MISMATCH_CHAR = '>'
BEST = '--best'
ROUT = '--refout'

rad_tags = File.open(rad_fasta_file)
rad_fasta_line = rad_tags.gets
rad_tag_name = '<unknown>'
puts "rec_seq = #{rec_seq}"
puts "rec_seq.size = #{rec_seq.size}"
puts 'validating RAD tags...'
while rad_fasta_line
  print '.'
  STDOUT.flush
  if /^>/.match(rad_fasta_line)
    rad_tag_name = rad_fasta_line
  else
    test_seq = rad_fasta_line[0, rec_seq.size]
    if !/^[ATCG]/.match(test_seq)
      puts 'MALFORMED FASTA FILE DETECTED'
      puts "TAG NAME: #{rad_tag_name}"
      exit 1
    elsif test_seq != rec_seq
      rad_tags.close
      puts 'RAD TAG DOES NOT MATCH SPECIFIED CUT SITE'
      puts "TAG NAME: #{rad_tag_name}"
      puts "TAG FIRST #{rec_seq_size} BASES: #{test_seq}"
      puts "SPECIFIED CUT SITE SE SEQ: #{rec_seq}"
      exit 1
    end
  end
  if rad_fasta_line.size > max_tag_length
    max_tag_length = rad_fasta_line.size
  end
  rad_fasta_line = rad_tags.gets
end
puts "\nRAD tags match"
puts

rad_tags.close

assembly_align_ary = Array.new

#bowtie args
rec_seq_size = rec_seq.size
puts "rec_seq.size = #{rec_seq.size}"
puts "rec_seq_size = #{rec_seq_size}"
MAX_MISMATCHES = 2

puts 'aligning sequences to reference(s)...'
assembly_scores.each { |assembly_score|
  print '.'
  STDOUT.flush
  contigs_fa_fn = assembly_score.name
  core_contigs_fa_fn = "#{assembly_score.name}_core.fasta"
  contigs_fa_fh = File.open(contigs_fa_fn)
  core_contigs_fa_fh = File.open(core_contigs_fa_fn, 'w')
  header='< default contig core header'
  while contigs_fa_fn_line = contigs_fa_fh.gets
    if /^>/.match(contigs_fa_fn_line)
      puts contigs_fa_fn_line
      header = contigs_fa_fn_line.strip()
    else
      if CONTIG_CORE_CAPTURE_RGX.match(contigs_fa_fn_line)
        contig_core = $1
        if !contig_core.nil? and contig_core.length > 0
          core_contigs_fa_fh.print "#{header}\n"
          core_contigs_fa_fh.print "#{contig_core}\n"
        end
      end
    end
  end
  contigs_fa_fh.close
  core_contigs_fa_fh.close
  bowtie_idx_name = Time.new.to_f.to_s.sub('.', '_')
  bowtie_core_contigs_idx_name = "#{Time.new.to_f.to_s.sub('.', '_')}_core_contigs"
  sleep(1) #assuring a new Time
  %x(bowtie-build #{contigs_fa_fn} #{bowtie_idx_name})
  %x(bowtie-build #{core_contigs_fa_fn} #{bowtie_core_contigs_idx_name})
  assem_dir_name = contigs_fa_fn.strip.insert(0, "RadiQual_Track_#{bowtie_idx_name}").gsub('.', '_').gsub('/', '-').gsub(' ', '_')
  assem_dir = "#{out_dir}/#{assem_dir_name}"
  assem_vid = "#{assem_dir}/#{bowtie_idx_name}"
  if File.exists? assem_dir
    unless File.directory? assem_dir
      assem_dir = "#{assem_dir}_dir"
    end
  else
    %x(mkdir -p #{assem_dir})
  end

  puts "sam 1"
  assembly_score.setCoreCutResult(%x(bowtie -a -n0 -l#{rec_seq_size} -c #{bowtie_core_contigs_idx_name} #{rec_seq} 2>&1))
  %x(bowtie -a -n0 -l#{rec_seq_size} -c #{bowtie_idx_name} #{rec_seq} --sam #{assem_vid}.sam)
  %x(samtools view -bS #{assem_vid}.sam > #{assem_vid}_cutsites.bam)

  puts "sam 2"
  assembly_score.setCutResult(%x(bowtie -a -n0 -l#{rec_seq_size} -c #{bowtie_idx_name} #{rec_seq} 2>&1))
  %x(bowtie -a -n0 -l#{rec_seq_size} -c #{bowtie_idx_name} #{rec_seq} --sam #{assem_vid}.sam)
  %x(samtools view -bS #{assem_vid}.sam > #{assem_vid}_cutsites.bam)

  puts "sam 3"
  puts "bowtie #{bowtie_idx_name} -n#{MAX_MISMATCHES} -l#{rec_seq_size} #{BEST} -f #{rad_fasta_file}"

  puts "sam 3.4"
  %x(bowtie #{bowtie_idx_name} -n#{MAX_MISMATCHES} -l#{rec_seq_size} #{BEST} -f #{rad_fasta_file} > validation.debug.stdout 2> validation.debug.stderr)

  puts "sam 3.5"
  assembly_score.setRadResult(%x(bowtie #{bowtie_idx_name} -n#{MAX_MISMATCHES} -l#{rec_seq_size} #{BEST} -f #{rad_fasta_file} 2>&1))

  puts "sam 3.6"
  %x(bowtie #{bowtie_idx_name} -n#{MAX_MISMATCHES} -l#{rec_seq_size} #{BEST} -f #{rad_fasta_file} --sam #{assem_vid}.sam)

  puts "sam 3.7"
  %x(samtools view -bS #{assem_vid}.sam > #{assem_vid}_radtags.bam)

  puts "sam 4"
  %x(samtools merge #{assem_vid}_merged.bam #{assem_vid}_cutsites.bam #{assem_vid}_radtags.bam)

  puts "sam 5"
  %x(samtools sort #{assem_vid}_merged.bam #{assem_vid}_merged.sorted)

  puts "sam 6"
  %x(samtools index #{assem_vid}_merged.sorted.bam)

  puts "sam 7"
  %x(rm -f #{bowtie_idx_name}.*)
  %x(rm -f #{bowtie_core_contigs_idx_name}.*)
  %x(rm -f #{assem_vid}.sam)
  %x(rm -f #{assem_vid}_cutsites.bam)
  %x(rm -f #{assem_vid}_radtags.bam)
  %x(rm -f #{assem_vid}_merged.bam)

  assembly_align = RadPack::AssemblyAlignment.new(assembly_score.name)
  tmp_locus_ary = Array.new
  assembly_score.cut_output.each_line { |line|
    unless /^#|Reported/.match(line)
      line_ary = line.split
      fr = line_ary[1]
      ref_strand_name = line_ary[2]
      offset = Integer(line_ary[3])
      seq = line_ary[4]
      cut_site_align = RadPack::CutSiteAlignment.new(ref_strand_name, seq, offset, fr)
      locus_align = RadPack::LocusAlignment.new(cut_site_align)
      tmp_locus_ary << locus_align
    end
  }

  assembly_score.rad_output.each_line { |line|
    unless /^#|Reported/.match(line)
      puts "LINE #{line} END_LINE"
      line_ary = line.split
      tag_name = line_ary[0]
      puts "TAG_NAME #{tag_name} END_TAG_NAME"
      fr = line_ary[1]
      puts "FR #{fr} END_FR"
      ref_strand_name = line_ary[2]
      puts "REF_STRAND_NAME #{ref_strand_name} END_REF_STRAND_NAME"
      offset = Integer(line_ary[3])
      puts "OFFSET #{offset} END_OFFSET"
      seq = line_ary[4]
      puts "SEQ #{seq} END_SEQ"
      rad_tag_align = RadPack::RadTagAlignment.new(tag_name, ref_strand_name, seq, offset, fr)
      puts "RAD_TAG_ALIGN #{rad_tag_align} END_RAD_TAG_ALIGN"
      locus_name = ref_strand_name
      puts "LOCUS_NAME #{locus_name} END_LOCUS_NAME"
      if fr == '+'
        x = offset
        locus_name += "-#{x}"
      elsif fr == '-'
        x = (offset + seq.size)
        locus_name += "-#{x}"
      else
        puts "ERROR: RAD ALIGNMENT DIRECTION '#{fr}' UNRECOGNIZABLE"
        exit 1
      end

      tmp_idx = tmp_locus_ary.index { |l| l.name == locus_name }
      puts "TMP_IDX #{tmp_idx} END_TMP_IDX"
      if !tmp_idx.nil?
        locus_align = tmp_locus_ary[tmp_idx]
        locus_align.addRadTagAlign(rad_tag_align)
        if /length_([0-9]+)_cov/.match(ref_strand_name)
          n = Integer($1)
        else
          puts 'WARNING: IRREGULAR CONTIG NAME. REVERTING TO CMD LINE TOOLS FOR PARSING'
          if %x(uname).strip() == "Linux"
            n = Integer(%x(echo "$(cat #{assembly_score.name})>" | tr -d '\n' | grep -Po '(?<=#{ref_strand_name}).+(?=>)' | wc | awk -F' ' '{print $3}'))
          elsif %x(uname).strip() == "Darwin"
            n = Integer(%x(echo "$(cat #{assembly_score.name})>" | tr -d '\n' | perl -nle 'print $1 if m{(?<=#{ref_strand_name}).+(?=>)}' | wc | awk -F' ' '{print $3}'))
          else
            puts "ERROR: '#{%x(uname)}' OS not supported"
            exit(1)
          end
        end
        locus_align.n = n
      else
        puts 'WARNING: RAD TAG MATCHED LOCUS WITHOUT CUTSITE'
        puts " RAD TAG: #{rad_tag_align.to_s}"
        puts " LOCUS: #{locus_name}"
      end
    end
  }
  tmp_locus_ary.each { |l|
    if l.rad_tag_align_ary.size < l.expected
      if l.offset < max_tag_length
        l.expected -= 1
      end
      if l.rad_tag_align_ary.size < l.expected
        if l.n == 0
          ref_strand_name = l.ref_strand_name
        end
        if /length_([0-9]+)_cov/.match(ref_strand_name)
          n = Integer($1)
        else
          puts 'WARNING: IRREGULAR CONTIG NAME. REVERTING TO CMD LINE TOOLS FOR PARSING'
          if %x(uname).strip() == "Linux"
            n = Integer(%x(echo "$(cat #{assembly_score.name})>" | tr -d '\n' | grep -Po '(?<=#{ref_strand_name}).+(?=>)' | wc | awk -F' ' '{print $3}'))
          elsif %x(uname).strip() == "Darwin"
            n = Integer(%x(echo "$(cat #{assembly_score.name})>" | tr -d '\n' | perl -nle 'print $1 if m{(?<=#{ref_strand_name}).+(?=>)}' | wc | awk -F' ' '{print $3}'))
          else
            puts "ERROR: '#{%x(uname)}' OS not supported"
            exit(1)
          end
        end
        l.n = n
      end
      if l.n - l.offset < max_tag_length
        l.expected -= 1
      end
    end
    #puts "l.offset: #{l.offset}"
    #puts "l.expected: #{l.expected}"
    #puts "l.n: #{l.n}"
    if l.expected > 0
      ref_strand_name = l.ref_strand_name
      idx = assembly_align.contig_ary.index { |c| c.name == ref_strand_name }
      if !idx.nil?
        contig_align = assembly_align.contig_ary[idx]
        contig_align.addLocusAlign(l)
      else
        contig_align = RadPack::ContigAlignment.new(ref_strand_name)
        contig_align.addLocusAlign(l)
        assembly_align.addContigAlign(contig_align)
      end
    end
  }
  assembly_align_ary << assembly_align
}
puts "\nsequences aligned"
puts

assembly_scores.sort! { |i, j|
  [j.getActOvrExpAlignments, j.compareRadTags] <=> [i.getActOvrExpAlignments, i.compareRadTags]
}

puts 'writing results to file system...'
summary_file = File.open("#{out_dir}/#{assem_score_file}", 'w')
summary_file.print "ASSEMBLY SCORE SUMMARIES\n"
summary_file.print "========================\n\n"
summary_file.puts
assembly_scores.each { |assembly_score|
  puts
  puts assembly_score.to_s
  puts
  summary_file.print assembly_score.to_f
  summary_file.puts
}
summary_file.close
alignment_file = File.open("#{out_dir}/#{assem_align_file}", 'w')
assembly_align_ary.sort! { |i, j|
  (j.expected > 0 ? (j.actual / j.expected) : 0) <=> (i.expected > 0 ? (i.actual / i.expected) : 0)
}
alignment_file.print "ALIGNMENT SUMMARIES\n"
alignment_file.print "===================\n\n"
alignment_file.puts
assembly_align_ary.each { |assembly_alignment|
  alignment_file.print assembly_alignment.to_s
  alignment_file.puts
}
alignment_file.close
puts 'finished'
