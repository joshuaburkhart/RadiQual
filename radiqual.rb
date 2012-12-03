#!/usr/bin/ruby

require 'time'
require 'optparse'

options ={}
optparse = OptionParser.new { |opts|
    opts.banner = <<-EOS
Usage: ruby rad_aligner.rb -c <cut site cohesive end sequence> -s <cut site sticky end sequence> -t </path/to/fasta/file/with/rad/tags> -o </path/to/output/dir> </path/to/fasta/file/with/contigs/1> [ ... </path/to/fasta/file/with/contigs/n>]

Example: ruby rad_aligner.rb -c C -s CGTAG -t ~/tmp/mock_rad_tags.fasta -o ~/tmp/ ~/tmp/mock_contigs.fasta
    EOS
    opts.on('-h','--help','Display this screen'){
        puts opts
        exit
    }
    options[:cohesive_end] = nil
    opts.on('-c','--cohesive SEQ','Cohesive End Sequence SEQ') { |seq|
        options[:cohesive_end] = seq
    }
    options[:sticky_end] = nil
    opts.on('-s','--sticky SEQ','Sticky End Sequence SEQ') { |seq|
        options[:sticky_end] = seq
    }
    options[:rad_tag_file] = nil
    opts.on('-t','--tag_file FILE','File Containing List of RAD Tags FILE') { |file|
        options[:rad_tag_file] = file
    }
    options[:out_dir] = "."
    opts.on('-o','--out_dir DIR','Directory to write output files to DIR') { |dir|
        options[:out_dir] = dir
    }
}

optparse.parse!
if(options[:cohesive_end].nil?)
    raise OptionParser::MissingArgument,"Cohesive End = \'#{options[:cohesive_end]}\'"
elsif(options[:sticky_end].nil?)
    raise OptionParser::MissingArgument,"Sticky End = \'#{options[:sticky_end]}\'"
elsif(options[:rad_tag_file].nil?)
    raise OptionParser::MissingArgument,"RAD Tag File = \'#{options[:rad_tag_file]}\'"
end

class AssemblyScore
    attr_accessor :name
    attr_accessor :aligned_cut_sites
    attr_accessor :aligned_rad_tags
    attr_accessor :cut_output
    attr_accessor :rad_output
    attr_accessor :rad_mismatches

    def initialize(name)
        @name = name
    end

    def setCutResult(cut_output)
        @cut_output = cut_output
        @cut_output.match(/d (\d+) a/)
        @aligned_cut_sites = Integer($1)
    end

    def setRadResult(rad_output)
        @rad_output = rad_output
        @rad_output.match(/t:\s+(\d+)\s+\(/)
                          @aligned_rad_tags = Integer($1)
                          @rad_mismatches = @rad_output.count(MISMATCH_CHAR)
    end

    def getActOvrExpAlignments
        if (!@aligned_rad_tags.nil? && !@aligned_cut_sites.nil? && @aligned_cut_sites != 0)
            return (1000.0 * @aligned_rad_tags) / (2000.0 * @aligned_cut_sites)
        else
            return -1
        end
    end

    def compareRadTags
        if(!@aligned_rad_tags.nil?)
            return @aligned_rad_tags
        else
            return -1
        end
    end

    def to_s
        alignments = getActOvrExpAlignments()
        if(alignments == -1)
            alignments = "NO CUT SITES ALIGNED TO REFERENCE"
        end
        "FILE NAME: #{@name}\nNUMBER OF CUT SITES ALIGNED: #{@aligned_cut_sites}\nNUMBER OF RAD TAGS ALIGNED: #{@aligned_rad_tags}\nRAD SNP MISMATCHES: #{@rad_mismatches}\nACTUAL / EXPECTED: #{alignments}\n"
    end

    def to_f
        self.to_s +
            "BOWTIE CUT SITE OUTPUT:\n--\n#{@cut_output}--\nBOWTIE RAD TAG OUTPUT:\n--\n#{@rad_output}--\n"
    end
end

CutSiteAlignment = Struct.new(:ref_strand_name,:seq,:offset,:fr)

RadTagAlignment = Struct.new(:name,:ref_strand_name,:seq,:offset,:fr)

class LocusAlignment
    attr_accessor :name
    attr_accessor :cut_site_align
    attr_accessor :offset
    attr_accessor :ref_strand_name
    attr_accessor :rad_tag_align_ary
    attr_accessor :actual
    attr_accessor :expected
    attr_accessor :n

    def initialize(cut_site_align=nil)
        if(cut_site_align != nil)
            @cut_site_align = cut_site_align
            @offset = cut_site_align[:offset]
            @ref_strand_name = cut_site_align[:ref_strand_name]
            @name = "#{@ref_strand_name}-#{@offset}"
        else
            @name = "<unknown>"
        end
        @rad_tag_align_ary = Array.new
        @actual = 0
        @expected = 2
        @n = 0 #reference size
    end

    def addRadTagAlign(rad_tag_align)
        @rad_tag_align_ary << rad_tag_align
        @actual += 1
    end

    def to_s
        msg = " "*2 + "LOCUS: #{@name} WITH #{@rad_tag_align_ary.size} RAD TAG(S):\n"
        if(!@cut_site_align.nil?)
            msg += " "*6 + "CUT SITE: REFERENCE_STRAND = #{@cut_site_align[:ref_strand_name]}, OFFSET = #{@cut_site_align[:offset]} (#{@cut_site_align[:fr]}), SEQUENCE = #{@cut_site_align[:seq]}\n"
        end
        @rad_tag_align_ary.each { |tag|
            msg += " "*6+"RAD TAG: NAME = #{tag[:name]}, OFFSET = #{tag[:offset]} (#{tag[:fr]}), SEQUENCE = #{tag[:seq]}\n"
        }
        return msg
    end

end

class ContigAlignment
    attr_accessor :name
    attr_accessor :locus_align_ary

    def initialize(name)
        @name = name
        @locus_align_ary = Array.new
    end

    def addLocusAlign(locus_align)
        @locus_align_ary << locus_align
    end

    def actual
        a = 0
        @locus_align_ary.each { |l|
            a += l.actual
        }
        return a
    end

    def expected
        e = 0
        @locus_align_ary.each { |l|
            e += l.expected
        }
        return e
    end

    #order loci by rad_tag_align_ary size
    def to_s
        @locus_align_ary.sort! { |i,j|
            j.rad_tag_align_ary.size <=> i.rad_tag_align_ary.size
        }
        msg  = " "*1+"CONTIG: NAME = #{@name}\n"
        msg += " "*1+"---------------" + "-"*@name.size + "\n"
        @locus_align_ary.each { |locus|
            msg += locus.to_s
        }
        return msg
    end
end

class AssemblyAlignment
    attr_accessor :name
    attr_accessor :contig_ary

    def initialize(name)
        @name = name
        @contig_ary = Array.new
    end

    def addContigAlign(contig)
        if(contig.locus_align_ary.size > 0)
            @contig_ary << contig
        end
    end

    def actual
        a = 0
        @contig_ary.each { |c|
            a += c.actual
        }
        return a
    end

    def expected
        e = 0
        @contig_ary.each { |c|
            e += c.expected
        }
        return e
    end

    #order contigs by locus_align_ary size
    def to_s
        @contig_ary.sort! { |i,j|
            j.locus_align_ary.size <=> i.locus_align_ary.size
        }
        msg  = "ASSEMBLY: NAME = #{@name}\n"
        msg += "=================" + "="*@name.size + "\n"
        @contig_ary.each { |contig|
            msg += "\n" + contig.to_s
        }
        return msg
    end
end

##########
# DRIVER #
##########

exec_id = Integer(Time.new)
ce_seq = options[:cohesive_end]
se_seq = options[:sticky_end]
rad_fasta_file = options[:rad_tag_file]
out_dir = options[:out_dir]
assem_score_file = "assembly_scores.#{exec_id}.txt"
assem_align_file = "assembly_aligns.#{exec_id}.txt"

puts "EXECUTION ID: #{exec_id}"
puts "COHESIVE END SEQ: #{ce_seq}"
puts "STICKY END SEQ: #{se_seq}"
puts "RAD FASTA FILE: #{rad_fasta_file}"
puts "OUTPUT DIR: #{out_dir}"
puts "ASSEMBLY FILES: "

assembly_scores = Array.new
ARGV.each { |assem_file|
    puts assem_file
    assembly_scores << AssemblyScore.new(assem_file)
}
puts

max_tag_length = -1;
MISMATCH_CHAR = ">"
BEST = "--best"
ROUT = "--refout"

rad_tags = File.open(rad_fasta_file)
rad_fasta_line = rad_tags.gets
rad_tag_name = "<unknown>"
cut_seq = "#{ce_seq}#{se_seq}"
se_seq_size = se_seq.size
puts "validating RAD tags..."
while(rad_fasta_line)
    print "."
    STDOUT.flush
    if(rad_fasta_line.match(/^>/))
        rad_tag_name = rad_fasta_line
    else
        test_seq = rad_fasta_line[0,se_seq_size]
        if(!test_seq.match(/^[ATCG]/))
            puts "MALFORMED FASTA FILE DETECTED"
            puts "TAG NAME: #{rad_tag_name}"
            exit 1
        elsif(test_seq != se_seq)
            rad_tags.close
            puts "RAD TAG DOES NOT MATCH SPECIFIED CUT SITE"
            puts "TAG NAME: #{rad_tag_name}"
            puts "TAG FIRST #{se_seq_size} BASES: #{test_seq}"
            puts "SPECIFIED CUT SITE SE SEQ: #{se_seq}"
            exit 1
        end
    end
    if(rad_fasta_line.size > max_tag_length)
        max_tag_length = rad_fasta_line.size
    end
    rad_fasta_line = rad_tags.gets
end
puts "\nRAD tags match"
puts

rad_tags.close

assembly_align_ary = Array.new

#bowtie args
cut_seq_size = cut_seq.size
MAX_MISMATCHES = 3

puts "aligning sequences to reference(s)..."
assembly_scores.each { |a|
    print "."
    STDOUT.flush
    contigs_fa_file = a.name
    bowtie_idx_name = Time.new.to_f.to_s.sub('.','_')
    sleep(1) #assuring a new Time
    %x(bowtie-build #{contigs_fa_file} #{bowtie_idx_name})
    assem_dir = contigs_fa_file.strip.sub(".","_").sub("/","-")
    assem_vid = "#{assem_dir}/#{bowtie_idx_name}" 
    %x(mkdir #{assem_dir})

    a.setCutResult(%x(bowtie -a -n0 -l#{cut_seq_size} -c #{bowtie_idx_name} #{cut_seq} 2>&1))
    %x(bowtie -a -n0 -l#{cut_seq_size} -c #{bowtie_idx_name} #{cut_seq} --sam #{assem_vid}.sam)
    %x(samtools view -bS #{assem_vid}.sam > #{assem_vid}_cutsites.bam)

    a.setRadResult(%x(bowtie #{bowtie_idx_name} -n#{MAX_MISMATCHES} -l#{se_seq_size} #{BEST} -f #{rad_fasta_file} 2>&1))
    %x(bowtie #{bowtie_idx_name} -n#{MAX_MISMATCHES} -l#{se_seq_size} #{BEST} -f #{rad_fasta_file} --sam #{assem_vid}.sam)
    %x(samtools view -bS #{assem_vid}.sam > #{assem_vid}_radtags.bam)

    %x(samtools merge #{assem_vid}_merged.bam #{assem_vid}_cutsites.bam #{assem_vid}_radtags.bam)
    %x(samtools sort #{assem_vid}_merged.bam #{assem_vid}_merged.sorted)
    %x(samtools index #{assem_vid}_merged.sorted.bam)
    %x(rm -f #{bowtie_idx_name}.*)

    assembly_align = AssemblyAlignment.new(a.name)
    tmp_locus_ary = Array.new
    a.cut_output.each_line { |line|
        if(!line.match(/^#|Reported/))
            line_ary = line.split
            fr = line_ary[1]
            ref_strand_name = line_ary[2]
            offset = Integer(line_ary[3])
            seq = line_ary[4]
            cut_site_align = CutSiteAlignment.new(ref_strand_name,seq,offset,fr)
            locus_align = LocusAlignment.new(cut_site_align)
            tmp_locus_ary << locus_align
        end
    }

    a.rad_output.each_line { |line|
        if(!line.match(/^#|Reported/))
            line_ary = line.split
            tag_name = line_ary[0]
            fr = line_ary[1]
            ref_strand_name = line_ary[2]
            offset = Integer(line_ary[3])
            seq = line_ary[4]
            rad_tag_align = RadTagAlignment.new(tag_name,ref_strand_name,seq,offset,fr)
            locus_name = ref_strand_name
            if(fr == "+")
                x = offset - ce_seq.size
                locus_name += "-#{x}"
            elsif(fr == "-")
                x = (offset + seq.size) - se_seq_size
                locus_name += "-#{x}"
            else
                puts "ERROR: RAD ALIGNMENT DIRECTION '#{fr}' UNRECOGNIZABLE"
                exit 1
            end

            tmp_idx = tmp_locus_ary.index { |l| l.name == locus_name}
            if(!tmp_idx.nil?)
                locus_align = tmp_locus_ary[tmp_idx]
                locus_align.addRadTagAlign(rad_tag_align)
                if(ref_strand_name.match(/length_([0-9]+)_cov/))
                    n = Integer($1)
                else
                    puts "WARNING: IRREGULAR CONTIG NAME. REVERTING TO CMD LINE TOOLS FOR PARSING"
                    n = Integer(%x(echo "$(cat #{a.name})>" | tr -d '\n' | grep -Po '(?<=#{ref_strand_name}).+(?=>)' | wc | awk -F' ' '{print $3}'))
                end
                locus_align.n = n
            else
                puts "WARNING: RAD TAG MATCHED LOCUS WITHOUT CUTSITE"
                puts " RAD TAG: #{rad_tag_align.to_s}"
                puts " LOCUS: #{locus_name}"

            end
        end
    }
    tmp_locus_ary.each { |l|
        if(l.rad_tag_align_ary.size < l.expected)
            if(l.offset < max_tag_length)
                l.expected -= 1
            end
            if(l.rad_tag_align_ary.size < l.expected)
                if(l.n == 0)
                    ref_strand_name = l.ref_strand_name
                    if(ref_strand_name.match(/length_([0-9]+)_cov/))
                        n = Integer($1)
                    else
                        puts "WARNING: IRREGULAR CONTIG NAME. REVERTING TO CMD LINE TOOLS FOR PARSING"
                        n = Integer(%x(echo "$(cat #{a.name})>" | tr -d '\n' | grep -Po '(?<=#{ref_strand_name}).+(?=>)' | wc | awk -F' ' '{print $3}'))
                    end
                    l.n = n
                end
                if(l.n - l.offset < max_tag_length)
                    l.expected -= 1
                end
            end
        end
        #puts "l.offset: #{l.offset}"
        #puts "l.expected: #{l.expected}"
        #puts "l.n: #{l.n}"
        if(l.expected > 0)
            ref_strand_name = l.ref_strand_name
            idx = assembly_align.contig_ary.index { |c| c.name == ref_strand_name}
            if(!idx.nil?)
                contig_align = assembly_align.contig_ary[idx]
                contig_align.addLocusAlign(l)
            else 
                contig_align = ContigAlignment.new(ref_strand_name)
                contig_align.addLocusAlign(l)
                assembly_align.addContigAlign(contig_align)
            end
        end
    }
    assembly_align_ary << assembly_align
}
puts "\nsequences aligned"
puts

assembly_scores.sort! { |i,j|
    [j.getActOvrExpAlignments,j.compareRadTags] <=> [i.getActOvrExpAlignments,i.compareRadTags]
}

puts "writing results to file system..."
summary_file = File.open("#{out_dir}/#{assem_score_file}",'w')
summary_file.print "ASSEMBLY SCORE SUMMARIES\n"
summary_file.print "========================\n\n"
summary_file.puts
assembly_scores.each { |a|
    puts
    puts a.to_s
    puts
    summary_file.print a.to_f
    summary_file.puts
}
summary_file.close
alignment_file = File.open("#{out_dir}/#{assem_align_file}",'w')
assembly_align_ary.sort! { |i,j|
    (j.actual / j.expected) <=> (i.actual / i.expected)
}
alignment_file.print "ALIGNMENT SUMMARIES\n"
alignment_file.print "===================\n\n"
alignment_file.puts
assembly_align_ary.each { |a|
    alignment_file.print a.to_s
    alignment_file.puts
}
alignment_file.close
puts "finished"
