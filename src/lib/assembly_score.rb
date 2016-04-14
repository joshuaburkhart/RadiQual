module RadPack
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
      if !@aligned_rad_tags.nil? && !@aligned_cut_sites.nil? && @aligned_cut_sites != 0
        (1000.0 * @aligned_rad_tags) / (2000.0 * @aligned_cut_sites)
      else
        -1
      end
    end

    def compareRadTags
      if !@aligned_rad_tags.nil?
        @aligned_rad_tags
      else
        -1
      end
    end

    def to_s
      alignments = getActOvrExpAlignments
      if alignments == -1
        alignments = "NO CUT SITES ALIGNED TO REFERENCE"
      end
      "FILE NAME: #{@name}\nNUMBER OF CUT SITES ALIGNED: #{@aligned_cut_sites}\nNUMBER OF RAD TAGS ALIGNED: #{@aligned_rad_tags}\nRAD SNP MISMATCHES: #{@rad_mismatches}\nACTUAL / EXPECTED: #{alignments}\n"
    end

    def to_f
      self.to_s +
          "BOWTIE CUT SITE OUTPUT:\n--\n#{@cut_output}--\nBOWTIE RAD TAG OUTPUT:\n--\n#{@rad_output}--\n"
    end
  end
end
