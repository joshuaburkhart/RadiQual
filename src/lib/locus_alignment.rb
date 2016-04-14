module RadPack
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
      if cut_site_align != nil
        @cut_site_align = cut_site_align
        @offset = cut_site_align[:offset]
        @ref_strand_name = cut_site_align[:ref_strand_name]
        @name = "#{@ref_strand_name}-#{@offset}"
      else
        @name = '<unknown>'
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
      msg = ' '*2 + "LOCUS: #{@name} WITH #{@rad_tag_align_ary.size} RAD TAG(S):\n"
      unless @cut_site_align.nil?
        msg += ' '*6 + "CUT SITE: REFERENCE_STRAND = #{@cut_site_align[:ref_strand_name]}, OFFSET = #{@cut_site_align[:offset]} (#{@cut_site_align[:fr]}), SEQUENCE = #{@cut_site_align[:seq]}\n"
      end
      @rad_tag_align_ary.each { |tag|
        msg += ' '*6+"RAD TAG: NAME = #{tag[:name]}, OFFSET = #{tag[:offset]} (#{tag[:fr]}), SEQUENCE = #{tag[:seq]}\n"
      }
      msg
    end

  end
end
