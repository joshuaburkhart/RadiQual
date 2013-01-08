module RadPack
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
end
