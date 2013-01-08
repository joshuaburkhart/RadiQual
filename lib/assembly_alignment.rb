module RadPack
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
end
