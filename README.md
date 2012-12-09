RadiQual
========

a project to assess genome assembly quality using RAD sequence alignment

Visualization Quick Start
-------------------------

1. Review assembly scores produced by RadiQual.
2. Open IGV.
3. Load previously produced genome assembly fasta file as 'genome.'

    File -> Import Genome

4. Load .bam file produced by RadiQual for loaded assembly.

    File -> Load from File

    *NOTE: Be sure the index (.bam.bai) file is in the same directory as the .bam file loaded.

5. Review the assembly alignments produced by RadiQual for the loaded genome assembly.

    *NOTE: It may be helpful to see all the bases in the aligned track.

    Context Click (on track) -> Show all bases
