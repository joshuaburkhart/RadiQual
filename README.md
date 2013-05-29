RadiQual
========

a project to assess genome assembly quality using RAD sequence alignment

Instructions
------------

1. Filter RAD tags using the correct command from Radiqual/filter_tags.txt.
   
        $ cat /home11/mmiller/Wyeomyia/output/ustacks_out/KC.tags.tsv | grep -P '^[0-9]+[\s]+[0-9]+[\s]+[0-9]+[\s]+[a-z]+[\s]+[0-9]+[\s]+[0-9a-zA-Z_-]+[\s]+[ACTG]+[\s]+' | awk -F' ' '{print $NF}' | uniq | sort > /home11/mmiller/Wyeomyia/output/ustacks_out/rad_tags.txt

2. Review http://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites or http://rebase.neb.com/rebase/rebase.html to find the correct restriction enzyme cut (recognition), cohesive end, and sticky end sequences.

3. Confirm actual restriction enzyme sequence matches expected.

        $ head -n1 rad_tags.txt && tail -n1 rad_tags.txt

4. Index filtered RAD tags with Radiqual/index_tags.sh.

        $ index_tags.sh rad_tags.txt

5. Ensure both bowtie and samtools are installed and in the $PATH as they are required for this program to run.

        $ qsub -IXl nodes=1:ppn=32 -q fatnodes

        $ module load bowtie

        $ bowtie --version

        $ module load samtools

        $ samtools

6. Review Radiqual/radiqual.rb options and example usage with -h.

        $ radiqual.rb -h

7. Execute Radiqual/radiqual.rb with correct options.

        *NOTE: This program has been tested using Velvet 1.2.07 contigs.fa files.

        $ radiqual.rb -c CC -s TGCAGG -t /home11/mmiller/Wyeomyia/output/ustacks_out/rad_tags.idx.txt -o /home11/mmiller/Wyeomyia/output/radiqual_out/ /home11/mmiller/Wyeomyia/output/velvet_out/velvet-wy_unfiltered_combined_reads_diginorm_paired.fastq.keep_k\=55_e\=19/contigs.fa

        *NOTE: Multiple genome assemblies can be compared as below.

        $ radiqual.rb -c CC -s TGCAGG -t /home11/mmiller/Wyeomyia/output/ustacks_out/rad_tags.idx.txt -o /home11/mmiller/Wyeomyia/output/radiqual_out/ /home11/mmiller/Wyeomyia/output/velvet_out/velvet-wy_unfiltered_combined_reads_diginorm_paired.fastq.keep_k\=55_e\=19/contigs.fa /home11/mmiller/Wyeomyia/output/velvet_out/velvet-wy_unfiltered_combined_reads_diginorm_paired.fastq.keep_k\=53_e\=19/contigs.fa /home11/mmiller/Wyeomyia/output/velvet_out/velvet-wy_unfiltered_combined_reads_diginorm_paired.fastq.keep_k\=31_e\=22/contigs.fa

Visualization Overview
----------------------

1. Review assembly scores produced by RadiQual.

2. Open IGV.

3. Load previously produced genome assembly fasta file as 'genome.'

        File -> Import Genome

4. Load .bam file produced by RadiQual for loaded assembly.

        *NOTE: Be sure the index (.bam.bai) file is in the same directory as the .bam file loaded.

        File -> Load from File

5. Review the assembly alignments produced by RadiQual for the loaded genome assembly.

        *NOTE: It may be helpful to see all the bases in the aligned track.

        Context Click (on track) -> Show all bases

Visualization Example Walkthrough
---------------------------------

The output is in "/home13/jburkhar/research/Streisfeld/out/align_out/radiqual_out/"

1. read the "assembly_scores.###.txt" file to see how the assemblies compare to each other (the "###" is just because Im lazy)

    a. read the file with less by typing "less assembly_scores.###.txt"
    
    b. press '/' (the forward slash key) to enter search mode
    
    c. type "FILE" and press 'enter' (or 'return' on a mac)
    
    d. press 'n' to go to the next occurrence (these are the different assemblies)
    
    e. hold 'shift' and press 'n' to go to the previous occurrence
    
    f. with this technique you can review every single alignment made against each assembly

2. view the chosen genome assembly

    a. open IGV
    
    b. load the fasta file for the assembly of your choice by navigating to the actual velvet output (this is going to be in "/home13/jburkhar/research/Streisfeld/out/velvet_out")
    
        download it from ACISS and import into IGV:

        File -> Import Genome

    c. load .bam file produced by RadiQual for loaded assembly by finding it in "/home13/jburkhar/research/Streisfeld/out/align_out/radiqual_out/" (youll have to navigate into a similar directory structure as the actual assembly... Im thinking about changing this bc it sucks. The file you want is inside a folder called "contigs.fa" and will be named something like "###_###.merged.sorted.bam". Youll also need the "###_###.merged.sorted.bam.bai" index file.

        download both from ACISS and import the .bam file into IGV: 

        File -> Load from File

3. read the "assembly_aligns.###.txt" file to see which contigs or alignments you might like to view using the same "less" instructions as before.. except searching for "ASSEMLBY" or "CONTIG" or "LOCUS" (IGV allows for individual contig visualization)

        It may be helpful to see all the bases in the aligned track:
    
        context click (right click) on track -> Show all bases

Got questions? Forget it!
    
Kidding, send any feedback to burkhart@cs.uoregon.edu so I can make improvements.
