# TTrichiura_Tubulin
Pipeline to analyze beta-tubulin gene of Trichuris Trichiura, and call nucleotide and amino acid mutations on coding and non-coding regions, separately


The Analysis is based on a reference sequence and the corresponding GFF3 annotation.

Performs:
 - Quality Control and Filter (Trimmomatic)
 - Alignment over reference sequence
 - Full-reference nucleotide consensus calling based on python/samtools and defined criteria
 - Extraction of nucleotide sequence over GFFv2 annotated coding regions
 - Translation of the coding nucleotide sequence
 - Variant calling at the nucleotide level (using viralrecon modules)
 - Amino acid variant calling over coding sequence


 Needs:
 - Trimmomatic
 - Bowtie2
 - samtools
 - gffread
 - seqkit
 