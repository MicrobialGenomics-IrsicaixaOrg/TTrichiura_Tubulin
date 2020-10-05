# ![MicrobialGenomics/TTrichiura_Tubulin](images/nf-core-TTrichiura_Tubulin_logo.png)

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline. Please click [here](https://raw.githack.com/MicrobialGenomics/TTrichiura_Tubulin/master/docs/html/multiqc_report.html) to see an example MultiQC report generated using the parameters defined in [this configuration file](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/blob/master/conf/test_full.config) to run the pipeline on [samples](https://zenodo.org/record/3735111) which were prepared from the [ncov-2019 ARTIC Network V1 amplicon set](https://artic.network/ncov-2019) and sequenced on the Illumina MiSeq platform in 301bp paired-end format.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Preprocessing](#Preprocessing)
    * [parallel-fastq-dump](#parallel-fastq-dump) - Download samples from SRA
    * [cat](#cat) - Merge re-sequenced FastQ files
    * [FastQC](#fastqc) - Raw read QC
    * [trimmomatic](#trimmomatic) - Adapter and quality trimming with trimmomatic
* [Variant calling](#variant-calling)
    * [Bowtie 2](#bowtie-2) - Read alignment relative to reference genome
    * [SAMtools](#samtools) - Sort, index and generate metrics for alignments
    * [iVar trim](#ivar-trim) - Primer sequence removal for amplicon data
    * [picard MarkDuplicates](#picard-markduplicates) - Duplicate read marking and removal
    * [picard CollectMultipleMetrics](#picard-collectmultiplemetrics) - Whole genome coverage and alignment metrics
    * [mosdepth](#mosdepth) - Whole-genome and amplicon coverage metrics
    * [VarScan 2, BCFTools, BEDTools](#varscan-2-bcftools-bedtools) *||* [iVar variants and iVar consensus](#ivar-variants-and-ivar-consensus) *||* [BCFTools and BEDTools](#bcftools-and-bedtools) - Variant calling and consensus sequence generation
        * [SnpEff and SnpSift](#snpeff-and-snpsift) - Genetic variant annotation and functional effect prediction
        * [QUAST](#quast) - Consensus assessment report
    * [BCFTools isec](#bcftools-isec) - Intersect variants across all callers
    * [Custom Consensus](#custom-consensus) - Consensus sequence generation with oun pipeline
    * [Codon frquency](#codfreq) - Codon frequency usage per samble

## Preprocessing

### parallel-fastq-dump

Please see the [usage docs](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/blob/master/docs/usage.md#supported-public-repository-ids) for a list of supported public repository identifiers and how to provide them to the pipeline. The final sample information for all identifiers is obtained from the ENA which provides direct download links for FastQ files as well as their associated md5sums. If a download link exists, the files will be downloaded by FTP otherwise they will be downloaded using [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump).

<details markdown="1">
  <summary>Output files</summary>

* `preprocess/sra/`
    * `sra_run_info.tsv`: Run information file for all samples to be downloaded from the ENA/SRA.
    * `*.fastq.gz`: Paired-end/single-end reads downloaded and extracted from the ENA/SRA.
* `preprocess/sra/md5/`
    * `*.md5`: Files containing `md5` sum for FastQ files downloaded from ENA/SRA.
* `preprocess/sra/log/`
    * `*.fastq_dump.log`: Log file generated from stdout whilst running `parallel-fastq-dump`.

> **NB:** Downloaded FastQ files will only be saved in the results directory if the `--save_sra_fastq` parameter is supplied.  

</details>

### cat

If multiple libraries/runs have been provided for the same sample in the input samplesheet (e.g. to increase sequencing depth) then these will be merged at the very beginning of the pipeline in order to have consistent sample naming throughout the pipeline. Please refer to the [usage docs](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/blob/dev/docs/usage.md#format) to see how to specify these samples in the input samplesheet.

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

<details markdown="1">
  <summary>Output files</summary>

* `preprocess/fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
* `preprocess/fastqc/zips/`
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain adapter sequence and regions of low quality. To see how your reads look after trimming please refer to the FastQC reports in the `preprocess/fastp/fastqc/` directory.

</details>

![MultiQC - FastQC per base sequence plot](images/mqc_fastqc_plot.png)

### trimmomatic

[`Trimmomatic`](http://www.usadellab.org/cms/index.php?page=trimmomatic) a flexible trimmer for Illumina Sequence Data.

<details markdown="1">
  <summary>Output files</summary>

* `preprocess/trimmomatic/`
    * `*.trimmomatic.html`: Trimming report in html format.
    * `*.trim.fastq.gz`: Paired-end/single-end trimmed reads.
    * `*.trim.fail.gz`: Unpaired trimmed reads (only for paired-end data).  
* `preprocess/trimmomatic/log/`
    * `*.trimmomatic.log`: Trimming log file.
* `preprocess/trimmomatic/fastqc/`:
    * `*.trim_fastqc.html`: FastQC report of the trimmed reads.
* `preprocess/trimmomatic/fastqc/zips/`
    * `*.trim_fastqc.zip`: Zip archive containing the FastQC report.

> **NB:** Post-trimmed FastQ files will only be saved in the results directory if the `--save_trimmed` parameter is supplied.

## Variant calling

A file called `summary_variants_metrics_mqc.tsv` containing a selection of read and variant calling metrics will be saved in the `variants/` results directory. The same metrics have also been added to the top of the MultiQC report.

### Bowtie 2

[Bowtie 2](http://bio-bwa.sourceforge.net/) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. Bowtie 2 supports gapped, local, and paired-end alignment modes.

<details markdown="1">
  <summary>Output files</summary>

* `variants/bam/`
    * `<SAMPLE>.bam`: Original BAM file created by Bowtie 2. Only present if `--save_align_intermeds` parameter is supplied.
* `variants/bam/log/`
    * `<SAMPLE>.bowtie2.log`: Bowtie 2 mapping log file.

</details>

![MultiQC - Bowtie2 alignment score plot](images/mqc_bowtie2_plot.png)

### SAMtools

Bowtie 2 BAM files are further processed with [SAMtools](http://samtools.sourceforge.net/) to sort them by coordinate, for indexing, as well as to generate read mapping statistics.

<details markdown="1">
  <summary>Output files</summary>

* `variants/bam/`
    * `<SAMPLE>.sorted.bam`: Coordinate sorted BAM file containing read alignment information.
    * `<SAMPLE>.sorted.bam.bai`: Index file for coordinate sorted BAM file.
* `variants/bam/samtools_stats/`
    * SAMtools `<SAMPLE>.sorted.bam.flagstat`, `<SAMPLE>.sorted.bam.idxstats` and `<SAMPLE>.sorted.bam.stats` files generated from the alignment files.

> **NB:** BAM files and their associated indices will only be saved in the results directory if the `--save_align_intermeds` parameter is supplied.

</details>

![MultiQC - SAMtools alignment scores plot](images/mqc_samtools_stats_plot.png)

### iVar trim

If the `--protocol amplicon` parameter is provided then [iVar](http://gensoft.pasteur.fr/docs/ivar/1.0/manualpage.html) is used to trim amplicon primer sequences from the aligned reads. iVar uses the primer positions supplied in `--amplicon_bed` to soft clip primer sequences from a coordinate sorted BAM file.

<details markdown="1">
  <summary>Output files</summary>

* `variants/bam/`
    * `<SAMPLE>.trim.sorted.bam`: Coordinate sorted BAM file after primer trimming.
    * `<SAMPLE>.trim.sorted.bam.bai`: Index file for coordinate sorted BAM file after primer trimming.
* `variants/bam/samtools_stats/`
    * SAMtools `<SAMPLE>.trim.flagstat`, `<SAMPLE>.trim.idxstats` and `<SAMPLE>.trim.stats` files generated from the primer trimmed alignment files.
* `variants/bam/log/`
    * `<SAMPLE>.trim.ivar.log`: iVar trim log file obtained from stdout.

> **NB:** Post-trimmed BAM files and their associated indices will only be saved in the results directory if the `--save_align_intermeds` parameter is supplied.

</details>

![MultiQC - iVar trim primer heatmap](images/mqc_ivar_trim_plot.png)

### picard MarkDuplicates

Unless you are using [UMIs](https://emea.illumina.com/science/sequencing-method-explorer/kits-and-arrays/umi.html) it is not possible to establish whether the fragments you have sequenced from your sample were derived via true biological duplication (i.e. sequencing independent template fragments) or as a result of PCR biases introduced during the library preparation. By default, the pipeline uses picard MarkDuplicates to *mark* the duplicate reads identified amongst the alignments to allow you to guage the overall level of duplication in your samples. However, you can also choose to remove any reads identified as duplicates via the `--filter_dups` parameter.

<details markdown="1">
  <summary>Output files</summary>

* `variants/bam/`
    * `<SAMPLE>.<SUFFIX>.sorted.bam`: Coordinate sorted BAM file after duplicate marking.
    * `<SAMPLE>.<SUFFIX>.sorted.bam.bai`: Index file for coordinate sorted BAM file after duplicate marking.
* `variants/bam/samtools_stats/`
    * SAMtools `<SAMPLE>.<SUFFIX>.flagstat`, `<SAMPLE>.<SUFFIX>.idxstats` and `<SAMPLE>.<SUFFIX>.stats` files generated from the duplicate marked alignment files.
* `variants/bam/picard_metrics/`
    * `<SAMPLE>.<SUFFIX>.MarkDuplicates.metrics.txt`: Metrics file from MarkDuplicates.

> **NB:** The value of `<SUFFIX>` in the output file names above will depend on the preceeding steps that were run in the pipeline. If `--protocol amplicon` is specified then this process will be run on the iVar trimmed alignments and the value of `<SUFFIX>` will be `trim.mkD`. However, if `--protocol metagenomic` is specified then the process will be run on the alignments obtained directly from Bowtie 2 and the value of `<SUFFIX>` will be `mkD`; where `mkD` is an abbreviation for MarkDuplicates.

</details>

![MultiQC - Picard MarkDuplicates metrics plot](images/mqc_picard_duplicates_plot.png)

### picard CollectMultipleMetrics

[picard-tools](https://broadinstitute.github.io/picard/command-line-overview.html) is a set of command-line tools for manipulating high-throughput sequencing data. We use picard-tools in this pipeline to obtain mapping and coverage metrics.

<details markdown="1">
  <summary>Output files</summary>

* `variants/bam/picard_metrics/`  
    * `<SAMPLE>.<SUFFIX>.CollectMultipleMetrics.*`: Alignment QC files from picard CollectMultipleMetrics in `*_metrics` textual format and plotted in `*.pdf` format.
    * `<SAMPLE>.<SUFFIX>.CollectWgsMetrics.coverage_metrics`: Coverage metrics file from CollectWgsMetrics.

> **NB:** The value of `<SUFFIX>` in the output file names above will depend on the preceeding steps that were run in the pipeline. If `--protocol amplicon` is specified then this process will be run on the iVar trimmed alignments and the value of `<SUFFIX>` will be `trim.mkD`. However, if `--protocol metagenomic` is specified then the process will be run on the alignments obtained directly from Bowtie 2 and the value of `<SUFFIX>` will be `mkD`; where `mkD` is an abbreviation for MarkDuplicates.

</details>

![MultiQC - Picard whole genome coverage plot](images/mqc_picard_wgs_coverage_plot.png)

![MultiQC - Picard insert size plot](images/mqc_picard_insert_size_plot.png)

### mosdepth

[mosdepth](mosdepth) is a fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing. mosdepth is used in this pipeline to obtain genome-wide coverage values in 200bp windows and for `--protocol amplicon` to obtain amplicon/region-specific coverage metrics. The results are then either rendered in MultiQC (genome-wide coverage) or are plotted using custom `R` scripts.

<details markdown="1">
  <summary>Output files</summary>

* `variants/bam/mosdepth/genome/`
    * `<SAMPLE>.<SUFFIX>.genome.mosdepth.global.dist.txt`: A distribution of proportion of bases covered at or above a given threshold for each chromosome and genome-wide.
    * `<SAMPLE>.<SUFFIX>.genome.mosdepth.region.dist.txt`: A distribution of proportion of bases covered at or above a given threshold for each chromosome and genome-wide.
    * `<SAMPLE>.<SUFFIX>.genome.mosdepth.summary.txt`: Summary metrics including mean, min and max coverage values.
    * `<SAMPLE>.<SUFFIX>.genome.per-base.bed.gz`: Per-base depth output genome-wide.
    * `<SAMPLE>.<SUFFIX>.genome.per-base.bed.gz.csi`: CSI index that can be used for tabix queries from above file.
    * `<SAMPLE>.<SUFFIX>.genome.regions.bed.gz`: Mean regional depth for 200bp windows genome-wide.
    * `<SAMPLE>.<SUFFIX>.genome.regions.bed.gz.csi`: CSI index that can be used for tabix queries from above file.
* `variants/bam/mosdepth/genome/plots/`
    * `all_samples.<SUFFIX>.genome.regions.coverage.tsv`: File aggregating genome-wide coverage values across all samples used for plotting.
    * `<SAMPLE>.<SUFFIX>.genome.regions.coverage.pdf`: Whole-genome coverage plot.
    * `<SAMPLE>.<SUFFIX>.genome.regions.coverage.tsv`: File containing coverage values for the above plot.
* `variants/bam/mosdepth/amplicon/`
    * `<SAMPLE>.<SUFFIX>.amplicon.mosdepth.global.dist.txt`: A distribution of proportion of bases covered at or above a given threshold for each chromosome and genome-wide.
    * `<SAMPLE>.<SUFFIX>.amplicon.mosdepth.region.dist.txt`: A distribution of proportion of bases covered at or above a given threshold for each chromosome and genome-wide.
    * `<SAMPLE>.<SUFFIX>.amplicon.mosdepth.summary.txt`: Summary metrics including mean, min and max coverage values.
    * `<SAMPLE>.<SUFFIX>.amplicon.per-base.bed.gz`: Per-base depth output genome-wide.
    * `<SAMPLE>.<SUFFIX>.amplicon.per-base.bed.gz.csi`: CSI index that can be used for tabix queries from above file.
    * `<SAMPLE>.<SUFFIX>.amplicon.regions.bed.gz`: Mean regional depth for individual amplicons genome-wide.
    * `<SAMPLE>.<SUFFIX>.amplicon.regions.bed.gz.csi`: CSI index that can be used for tabix queries from above file.
    * `<SAMPLE>.<SUFFIX>.amplicon.thresholds.bed.gz`: Threshold output to indicate how many bases in each region are covered at given thresholds.
    * `<SAMPLE>.<SUFFIX>.amplicon.thresholds.bed.gz.csi`: CSI index that can be used for tabix queries from above file.
* `variants/bam/mosdepth/amplicon/plots/`
    * `all_samples.<SUFFIX>.amplicon.regions.coverage.tsv`: File aggregating per-amplicon coverage values across all samples used for plotting.
    * `all_samples.<SUFFIX>.amplicon.regions.heatmap.pdf`: Heatmap showing per-amplicon coverage across all samples.
    * `<SAMPLE>.<SUFFIX>.amplicon.regions.coverage.pdf`: Bar plot showing per-amplicon coverage for an individual sample.
    * `<SAMPLE>.<SUFFIX>.amplicon.regions.coverage.tsv`: File containing per-amplicon coverage values for the above plot.

> NB: The value of `<SUFFIX>` in the output file names above will depend on the preceeding steps that were run in the pipeline. If `--protocol amplicon` is specified then this process will be run on the iVar trimmed alignments and the value of `<SUFFIX>` will be `trim.mkD`. However, if `--protocol metagenomic` is specified then the process will be run on the alignments obtained directly from Bowtie 2 and the value of `<SUFFIX>` will be `mkD`; where `mkD` is an abbreviation for MarkDuplicates.

</details>

![R - Sample genome-wide coverage plot](images/r_genome_coverage.png)

<p markdown="1" align="center">
  <img src="images/r_amplicon_barplot.png" alt="R - Sample per-amplicon coverage plot">
</p>

### VarScan 2, BCFTools, BEDTools

[VarScan 2](http://dkoboldt.github.io/varscan/) is a platform-independent software tool to detect variants in NGS data. In this pipeline, VarScan 2 is used in conjunction with SAMtools in order to call both high and low frequency variants.

[BCFtools](http://samtools.github.io/bcftools/bcftools.html) is a set of utilities that manipulate variant calls in [VCF](https://vcftools.github.io/specs.html) and its binary counterpart BCF format. BCFTools is used in the variant calling and *de novo* assembly steps of this pipeline to obtain basic statistics from the VCF output. It is also used in the VarScan 2 variant calling branch of the pipeline to generate a consensus sequence by integrating high frequency variant calls into the reference genome.

[BEDTools](https://bedtools.readthedocs.io/en/latest/) is a swiss-army knife of tools for a wide-range of genomics analysis tasks. In this pipeline we use `bedtools genomecov` to compute the per-base mapped read coverage in bedGraph format, and `bedtools maskfasta` to mask sequences in a Fasta file based on intervals defined in a feature file. This may be useful for creating your own masked genome file based on custom annotations or for masking all but your target regions when aligning sequence data from a targeted capture experiment.

<details markdown="1">
  <summary>Output files</summary>

* `variants/varscan2/`
    * `<SAMPLE>.vcf.gz`: Low frequency variants VCF file.
    * `<SAMPLE>.vcf.gz.tbi`: Low frequency variants VCF index file.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.vcf.gz`: High frequency variants VCF file.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.vcf.gz.tbi`: High frequency variants VCF index file.
* `variants/varscan2/consensus/`
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.consensus.fa`: Consensus Fasta file generated by integrating the high frequency variants called by VarScan into the reference genome.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.consensus.masked.fa`: Masked consensus Fasta file.
* `variants/varscan2/consensus/base_qc/`
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.ACTG_density.pdf`: Plot showing density of ACGT bases within the consensus sequence.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.base_counts.pdf`: Plot showing frequency and percentages of all bases in consensus sequence.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.base_counts.tsv`: File containing frequency and percentages of all bases in consensus sequence.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.N_density.pdf`: Plot showing density of N bases within the consensus sequence.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.N_run.tsv`: File containing start positions and width of N bases in consensus sequence.
* `variants/varscan2/bcftools_stats/`
    * `<SAMPLE>.bcftools_stats.txt`: Statistics and counts obtained from low frequency variants VCF file.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.bcftools_stats.txt`: Statistics and counts obtained from high frequency variants VCF file.
* `variants/varscan2/log/`
    * `<SAMPLE>.varscan2.log`: Log file generated from stderr by VarScan 2.
* `variants/bam/mpileup/`
    * `<SAMPLE>.<SUFFIX>.mpileup`: mpileup files summarize all the data from aligned reads at a given genomic position. Each row of the mpileup file gives similar information to a single vertical column of reads as visualised in IGV.

> **NB:** The value of `<MAX_ALLELE_FREQ>` in the output file names above is determined by the `--max_allele_freq` parameter (Default: 0.8).  
> **NB:** Output mpileup files will only be saved in the  directory if the `--save_mpileup` parameter is supplied. The naming convention for these files will depend on the preceeding steps that were run in the pipeline as described in the paragraph explaining the value of `<SUFFIX>` in the section above.

</details>

![MultiQC - VarScan 2 variants called plot](images/mqc_varscan2_plot.png)

### iVar variants and iVar consensus

[iVar](https://github.com/andersen-lab/ivar/blob/master/docs/MANUAL.md) is a computational package that contains functions broadly useful for viral amplicon-based sequencing. We use iVar in this pipeline to [trim primer sequences](#ivar-trim) for amplicon input data as well as to call variants and for consensus sequence generation.

<details markdown="1">
  <summary>Output files</summary>

* `variants/ivar/`
    * `<SAMPLE>.tsv`: Low frequency variants in TSV format.
    * `<SAMPLE>.vcf.gz`: Low frequency variants VCF file.
    * `<SAMPLE>.vcf.gz.tbi`: Low frequency variants VCF index file.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.vcf.gz`: High frequency variants VCF file.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.vcf.gz.tbi`: High frequency variants VCF index file.
* `variants/ivar/consensus/`
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.consensus.fa`: Consensus Fasta file generated by iVar at the frequency threshold set by the `--max_allele_freq` parameter.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.consensus.qual.txt`: File with the average quality of each base in the consensus sequence.
* `variants/ivar/consensus/base_qc/`
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.ACTG_density.pdf`: Plot showing density of ACGT bases within the consensus sequence.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.base_counts.pdf`: Plot showing frequency and percentages of all bases in consensus sequence.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.base_counts.tsv`: File containing frequency and percentages of all bases in consensus sequence.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.N_density.pdf`: Plot showing density of N bases within the consensus sequence.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.N_run.tsv`: File containing start positions and width of N bases in consensus sequence.
* `variants/ivar/log/`
    * `<SAMPLE>.variant.counts.log`: Variant counts for low frequency variants.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.variant.counts.log`: Variant counts for high frequency variants.
* `variants/ivar/bcftools_stats/`
    * `<SAMPLE>.bcftools_stats.txt`: Statistics and counts obtained from low frequency variants VCF file.
    * `<SAMPLE>.AF<MAX_ALLELE_FREQ>.bcftools_stats.txt`: Statistics and counts obtained from high frequency variants VCF file.

</details>

![MultiQC - iVar variants called plot](images/mqc_ivar_variants_plot.png)

### BCFTools and BEDTools

[BCFtools](http://samtools.github.io/bcftools/bcftools.html) can be used to call variants directly from BAM alignment files. The functionality to call variants with BCFTools in this pipeline was inspired by work carried out by [Conor Walker](https://github.com/conorwalker/covid19/blob/3cb26ec399417bedb7e60487415c78a405f517d6/scripts/call_variants.sh). In contrast to VarScan 2 and iVar, the original variant calls obtained by BCFTools are not filtered further by a higher allele frequency. It seems that the default calls obtained by BCFTools appear to be comparable with the high frequency variants generated by VarScan 2 and iVar.

<details markdown="1">
  <summary>Output files</summary>

* `variants/bcftools/`
    * `<SAMPLE>.vcf.gz`: Variants VCF file.
    * `<SAMPLE>.vcf.gz.tbi`: Variants VCF index file.
* `variants/bcftools/consensus/`
    * `<SAMPLE>.consensus.fa`: Consensus Fasta file generated by integrating the variants called by BCFTools into the reference genome.
    * `<SAMPLE>.consensus.masked.fa`: Masked consensus Fasta file.
* `variants/bcftools/consensus/base_qc/`
    * `<SAMPLE>.ACTG_density.pdf`: Plot showing density of ACGT bases within the consensus sequence.
    * `<SAMPLE>.base_counts.pdf`: Plot showing frequency and percentages of all bases in consensus sequence.
    * `<SAMPLE>.base_counts.tsv`: File containing frequency and percentages of all bases in consensus sequence.
    * `<SAMPLE>.N_density.pdf`: Plot showing density of N bases within the consensus sequence.
    * `<SAMPLE>.N_run.tsv`: File containing start positions and width of N bases in consensus sequence.
* `variants/bcftools/bcftools_stats/`
    * `<SAMPLE>.bcftools_stats.txt`: Statistics and counts obtained from VCF file.

</details>

![MultiQC - BCFTools variant counts](images/mqc_bcftools_stats_plot.png)

### SnpEff and SnpSift

[SnpEff](http://snpeff.sourceforge.net/SnpEff.html) is a genetic variant annotation and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

[SnpSift](http://snpeff.sourceforge.net/SnpSift.html) annotates genomic variants using databases, filters, and manipulates genomic annotated variants. After annotation with SnpEff, you can use SnpSift to help filter large genomic datasets in order to find the most significant variants.

<details markdown="1">
  <summary>Output files</summary>

* `variants/<CALLER>/snpeff/`
    * `*.snpEff.csv`: Variant annotation csv file.
    * `*.snpEff.genes.txt`: Gene table for annotated variants.
    * `*.snpEff.summary.html`: Summary html file for variants.
    * `*.snpEff.vcf.gz`: VCF file with variant annotations.
    * `*.snpEff.vcf.gz.tbi`: Index for VCF file with variant annotations.
    * `*.snpSift.table.txt`: SnpSift summary table.

> **NB:** The value of `<CALLER>` in the output directory name above is determined by the `--callers` parameter (Default: 'varscan2,ivar,bcftools'). If applicable, you will have two sets of files where the file name prefix will be `<SAMPLE>` for low-frequency variants and `<SAMPLE>.AF<MAX_ALLELE_FREQ>` for high frequency variants.

</details>

![MultiQC - SnpEff annotation counts](images/mqc_snpeff_plot.png)

### QUAST

[QUAST](http://bioinf.spbau.ru/quast) is used to generate a single report with which to evaluate the quality of the consensus sequence across all of the samples provided to the pipeline. The HTML results can be opened within any browser (we recommend using Google Chrome). Please see the [QUAST output docs](http://quast.sourceforge.net/docs/manual.html#sec3) for more detailed information regarding the output files.

<details markdown="1">
  <summary>Output files</summary>

* `variants/<CALLER>/quast/AF<MAX_ALLELE_FREQ>/`
    * `report.html`: Results report in HTML format. Also available in various other file formats i.e. `report.pdf`, `report.tex`, `report.tsv` and `report.txt`.

> **NB:** The value of `<CALLER>` in the output directory name above is determined by the `--callers` parameter (Default: 'varscan2,ivar,bcftools') and the value of `<MAX_ALLELE_FREQ>` is determined by the `--max_allele_freq` parameter (Default: 0.8).

</details>

### BCFTools isec

[BCFTools isec](http://samtools.github.io/bcftools/bcftools.html#isec) can be used to intersect the variant calls generated by the 3 different callers used in the pipeline. This permits a quick assessment of how consistently a particular variant is being called using different algorithms and to prioritise the investigation of the variants.

<details markdown="1">
  <summary>Output files</summary>

* `variants/intersect/<SAMPLE>/`
    * `*.vcf.gz`: VCF file containing variants common to at least 2/3 callers. There will be one file for each caller - see `README.txt` for details.
    * `*.vcf.gz.tbi`: Index for VCF file.
    * `README.txt`: File containing command used and file name mappings.
    * `sites.txt`: List of variants common to at least 2/3 callers in textual format. The last column indicates presence (1) or absence (0) amongst the 3 different callers.

> **NB:** This process will only be executed when all 3 variant callers are specified to run, as is by default i.e. `--callers varscan2,ivar,bcftools`.

</details>


## Workflow reporting and genomes

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from FastQC, fastp, Cutadapt, Bowtie 2, Kraken 2, VarScan 2, iVar, samtools flagstat, samtools idxstats, samtools stats, picard CollectMultipleMetrics and CollectWgsMetrics, BCFTools, SnpEff and QUAST.

The default [`multiqc config file`](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/blob/master/assets/multiqc_config.yaml) has been written in a way in which to structure these QC metrics to make them more interpretable in the final report.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

Please click [here](https://raw.githack.com/MicrobialGenomics/TTrichiura_Tubulin/master/docs/html/multiqc_report.html) to see an example MultiQC report generated using the parameters defined in [this configuration file](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/blob/master/conf/test_full.config) to run the pipeline on [samples](https://zenodo.org/record/3735111) which were prepared from the [ncov-2019 ARTIC Network V1 amplicon set](https://artic.network/ncov-2019) and sequenced on the Illumina MiSeq platform in 301bp paired-end format.

<details markdown="1">
  <summary>Output files</summary>

* `multiqc/`  
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

### Reference genome files

A number of genome-specific files are generated by the pipeline because they are required for the downstream processing of the results. If the `--save_reference` parameter is provided then the Bowtie 2 alignment indices, BLAST and Kraken 2 databases downloaded/generated by the pipeline will be saved in the `genome/` directory. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build a Kraken 2 database for the host genome. This can be quite a time-consuming process and it permits their reuse for future runs of the pipeline or for other purposes.

<details markdown="1">
  <summary>Output files</summary>

* `genome/`  
    * `BlastDB/`: BLAST database for viral genome.
    * `Bowtie2Index/`: Bowtie 2 index for viral genome.
    * `kraken2_<KRAKEN2_DB_NAME>/`: Kraken 2 database for host genome.
    * `SnpEffDB/`: SnpEff database for viral genome.
    * `snpeff.config`: SnpEff config file for viral genome.
    * Unzipped genome fasta file for viral genome
    * Unzipped genome annotation GFF file for viral genome

</details>

### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">
  <summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
    * Documentation for interpretation of results in HTML format: `results_description.html`.

</details>
