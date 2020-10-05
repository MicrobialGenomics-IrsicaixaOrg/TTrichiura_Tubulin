# MicrobialGenomics/TTrichiura_Tubulin

**Pipeline to analyze beta-tubulin gene of Trichuris Trichiura, and call nucleotide and amino acid mutations on coding and non-coding regions, separately**.
Performs:
- Quality Control and Filter (Trimmomatic)
- Alignment over reference sequence
- Full-reference nucleotide consensus calling based on python/samtools and defined criteria
- Extraction of nucleotide sequence over GFFv2 annotated coding regions
- Translation of the coding nucleotide sequence
- Variant calling at the nucleotide level (using viralrecon modules)
- Amino acid variant calling over coding sequence

[![GitHub Actions CI Status](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/workflows/nf-core%20CI%20amplicon/badge.svg)](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/actions)
[![Actions Status](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/workflows/nf-core%20linting/badge.svg)](https://github.com/MicrobialGenomics/TTrichiura_Tubulin/actions)

[![Docker](https://img.shields.io/docker/automated/microbialgenomics/TTrichiura_Tubulin.svg)](https://hub.docker.com/r/microbialgenomics/TTrichiura_Tubulin)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

## Pipeline additions

1. Adapter trimming with [`Trimmomatic`](http://www.usadellab.org/cms/index.php?page=trimmomatic)
2. Codon Frequency with [`CodFrq`](https://github.com/hivdb/codfreq)
3. Custom consensus with [`consensusSequence_v2.py`](bin/consensusSequence_v2.py)

## Introduction

**nfcore/TTrichiura_Tubulin** is a bioinformatics analysis pipeline used to perform low-frequency variant calling for beta tubulin gene of ttrichura. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with Docker containers making installation trivial and results highly reproducible. Furthermore, automated continuous integration tests that run the pipeline on a full-sized dataset using AWS cloud ensure that the code is stable.

## Pipeline summary

1. Download samples via SRA, ENA or GEO ids ([`ENA FTP`](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html), [`parallel-fastq-dump`](https://github.com/rvalieris/parallel-fastq-dump); *if required*)
2. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html); *if required*)
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Adapter trimming ([`fastp`](https://github.com/OpenGene/fastp) or [`trimmomatic`](http://www.usadellab.org/cms/index.php?page=trimmomatic))
5. Variant calling
    1. Read alignment ([`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
    2. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3. Primer sequence removal ([`iVar`](https://github.com/andersen-lab/ivar); *amplicon data only*)
    4. Duplicate read marking ([`picard`](https://broadinstitute.github.io/picard/); *removal optional*)
    5. Alignment-level QC ([`picard`](https://broadinstitute.github.io/picard/), [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    6. Genome-wide and amplicon coverage QC plots ([`mosdepth`](https://github.com/brentp/mosdepth/))
    7. Choice of multiple variant calling and consensus sequence generation routes ([`VarScan 2`](http://dkoboldt.github.io/varscan/), [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/) *||* [`iVar variants and consensus`](https://github.com/andersen-lab/ivar) *||* [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/))
        * Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
        * Consensus assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
    8. Intersect variants across callers ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html))
    9. Custom consensus ([`consensusSequence_v2.py`](bin/consensusSequence_v2.py))
    10. Codon frequency calling ([`codfrq`](https://github.com/hivdb/codfreq))
6. Present QC and visualisation for raw read, alignment, assembly and variant calling results ([`MultiQC`](http://multiqc.info/))

> **NB:** The pipeline has a number of options to allow you to run only specific aspects of the workflow if you so wish.
For example, you can skip all of the assembly steps with the `--skip_assembly` parameter.
See the [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Pipeline reporting

Numerous QC and reporting steps are included in the pipeline in order to collate a full summary of the analysis within a single [MultiQC](https://multiqc.info/) report. You can see [an example MultiQC report here](https://raw.githack.com/MicrobialGenomics/viralrecon/master/docs/html/multiqc_report.html), generated using the parameters defined in [this configuration file](https://github.com/nf-core/viralrecon/blob/master/conf/test_full.config). The pipeline was run with [these samples](https://zenodo.org/record/3735111), prepared from the [ncov-2019 ARTIC Network V1 amplicon set](https://artic.network/ncov-2019) and sequenced on the Illumina MiSeq platform in 301bp paired-end format.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run MicrobialGenomics/TTrichiura_Tubulin -profile test,<docker/conda>
    ```

4. Start running your own analysis!

    * Typical command for shotgun analysis:

        ```bash
        nextflow run  MicrobialGenomics/TTrichiura_Tubulin \
            --input samplesheet.csv \
            -profile <docker/conda>
        ```

    * Typical command for amplicon analysis:

        ```bash
        nextflow run  MicrobialGenomics/TTrichiura_Tubulin \
            --input samplesheet.csv \
            --protocol amplicon \
            --amplicon_bed ./nCoV-2019.artic.V3.bed \
            --skip_assembly \
            -profile <docker/conda>
        ```

See the [usage documentation](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The  MicrobialGenomics/TTrichiura_Tubulin pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](docs/usage.md#reference-genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)
