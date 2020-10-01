# ![nf-core/tt](docs/images/nf-core-tt_logo.png)

**Pipeline to analyze beta-tubulin gene of Trichuris Trichiura, and call nucleotide and amino acid mutations on coding and non-coding regions, separately**.

[![GitHub Actions CI Status](https://github.com/nf-core/tt/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/tt/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/tt/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/tt/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/tt.svg)](https://hub.docker.com/r/nfcore/tt)

## Introduction

MicrobialGenomics/TTrichiura_Tubulin is a bioinformatic pipeline used to perform intra-host/low-frequency variant calling at genomic and protein level for beta tubulin gene of Trichuris trichuria. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline summary

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html); *if required*)
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter trimming [`trimmomatic`](http://www.usadellab.org/cms/index.php?page=trimmomatic))
4. Variant calling
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
    9. Custom consensus ([consensusSequence_v2.py](bin/consensusSequence_v2.py))
    9. Codon frequency calling ([`codfrq`](https://github.com/hivdb/codfreq))
5. Present QC and visualisation for raw read, alignment, assembly and variant calling results ([`MultiQC`](http://multiqc.info/))

> **NB:** The pipeline has a number of options to allow you to run only specific aspects of the workflow if you so wish.
See the [usage docs](docs/usage.md) for all of the available options when running the pipeline.

