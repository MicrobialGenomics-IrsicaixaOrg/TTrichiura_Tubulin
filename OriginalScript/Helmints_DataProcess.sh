#!/bin/bash
# MNJ, Sep 2020.

ReferenceFile=/Users/mnoguera/Documents/Work/Projects/Helmints_2020/Sequencing/Reference/AF034219.fasta
#### We have done whatever possible to trim adapters, final control for non-assembled reads,
for file in *_R1_*fastq.gz
do
    name=${file%%_S*}
#### Trim low quality Reads & adapters (Nextera in this case)
    trimmomatic PE $file ${file%%_R1_001.fastq.gz}_R2_001.fastq.gz -summary prova.txt -baseout $name LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50
### Aligns against reference
    bowtie2 -p 4 -x $ReferenceFile -1 ${name}_1P -2 ${name}_2P -S ${name}.sam --very-sensitive-local
### samtools sort and index
    ~/bin/own_samSortAndIndex.sh ${name}.sam
### Calculate mean depth of coverage
	samtools depth ${name}.sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${name}_depthOfCoverage.txt
### Extract consensus from file
	python own_consensusSequence_v2.py ${name}.sorted.bam 15 > ${name}_consensus.fasta


#samtools mpileup -q 20 -Q 30 -uf $ReferenceFile ${name}.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq -Q 20 -d 20 > ${name}_consensus.fastq
#sed -n '1~4s/^@/>/p;2~4p' ${name}_consensus.fastq > ${name}_consensus.fasta
#seqkit fq2fa ${name}_consensus.fastq > ${name}_consensus.fasta # produces a multi-line fasta which is not ok for gffread
#seqtk seq -l 0 ${name}_consensus.fasta > tmp
#	mv tmp ${name}_consensus.fasta


### Need to rename the consensus in order to use GFF3 annotation (need of same name)
	echo ">AF034219.1" > tmp ### Use Reference accession to process GFF3
	fgrep -v ">" ${name}_consensus.fasta >> tmp
	mv tmp ${name}_consensus.fasta
### rm previous fasta index if any, it may interfere
	rm ${name}_consensus.fasta.fai

### Extract coding region of consensus
	gffread -W -x ${name}_consensus_cds.fasta  -g ${name}_consensus.fasta ${ReferenceFile%%.fasta}.gff3

### Recover original sample name
	echo ">${name}" > tmp ### Recover Sample Name
    fgrep -v ">" ${name}_consensus.fasta >> tmp
    mv tmp ${name}_consensus.fasta
### Translate coding region
	seqkit translate --quiet ${name}_consensus_cds.fasta -o ${name}_consensus_cds.faa
### Clean-up
	rm ${name}_1P ${name}_2P ${name}_1U ${name}_2U ${name}.sam ${name}.bam ${name}_consensus.fasta.fai
done
