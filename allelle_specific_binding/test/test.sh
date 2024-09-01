#!/bin/bash

module load bzip2

SAMPLE="3782-BW-64"
BAM="/mnt/isilon/sfgi/dudekm/raw_data/brandon_liver_ATAC/chr_renamed_bams/${SAMPLE}.bam"

BEDTOOLS="/home/dudekmf/local/bin/bedtools/bedtools"

# python ../splitSNP.py $BAM out/${SAMPLE}_rs12774423 chr10:718264:G:A

$BEDTOOLS bamtobed -i bams/unpaired_test.bam -bedpe

echo ""

$BEDTOOLS bamtobed -i bams/paired_test.bam -bedpe
