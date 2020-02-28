#!/bin/bash


module load GATK/3.8-0

fam=$1
#sbatch --cpus-per-task=16 --mem=64g --time=72:00:00 mutect.sh 2053
#https://software.broadinstitute.org/gatk/documentation/article?id=9183
#some files need ~45 hours
#July 2019: 2 exome samples finished in < 14 hours and the 3rd finished <18 hours.

GATK -m 64g MuTect2 \
     -threads 16 \
     -I:tumor ../sample_bam/$1aff.recalibrated.bam \
     -I:normal ../sample_bam/$1unaff.recalibrated.bam \
     --output_mode EMIT_VARIANTS_ONLY \
     -o $1_mutect2.vcf.gz \
	 --interval_padding 200 \
     -L /data/OGVFB/OGL_NGS/bed/xgen-exome-research-panel-targets-nochr-sorted.bed \
     -R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta  