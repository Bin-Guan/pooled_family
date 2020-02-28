#!/bin/bash

#Continued from mutect vcf.gz
#sbatch --mem=16g --cpus-per-task=1 ~/git/pooled_family/mutect_annovar_intervar.sh VCF
set -e

VCF=$1 #vcf.gz

#vt
module load vt/0.577

zcat $VCF \
| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
| vt decompose -s - \
| vt normalize -r /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - \
> ${VCF%.vcf.gz}.vt.vcf
	
module load intervar/2.1.2

convert2annovar.pl -format vcf4old ${VCF%.vcf.gz}.vt.vcf -includeinfo --outfile ${VCF%.vcf.gz}.avinput

table_annovar.pl ${VCF%.vcf.gz}.avinput $ANNOVAR_DATA/hg19 -buildver hg19 -remove \
	-out ${VCF%.vcf.gz} \
	--protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,dbnsfp33a,clinvar_20190305,gnomad_genome,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene,refGeneWithVer,popfreq_max_20150413,gnomad_exome,spidex,avsnp150,spliceai_filtered \
	-operation g,f,f,f,f,f,f,f,f,r,g,g,g,f,f,f,f,f \
	--argument '-hgvs',,,,,,,,,,,,'-splicing 50 -hgvs',,,,, --polish -nastring . --thread $SLURM_CPUS_PER_TASK --otherinfo

InterVar -i ${VCF%.vcf.gz}.avinput \
	--input_type=AVinput -d $ANNOVAR_DATA/hg19 \
	-o ${VCF%.vcf.gz} --skip_annovar --evidence_file=/data/OGVFB/OGL_NGS/OGL.variants.evidence.txt
	

sed "1 s/"Otherinfo"/"VCF_CHROM\\tVCF_POS\\tVCF_ID\\tVCF_REF\\tVCF_ALT\\tVCF_QUAL\\tVCF_FILTER\\tVCF_INFO\\tVCF_FORMAT\\tVCF_TUMOR\\tVCF_NORMAL"/" ${VCF%.vcf.gz}.hg19_multianno.txt > ${VCF%.vcf.gz}.annovar.txt

module load R/3.5.2
# #$1 = intervar ouput, $2 = annovar output, $3 = intervar_annovar_output for vcfanno (must be txt file), $4 OGLv1_panel_DxORcandidate

Rscript ~/git/pooled_family/annovar_intervar_mutect_v1.6.R ${VCF%.vcf.gz}.hg19_multianno.txt.intervar ${VCF%.vcf.gz}.annovar.txt ${VCF%.vcf.gz}.annovar.intervar.filtered.txt /data/OGVFB/OGL_NGS/variant_prioritization/data/OGLv1_panel_DxORcandidate.tsv



