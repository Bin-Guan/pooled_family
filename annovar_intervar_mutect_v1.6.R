### sort intervar output according to PVS/PS/PM/PP ###
### after InterVar and Annovar to prepare for vcfanno
## PVS = 8
## PS  = 6
## PM  = 3
## PP  = 1
## BA = -5
## BS = -3
## BP = -1
## if PVS = 0 & dbscSNV's ada > 0.8 and rf>0.5, Score=+3 (equals weight of PM)
## if PVS = 0 & frameshit or stop_gain or Start_loss (?intervar), Score=+3 (equals weight of PM)
## if PVS = 0 & |dpsi_max_tissue+dpsi_zscore| > 5, score=+3; > 2.5, score=+1 make a histogram of dpsi of WGS set and determine the cut-off
## if PVS = 0 & spliceai_rank >0.8, score=+8; >0.5, score=+6; >0.2, score=+3; >0.15, score=+1; make a histogram of splicai score and determine the cut-off
## Intervar - select one gene for each variant

args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)

# change the file name below
intervar <- read.delim(args[1], sep = "\t", header = TRUE, na.strings = c("."),
                       colClasses = c("factor","integer","integer","character","character","character","character","character","character",
										"character","character","character","character","character","numeric","numeric","numeric","numeric",
										"numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","numeric",
										"character","integer","character","character","character","character") )

#does the read.delim evaluate all rows for the types? seems to be yes.
annovar <- read.delim(args[2], sep = "\t", header = TRUE, na.strings = c("."),
                      colClasses = c("factor","integer","integer","character","character","character","character","character","character","character",
                                     "numeric","numeric","character","numeric","numeric","character","numeric","numeric","character","numeric",
                                     "numeric","character","numeric","numeric","character","numeric","numeric","character","numeric","numeric",
                                     "character","numeric","numeric","character","numeric","numeric","character","numeric","numeric","numeric",
                                     "numeric","character","numeric","numeric","character","numeric","numeric","character","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","numeric","character","character","numeric","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","character",
                                     "character","character","character","character","numeric","numeric","numeric","numeric","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","character","character","character","character","character","character",
                                     "character","character","character","character","character","character","character","character","character","character",
                                     "character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                     "numeric","numeric","numeric","character","character","factor","integer","character","character","character",
                                     "numeric","character","character","character","character","character") ) %>%
  separate('VCF_TUMOR', c('TUMOR_GT','TUMOR_AD','TUMOR_AF','TUMOR_4','TUMOR_5','TUMOR_6','TUMOR_7','TUMOR_8','TUMOR_9'), sep = ':', remove = FALSE, convert = TRUE) %>% 
  separate('VCF_NORMAL', c('NORMAL_GT','NORMAL_AD','NORMAL_AF','NORMAL_4','NORMAL_5','NORMAL_6','NORMAL_7','NORMAL_8','NORMAL_9'), sep = ':', remove = FALSE, convert = TRUE)

OGLv1_gene_class <- read.delim(args[4], sep = "\t", header = T, colClasses = c("character","character", "character") ) 

annovar_OGLv1_class <- left_join(annovar, OGLv1_gene_class, by = c("Gene.refGene" = "gene"))

#annovar_OGLv1_class <- merge(x = annovar, y = OGLv1_gene_class, 
#                             by.x = c("Ref.Gene"), by.y = c("gene"), all.x = TRUE, 
#                            sort = FALSE, suffixes = c(".annovar", ".intervar"), no.dups = TRUE,
#                             incomparables = NULL) 

#pick one annotation for each variant that is of highest Priority score
intervar_for_sorting <- intervar %>% 
  separate('InterVar..InterVar.and.Evidence', 
           c('InterVarInterpretation','PVS1','PS','PM','PP','BA1','BS','BP'), 
           sep = '\\]{0,1} [A-Z]{2,3}\\d{0,1}\\=\\[{0,1}', remove = FALSE, convert = TRUE) %>% 
  separate(PS, c('PS1','PS2','PS3','PS4','PS5'), sep = ',', convert = TRUE) %>% 
  separate(PM, c('PM1','PM2','PM3','PM4','PM5','PM6','PM7'), sep =',', convert = TRUE) %>% 
  separate(PP, c('PP1','PP2','PP3','PP4','PP5','PP6'), sep = ',', convert = TRUE) %>% 
  separate(BS, c('BS1','BS2','BS3','BS4','BS5'), sep = ',', convert = TRUE) %>% 
  separate(BP, c('BP1','BP2','BP3','BP4','BP5','BP6','BP7',"BP8"), sep = ',', convert = TRUE) %>% 
  mutate(BP8 = str_sub(BP8, 2, 2)) %>% 
  mutate(BP8 = as.integer(BP8)) %>%
  mutate(dbscSNV_ADA_SCORE_temp = dbscSNV_ADA_SCORE) %>% 
  mutate(dbscSNV_RF_SCORE_temp = dbscSNV_RF_SCORE) %>% 
  replace_na(list(dbscSNV_ADA_SCORE_temp = 0, dbscSNV_RF_SCORE_temp = 0)) %>%
  mutate(dbscSNV_value = ifelse((PVS1 == 0 & dbscSNV_ADA_SCORE_temp>0.8 & dbscSNV_RF_SCORE_temp>0.5), 1, 0)) %>% 
  mutate(truncating = ifelse((PVS1 == 0 & grepl("^frameshift|stop", ExonicFunc.refGene, ignore.case = TRUE)), 1, 0)) %>% 
  mutate(Priority.Score = (PVS1*8+(PS1+PS2+PS3+PS4+PS5)*6+(PM1+PM2+PM3+PM4+PM5+PM6+PM7+dbscSNV_value+truncating)*3+(PP1+PP2+PP3+PP4+PP5+PP6)-BA1*5-(BS1+BS2+BS3+BS4+BS5)*3-(BP1+BP2+BP3+BP4+BP5+BP6+BP7+BP8))) %>% 
  unite("variantkey", X.Chr:Alt, sep = "_", remove = FALSE ) %>%
  group_by(variantkey) %>%
  slice(which.max(Priority.Score))
  
annovar_inter <- merge(x = annovar_OGLv1_class, y = intervar_for_sorting, 
                       by.x = c("Chr", "Start", "End", "Ref", "Alt"), by.y = c("X.Chr", "Start", "End", "Ref", "Alt"), all.x = TRUE, 
                       sort = FALSE, suffixes = c(".annovar", ".intervar"), no.dups = TRUE,
                       incomparables = NULL) %>%
  mutate(dpsi_max_tissue_temp = dpsi_max_tissue) %>% 
  mutate(dpsi_zscore_temp = dpsi_zscore) %>% 
  replace_na(list(dpsi_max_tissue_temp = 0, dpsi_zscore_temp = 0 )) %>%
  mutate(dpsi_value = ifelse(PVS1 == 0 & abs(dpsi_max_tissue_temp + dpsi_zscore_temp) > 5, 3, (ifelse((PVS1 == 0 & abs(dpsi_max_tissue_temp + dpsi_zscore_temp) > 2.5), 1, 0)))  ) %>%
  mutate(Priority.Score = Priority.Score + dpsi_value ) %>% 
  mutate(PopFreqMax = if_else(is.na(PopFreqMax), 0, PopFreqMax)) %>%
  mutate(gnomAD_genome_ALL = if_else(is.na(gnomAD_genome_ALL), 0, gnomAD_genome_ALL)) %>% 
  filter(PopFreqMax < 0.05)  %>% 
  filter(gnomAD_genome_ALL < 0.05) %>% 
  filter(Priority.Score > -3) %>% 
  filter(TUMOR_AF > 0.10) %>% 
  select(c("Chr", "Start", "End", "Ref", "Alt", "VCF_FILTER", 'TUMOR_GT','TUMOR_AD','TUMOR_AF','NORMAL_GT','NORMAL_AD','NORMAL_AF',"Priority.Score","panel_class","Ref.Gene",  
		"Func.refGene.intervar", "Gene.refGeneWithVer", "GeneDetail.refGeneWithVer", "ExonicFunc.refGeneWithVer", "AAChange.refGeneWithVer",
		"clinvar..Clinvar", "InterVar..InterVar.and.Evidence", "PopFreqMax", "gnomAD_exome_ALL", "gnomAD_genome_ALL", 
		"Freq_esp6500siv2_all", "Freq_1000g2015aug_all","dbscSNV_ADA_SCORE.intervar", "dbscSNV_RF_SCORE.intervar", "dpsi_max_tissue", "dpsi_zscore", 
    "spliceai_filtered", "SIFT_score.intervar", "MetaSVM_score.intervar", "CADD_raw.intervar", "CADD_phred.intervar", "GERP.._RS.intervar", "phyloP46way_placental", 
		"Func.refGeneWithVer", "ExonicFunc.refGene.intervar", "avsnp150", "Interpro_domain.intervar"), everything()) %>% 
  replace_na(list(ID = ".")) %>% 
  rename(Chr_annovar = Chr, Start_annovar = Start, Ref_annovar = Ref, Alt_annovar = Alt, End_annovar = End)

write_tsv(annovar_inter, file.path('.', args[3]))
