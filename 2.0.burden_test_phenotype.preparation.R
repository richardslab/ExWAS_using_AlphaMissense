library(data.table)
library(dplyr)
library(fastDummies)
library(stringr)

##load EUR list
ukb.eur_list<-fread("ukb.eurIDsPCA.plink.txt")

### load common and rare var PCs
EUR_common_PCs<-fread("eur_commonPCA.txt.eigenvec")
EUR_rare_PCs<-fread("eur_rarePCA.txt.eigenvec")

############################## prepare continuous traits
### load Phenotype: eBMD
BMD_data<-read.table("UKBiobankv2-Qus-Working-Data-Qc.withArrayAndBatch.BMD.txt",header=T)
BMD_data_sub<-BMD_data%>%filter(BATCH !="-NA" & GTARRAY!="-NA")%>%select(FID,IID,AGE, SEX, ZBMD, GTARRAY)

##make the batch variable dummy variable
BMD_data_batch_variable<-BMD_data%>%filter(BATCH !="-NA" & GTARRAY!="-NA")%>%select(FID,IID,BATCH)
BMD_data_batch_variable <- dummy_cols(BMD_data_batch_variable, select_columns = 'BATCH')

##make the center variable dummy variable
BMD_data_center_variable<-BMD_data%>%filter(BATCH !="-NA" & GTARRAY!="-NA")%>%select(FID,IID,CENTRE)
BMD_data_center_variable <- dummy_cols(BMD_data_center_variable, select_columns = 'CENTRE')

## load other phenotypes for GWAS
UKB_continu_traits<-read.csv("UKB_continuous_traits_for_ExWAS_participant.csv")

### merge data
UKB_EUR_continuous_trait_together<-left_join(EUR_rare_PCs, EUR_common_PCs, by=c("FID"="#FID","IID"="IID"))
UKB_EUR_continuous_trait_together1<-left_join(UKB_EUR_continuous_trait_together, BMD_data_sub, by=c("FID"="FID","IID"="IID"))
UKB_EUR_continuous_trait_together2<-left_join(UKB_EUR_continuous_trait_together1,UKB_continu_traits, by=c("FID"="eid"))
UKB_EUR_continuous_trait_together3<-left_join(UKB_EUR_continuous_trait_together2,BMD_data_batch_variable, by=c("FID"="FID","IID"="IID"))
UKB_EUR_continuous_trait_together4<-left_join(UKB_EUR_continuous_trait_together3,BMD_data_center_variable, by=c("FID"="FID","IID"="IID"))

##inverse-rank based normalization
UKB_EUR_continuous_trait_together4$IRNT_TG<-qnorm((rank(UKB_EUR_continuous_trait_together4$p30870_i0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p30870_i0)))
UKB_EUR_continuous_trait_together4$IRNT_SBP<-qnorm((rank(UKB_EUR_continuous_trait_together4$p4080_i0_a0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p4080_i0_a0)))
UKB_EUR_continuous_trait_together4$IRNT_DBP<-qnorm((rank(UKB_EUR_continuous_trait_together4$p4079_i0_a0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p4079_i0_a0)))
UKB_EUR_continuous_trait_together4$IRNT_height<-qnorm((rank(UKB_EUR_continuous_trait_together4$p50_i0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p50_i0)))
UKB_EUR_continuous_trait_together4$IRNT_LDL<-qnorm((rank(UKB_EUR_continuous_trait_together4$p30780_i0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p30780_i0)))
UKB_EUR_continuous_trait_together4$IRNT_biliru<-qnorm((rank(UKB_EUR_continuous_trait_together4$p30660_i0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p30660_i0)))
UKB_EUR_continuous_trait_together4$IRNT_glu<-qnorm((rank(UKB_EUR_continuous_trait_together4$p30740_i0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p30740_i0)))
UKB_EUR_continuous_trait_together4$IRNT_RBC<-qnorm((rank(UKB_EUR_continuous_trait_together4$p30010_i0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p30010_i0)))
UKB_EUR_continuous_trait_together4$IRNT_Ca<-qnorm((rank(UKB_EUR_continuous_trait_together4$p30680_i0,na.last="keep", ties.method="average")-0.5)/sum(!is.na(UKB_EUR_continuous_trait_together4$p30680_i0)))

###select needed variables"
UKB_EUR_continuous_trait_together4_sub<-UKB_EUR_continuous_trait_together4%>%filter(p31 == p22001 & is.na(p21022)==F & is.na(p22001)==F)%>%##keep individuals with the same reported sex and genetic sex; remove individuals without age at recruitment time
  rowwise()%>%mutate(biological_sex=ifelse(p22001 == "Female",0,ifelse(p22001 == "Male",1,NA)),
                     Age_square=p21022*p21022,
                     sex_age_interaction=biological_sex*p21022,
                     sex_age_sq_interaction=biological_sex*Age_square)%>%
  select(FID, IID, starts_with("rarePC"), starts_with("PC"), biological_sex, p21022, Age_square, sex_age_interaction, sex_age_sq_interaction, ZBMD, IRNT_TG, IRNT_SBP, IRNT_DBP, IRNT_height, IRNT_LDL, IRNT_biliru,IRNT_glu,IRNT_RBC, IRNT_Ca)
colnames(UKB_EUR_continuous_trait_together4_sub)<-c("FID", "IID","rarePC1","rarePC2","rarePC3","rarePC4","rarePC5",
                                                    "rarePC6","rarePC7","rarePC8","rarePC9","rarePC10","rarePC11",
                                                    "rarePC12","rarePC13","rarePC14","rarePC15","rarePC16","rarePC17",
                                                    "rarePC18","rarePC19","rarePC20","PC1","PC2","PC3",
                                                    "PC4","PC5","PC6","PC7","PC8","PC9",
                                                    "PC10","biological_sex","Age","Age_square","sex_age_interaction","sex_age_sq_interaction",
                                                    "ZBMD","IRNT_TG","IRNT_SBP","IRNT_DBP","IRNT_height","IRNT_LDL",
                                                    "IRNT_biliru","IRNT_glu","IRNT_RBC","IRNT_Ca")

###remove withdraw
with_draw_list<-read.csv("w27449_2023-04-25.csv",header=F)
UKB_EUR_continuous_trait_together4_sub_no_withdraw<-UKB_EUR_continuous_trait_together4_sub%>%filter(!(FID %in% with_draw_list$V1))
write.table(UKB_EUR_continuous_trait_together4_sub_no_withdraw, "UKB_eur_continuous_phenotypes.txt",
            sep="\t",row.names=F, quote=F)

######################### process binary traits #########################
mendelvar_disease_genes<-fread("mendelvar_disease.txt")

##check Phecode and ICD10 code matching table
Phecode_map_icd10<-read.csv("Phecode_map_v1_2_icd10cm_beta.csv")

## exclude some phecodes
Phecode_map_icd10_filtered<-Phecode_map_icd10%>%filter(!(exclude_name %in% c("injuries & poisonings","","NULL")))

#check the numebr of phecode categories
phecode_list<-unique(Phecode_map_icd10_filtered$phecode)

Phecode_map_icd10_filtered_annot<-Phecode_map_icd10_filtered%>%select(phecode, phecode_str)%>%rowwise()%>%
mutate(phecode_modi=paste0("X",phecode),
       phecode_modi_cat=str_split(phecode_modi, "\\.")[[1]][1])%>%select(phecode_modi, phecode_str, phecode_modi_cat)%>%
  distinct()

Phecode_map_icd10_filtered_annot_icd10<-Phecode_map_icd10_filtered%>%
select(phecode, phecode_str, icd10cm)%>%rowwise()%>%mutate(phecode_modi=paste0("X",phecode),
                                                           phecode_modi_cat=str_split(phecode_modi, "\\.")[[1]][1])%>%
select(phecode_modi, phecode_str, phecode_modi_cat,icd10cm)%>%
distinct()

##### load phecode counts data for all UKB EUR individuals (phecode set 2)
UKB_EUR_phecode_case_ALL1<-read.csv("UKB_all_individuals_with_phecode_diseases_set1.csv", row.names=1)
UKB_EUR_phecode_case_ALL2<-read.csv("UKB_all_individuals_with_phecode_diseases_set2.csv", row.names=1)
UKB_EUR_phecode_case_ALL3<-read.csv("UKB_all_individuals_with_phecode_diseases_set3.csv", row.names=1)
UKB_EUR_phecode_case_ALL4<-read.csv("UKB_all_individuals_with_phecode_diseases_set4.csv", row.names=1)

UKB_EUR_phecode_case_ALL<-left_join(UKB_EUR_phecode_case_ALL1, UKB_EUR_phecode_case_ALL2, by=c("eid"="eid"))
UKB_EUR_phecode_case_ALL<-left_join(UKB_EUR_phecode_case_ALL, UKB_EUR_phecode_case_ALL3, by=c("eid"="eid"))
UKB_EUR_phecode_case_ALL<-left_join(UKB_EUR_phecode_case_ALL, UKB_EUR_phecode_case_ALL4, by=c("eid"="eid"))
## get the sum of cases for each disease
sums_cases_df <- data.frame(
  Column = names(UKB_EUR_phecode_case_ALL[,-1]),
  Sum = colSums(UKB_EUR_phecode_case_ALL[,-1])
)
sums_cases_df1<-left_join(sums_cases_df, Phecode_map_icd10_filtered_annot, by=c("Column"="phecode_modi"))

########################### obtain the traits with enough cases and mendel genes ###########################
##keep those phecodes with at least 10000 cases
sums_cases_df1_case_filter<-sums_cases_df1%>%filter(Sum >=10000)

DOID_for_top_phecodes<-read.csv("DOID_for_top_phecodes.csv")
DOID_for_top_phecodes1<-DOID_for_top_phecodes%>%rowwise()%>%mutate(DO_index=str_split(DOIID,":")[[1]][2])

###keep those with MendelVar Genes
mendelvar_genes_ALL<-NULL
mendelvar_genes_ALL_counts<-NULL
for (i in 1:length(unique(DOID_for_top_phecodes1$phecode))){
  ##traits
  phecode_test=unique(DOID_for_top_phecodes1$phecode)[i]
  ##DOI list
  DOID_for_top_phecodes1_sub<-DOID_for_top_phecodes1%>%filter(phecode ==phecode_test)
  ##get the number of genes for each traits
  mendelvar_genes_temp<-mendelvar_disease_genes%>%filter(do %in% DOID_for_top_phecodes1_sub$DO_index)%>%select(hgnc_gene_name,ensg_gene_name)%>%distinct()
  mendelvar_genes_temp$trait=phecode_test
  n_genes<-length(unique(mendelvar_genes_temp$hgnc_gene_name))
  mendelvar_genes_ALL<-rbind(mendelvar_genes_ALL,mendelvar_genes_temp)
  mendelvar_genes_ALL_counts<-rbind(mendelvar_genes_ALL_counts,data.frame(phecode_col=phecode_test,n_genes=n_genes))
}

##select traits that with at least 1 mendel genes
mendelvar_genes_ALL_counts_with_at_least_1_mendelgene<-mendelvar_genes_ALL_counts%>%filter(n_genes>0)%>%rowwise()%>%
  mutate(phecode_str=Phecode_map_icd10_filtered_annot[which(Phecode_map_icd10_filtered_annot$phecode_modi == phecode_col),]$phecode_str,
         phencode_cat=Phecode_map_icd10_filtered_annot[which(Phecode_map_icd10_filtered_annot$phecode_modi == phecode_col),]$phecode_modi_cat)


########################### get the phecode category for phecodes with enough cases and mendel genes ###########################
## get the phecode with the same phecode categories
sums_cases_df1_sub<-sums_cases_df1%>%filter(phecode_modi_cat %in% mendelvar_genes_ALL_counts_with_at_least_1_mendelgene$phencode_cat)
sums_cases_df1_sub_triats_with_over10000_cases<-sums_cases_df1_sub%>%filter(Sum>10000)

############# Prepare binary phenotype files #############
UKB_EUR_together4_sub_no_withdraw_covar<-UKB_EUR_continuous_trait_together4_sub_no_withdraw%>%select(FID:sex_age_sq_interaction)
UKB_EUR_phecode_case_ALL_sub<-UKB_EUR_phecode_case_ALL%>%select(eid, any_of(mendelvar_genes_ALL_counts_with_at_least_1_mendelgene$phecode_col))
UKB_EUR_binary_trait_together4_sub_no_withdraw<-left_join(UKB_EUR_together4_sub_no_withdraw_covar, UKB_EUR_phecode_case_ALL_sub, by=c("FID"="eid"))
colnames(UKB_EUR_binary_trait_together4_sub_no_withdraw)<-c("FID","IID","rarePC1","rarePC2","rarePC3","rarePC4","rarePC5",
                                                            "rarePC6","rarePC7","rarePC8","rarePC9","rarePC10","rarePC11",
                                                            "rarePC12","rarePC13","rarePC14","rarePC15","rarePC16","rarePC17",
                                                            "rarePC18","rarePC19","rarePC20","PC1","PC2","PC3",
                                                            "PC4","PC5","PC6","PC7","PC8","PC9",
                                                            "PC10","biological_sex","Age","Age_square","sex_age_interaction","sex_age_sq_interaction",
                                                            "hypertension","Hypercholesterolemia","Diaphragmatic_hernia","Osteoarthritis_localized","Cataract","T2D",
                                                            "Obesity","Osteoarthrosis_NOS","Major_depressive_disorder","Malignant_neoplasm","Hypothyroidism","Acute_renal_failure",
                                                            "Atrial_fibrillation","Cancer_of_prostate","Breast_cancer")

write.table(UKB_EUR_binary_trait_together4_sub_no_withdraw, "UKB_eur_binary_phenotypes.txt",
            sep="\t",row.names=F, quote=F)

## create a list of outcomes for burden test
continuous_trait_list<-data.frame(trait=c("ZBMD","IRNT_TG","IRNT_SBP","IRNT_DBP","IRNT_height","IRNT_LDL",
                                  "IRNT_biliru","IRNT_glu","IRNT_RBC","IRNT_Ca"))
write.table(continuous_trait_list, "UKB_continuous_trait_list.txt",
            sep="\t",row.names=F, quote=F)

binary_trait_list<-data.frame(trait=c("hypertension","Hypercholesterolemia","Diaphragmatic_hernia","Osteoarthritis_localized","Cataract","T2D",
                                      "Major_depressive_disorder","Hypothyroidism","Acute_renal_failure",
                                      "Atrial_fibrillation","Cancer_of_prostate","Breast_cancer"))

write.table(binary_trait_list, "UKB_binary_trait_list.txt",
            sep="\t",row.names=F, quote=F)
