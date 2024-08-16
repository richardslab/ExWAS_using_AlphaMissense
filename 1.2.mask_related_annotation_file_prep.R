library(data.table)
library(dplyr)

############# aggregate missense annotation 5 methods #############
ALL_5_methods_annot_all<-NULL

for (i in 1:22){
  ##SIFT
  SIFT_annot_tmp<-fread(paste0("/Annotation/tmp/tmp/sift.chr",i,".del.prelim.txt"),header=F)
  ##LRT
  LRT_annot_tmp<-fread(paste0("/Annotation/tmp/tmp/LRT.chr",i,".del.prelim.txt"),header=F)
  ##MutationTaster
  MutationTaster_annot_tmp<-fread(paste0("/Annotation/tmp/tmp/mutationTaster.chr",i,".del.prelim.txt"),header=F)
  ##pphenHVAR
  pphenHVAR_annot_tmp<-fread(paste0("/Annotation/tmp/tmp/pphenHVAR.chr",i,".del.prelim.txt"),header=F)
  ##pphenHDIV
  pphenHDIV_annot_tmp<-fread(paste0("/Annotation/tmp/tmp/pphenHDIV.chr",i,".del.prelim.txt"),header=F)
  ##get a list of all unique pairs of variant-gene pairs
  SIFT_annot_tmp_sub<-SIFT_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(SIFT_annot_tmp_sub)<-c("V1","V2","V3_SIFT")
  LRT_annot_tmp_sub<-LRT_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(LRT_annot_tmp_sub)<-c("V1","V2","V3_LRT")
  MutationTaster_annot_tmp_sub<-MutationTaster_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(MutationTaster_annot_tmp_sub)<-c("V1","V2","V3_MT")
  pphenHVAR_annot_tmp_sub<-pphenHVAR_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(pphenHVAR_annot_tmp_sub)<-c("V1","V2","V3_HVAR")
  pphenHDIV_annot_tmp_sub<-pphenHDIV_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(pphenHDIV_annot_tmp_sub)<-c("V1","V2","V3_HDIV")
  ALL_5_methods_annot<-rbind(SIFT_annot_tmp_sub[,1:2], LRT_annot_tmp_sub[,1:2], MutationTaster_annot_tmp_sub[,1:2], pphenHVAR_annot_tmp_sub[,1:2],pphenHDIV_annot_tmp_sub[,1:2])%>%distinct()
  ##create missense.5in5 and missense.1in5
  deli_5_methods_check <- ALL_5_methods_annot %>%left_join(SIFT_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>% left_join(LRT_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>% left_join(MutationTaster_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>%
  left_join(pphenHVAR_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>%left_join(pphenHDIV_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>% rowwise()%>%
  mutate(deli_5_methods_check = case_when(is.na(V3_SIFT) | is.na(V3_LRT) | is.na(V3_MT) | is.na(V3_HVAR) | is.na(V3_HDIV) ~ "missense.1in5",
    TRUE ~ "missense.5in5"
  )) %>%select(-starts_with("V3_"))
  write.table(deli_5_methods_check, paste0("/Annotation/tmp/missense.5methods.chr",i,".txt"), sep="\t", quote=F, row.names=F, col.names=F)
  #deli_5_methods_check<-read.table(paste0("/Annotation/tmp/missense.5methods.chr",i,".txt"), sep="\t")
  ALL_5_methods_annot_all<-rbind(ALL_5_methods_annot_all, deli_5_methods_check)
}
write.table(ALL_5_methods_annot_all, "/annotation_file/regenie.anno.file.5methods_missense.txt", quote=F,row.names = F, col.names =F, sep="\t")                                                                  


### only likely pathogeneic (alphamissense)
annt_file_alpha_miss<-NULL
for (i in 1:22){
  ##Alphamissense
  alphamiss_annot_tmp<-fread(paste0("/Annotation/tmp/tmp/alphamiss.chr",i,".del.prelim.txt"),header=F)
  alphamiss_annot_tmp_sub<-alphamiss_annot_tmp%>%select(V1,V2,V3)
  annt_file_alpha_miss<-rbind(annt_file_alpha_miss,alphamiss_annot_tmp_sub)
}
write.table(annt_file_alpha_miss, "/annotation_file/regenie.anno.file.alphamiss.likely.patho.txt", quote=F,row.names = F, col.names =F, sep="\t")

### likely pathogeneic (alphamissense) and pLoF
annt_file_pLOF<-NULL
for (i in 1:22){
  ##Alphamissense
  pLOF_annot_tmp<-fread(paste0("/Annotation/tmp/pLoF.chr",i,".txt"),header=F)
  annt_file_pLOF<-rbind(annt_file_pLOF,pLOF_annot_tmp)
}
write.table(annt_file_pLOF, "/annotation_file/regenie.anno.file.pLOF.txt", quote=F,row.names = F, col.names =F, sep="\t")


######################## create annotation file used for burden test ########################
annt_file_pLOF<-fread("/annotation_file/regenie.anno.file.pLOF.txt", header=F)
nrow(annt_file_pLOF)
annt_file_alpha_miss<-fread("/annotation_file/regenie.anno.file.alphamiss.likely.patho.txt", header=F)

annt_file<-rbind(annt_file_pLOF, annt_file_alpha_miss,ALL_5_methods_annot_all)
annt_file_alpha_miss_pLOF<-annt_file%>%filter(V3 =="likely_pathogenic" | V3 =="pLoF")%>%group_by(V1, V2) %>%
  arrange(match(V3, c("pLoF", "likely_pathogenic"))) %>%
  slice(1)%>%select(V1,V2,V3)%>%distinct()

write.table(annt_file_alpha_miss_pLOF, "/annotation_file/regenie.anno.file.pLOF.alphamiss.likely.patho.txt", quote=F,row.names = F, col.names =F, sep="\t")

### pLoF,missense.1in5, missense.5in5
annt_file_pLOF_missense<-annt_file%>%filter(V3 !="likely_pathogenic")%>%group_by(V1, V2) %>%
  arrange(match(V3, c("pLoF","missense.5in5","missense.1in5"))) %>%
  slice(1)%>%select(V1,V2,V3)%>%distinct()
write.table(annt_file_pLOF_missense, "/annotation_file/regenie.anno.file.pLOF.missense.txt", quote=F,row.names = F, col.names =F, sep="\t")

### pLoF, missense.1in5, missense.5in5, alphamissense
annt_file_pLOF_missense_and_alpha_missense<-annt_file%>%group_by(V1, V2) %>%
  arrange(match(V3, c("pLoF","likely_pathogenic","missense.5in5","missense.1in5"))) %>%
  slice(1)%>%select(V1,V2,V3)%>%distinct()
write.table(annt_file_pLOF_missense_and_alpha_missense, "/annotation_file/regenie.anno.file.pLOF.missense.alphamissense.txt", quote=F,row.names = F, col.names =F, sep="\t")


################################################################################
############# aggregate missense annotation 5 methods (canonical) #############
ALL_5_methods_annot_all<-NULL

for (i in 1:22){
  ##SIFT
  SIFT_annot_tmp<-fread(paste0("/Annotation/tmp_CANONICAL/tmp/sift.chr",i,".del.prelim.txt"),header=F)
  ##LRT
  LRT_annot_tmp<-fread(paste0("/Annotation/tmp_CANONICAL/tmp/LRT.chr",i,".del.prelim.txt"),header=F)
  ##MutationTaster
  MutationTaster_annot_tmp<-fread(paste0("/Annotation/tmp_CANONICAL/tmp/mutationTaster.chr",i,".del.prelim.txt"),header=F)
  ##pphenHVAR
  pphenHVAR_annot_tmp<-fread(paste0("/Annotation/tmp_CANONICAL/tmp/pphenHVAR.chr",i,".del.prelim.txt"),header=F)
  ##pphenHDIV
  pphenHDIV_annot_tmp<-fread(paste0("/Annotation/tmp_CANONICAL/tmp/pphenHDIV.chr",i,".del.prelim.txt"),header=F)
  ##get a list of all unique pairs of variant-gene pairs
  SIFT_annot_tmp_sub<-SIFT_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(SIFT_annot_tmp_sub)<-c("V1","V2","V3_SIFT")
  LRT_annot_tmp_sub<-LRT_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(LRT_annot_tmp_sub)<-c("V1","V2","V3_LRT")
  MutationTaster_annot_tmp_sub<-MutationTaster_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(MutationTaster_annot_tmp_sub)<-c("V1","V2","V3_MT")
  pphenHVAR_annot_tmp_sub<-pphenHVAR_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(pphenHVAR_annot_tmp_sub)<-c("V1","V2","V3_HVAR")
  pphenHDIV_annot_tmp_sub<-pphenHDIV_annot_tmp%>%dplyr::select(V1,V2,V3)
  colnames(pphenHDIV_annot_tmp_sub)<-c("V1","V2","V3_HDIV")
  ALL_5_methods_annot<-rbind(SIFT_annot_tmp_sub[,1:2], LRT_annot_tmp_sub[,1:2], MutationTaster_annot_tmp_sub[,1:2], pphenHVAR_annot_tmp_sub[,1:2],pphenHDIV_annot_tmp_sub[,1:2])%>%distinct()
  ##create missense.5in5 and missense.1in5
  deli_5_methods_check <- ALL_5_methods_annot %>%left_join(SIFT_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>% left_join(LRT_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>% left_join(MutationTaster_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>%
  left_join(pphenHVAR_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>%left_join(pphenHDIV_annot_tmp_sub, by = c("V1"="V1","V2"="V2")) %>% rowwise()%>%
  mutate(deli_5_methods_check = case_when(is.na(V3_SIFT) | is.na(V3_LRT) | is.na(V3_MT) | is.na(V3_HVAR) | is.na(V3_HDIV) ~ "missense.1in5",
    TRUE ~ "missense.5in5"
  )) %>%select(-starts_with("V3_"))
  write.table(deli_5_methods_check, paste0("/Annotation/tmp_CANONICAL/missense.5methods.chr",i,".txt"), sep="\t", quote=F, row.names=F, col.names=F)
  ALL_5_methods_annot_all<-rbind(ALL_5_methods_annot_all, deli_5_methods_check)
}
write.table(ALL_5_methods_annot_all, "/annotation_file/regenie.anno.file.5methods_missense_CANONICAL.txt", quote=F,row.names = F, col.names =F, sep="\t")                                                                  


### only likely pathogenic (alphamissense)
annt_file_alpha_miss<-NULL
for (i in 1:22){
  ##Alphamissense
  alphamiss_annot_tmp<-fread(paste0("/Annotation/tmp_CANONICAL/tmp/alphamiss.chr",i,".del.prelim.txt"),header=F)
  alphamiss_annot_tmp_sub<-alphamiss_annot_tmp%>%select(V1,V2,V3)
  annt_file_alpha_miss<-rbind(annt_file_alpha_miss,alphamiss_annot_tmp_sub)
}
write.table(annt_file_alpha_miss, "/annotation_file/regenie.anno.file.alphamiss.likely.patho_CANONICAL.txt", quote=F,row.names = F, col.names =F, sep="\t")

### likely pathogeneic (alphamissense) and pLoF
annt_file_pLOF<-NULL
for (i in 1:22){
  ##Alphamissense
  pLOF_annot_tmp<-fread(paste0("/Annotation/tmp_CANONICAL/pLoF.chr",i,".txt"),header=F)
  annt_file_pLOF<-rbind(annt_file_pLOF,pLOF_annot_tmp)
}
write.table(annt_file_pLOF, "/annotation_file/regenie.anno.file.pLOF_CANONICAL.txt", quote=F,row.names = F, col.names =F, sep="\t")


######################## create annotation file used for burden test ########################
annt_file_pLOF<-fread("/annotation_file/regenie.anno.file.pLOF_CANONICAL.txt", header=F)
annt_file_alpha_miss<-fread("/annotation_file/regenie.anno.file.alphamiss.likely.patho_CANONICAL.txt", header=F)
ALL_5_methods_annot_all<-fread("/annotation_file/regenie.anno.file.5methods_missense_CANONICAL.txt", header=F)

annt_file<-rbind(annt_file_pLOF, annt_file_alpha_miss,ALL_5_methods_annot_all)

annt_file_alpha_miss_pLOF<-annt_file%>%filter(V3 =="likely_pathogenic" | V3 =="pLoF")%>%group_by(V1, V2) %>%
  arrange(match(V3, c("pLoF", "likely_pathogenic"))) %>%
  slice(1)%>%select(V1,V2,V3)%>%distinct()

nrow(annt_file_alpha_miss_pLOF) == nrow(annt_file_alpha_miss_pLOF%>%select(V1,V2)%>%distinct())

write.table(annt_file_alpha_miss_pLOF, "/annotation_file/regenie.anno.file.pLOF.alphamiss.likely.patho_CANONICAL.txt", quote=F,row.names = F, col.names =F, sep="\t")

### pLoF,missense.1in5, missense.5in5
annt_file_pLOF_missense<-annt_file%>%filter(V3 !="likely_pathogenic")%>%group_by(V1, V2) %>%
  arrange(match(V3, c("pLoF","missense.5in5","missense.1in5"))) %>%
  slice(1)%>%select(V1,V2,V3)%>%distinct()
write.table(annt_file_pLOF_missense, "/annotation_file/regenie.anno.file.pLOF.missense_CANONICAL.txt", quote=F,row.names = F, col.names =F, sep="\t")

### pLoF, missense.1in5, missense.5in5, alphamissense
annt_file_pLOF_missense_and_alpha_missense<-annt_file%>%group_by(V1, V2) %>%
  arrange(match(V3, c("pLoF","likely_pathogenic","missense.5in5","missense.1in5"))) %>%
  slice(1)%>%select(V1,V2,V3)%>%distinct()

write.table(annt_file_pLOF_missense_and_alpha_missense, "/annotation_file/regenie.anno.file.pLOF.missense.alphamissense_CANONICAL.txt", quote=F,row.names = F, col.names =F, sep="\t")
