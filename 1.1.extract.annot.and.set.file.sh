pathTmp=/tmp/
mkdir -p "${pathTmp}"tmp

#path to the regenie inputs
pathReg=/annotation_file/

for x in {1..22}; do
  #path to the VEP annotation outputs
  pathAnnot=/Annotation/finalAnnot_${x}_dbNSFP_and_AlphaMissense.txt
  #--set-list for regenie below

  awk '(/IMPACT=HIGH/ || /missense_variant/)' $pathAnnot | \
    awk '{ print $4, $1 }' | sort -u -k 1.5,1 -k 2.6,2 > "${pathTmp}"chr${x}.var.gene.txt
  awk '{ print $1 }' "${pathTmp}"chr${x}.var.gene.txt | sort | uniq | \
    awk -v a="$x" '$(NF+1) = a' | \
    awk '$(NF+1) = (FNR FS $(NF+1))' > "${pathTmp}"chr${x}.gene.txt

  gene=($(awk '{print $1}' "${pathTmp}"chr${x}.gene.txt))
  for ((i=0;i<${#gene[@]};++i)); do awk -v a="${gene[i]}" '{if ($1==a) print $2}' "${pathTmp}"chr${x}.var.gene.txt | \
    uniq | \
    awk -f transpose.awk | \
    tr -s '[:blank:]' ","|
    awk 'BEGIN{FS=",";OFS=","} { $1=$1; print $0 }' ; done | \
    paste -d'\0'  "${pathTmp}"chr${x}.gene.txt -  > "${pathTmp}"regenie.set.list.chr${x}.txt

  rm "${pathTmp}"chr${x}.*


#pLOF i.e. VEP high impact
  awk '/IMPACT=HIGH/' $pathAnnot | \
    awk '{ print $1, $4 }' | \
    sort -u -k 1.5,1 -k 2.6,2 | \
    awk '$(NF+1) = "pLoF"' > "${pathTmp}"pLoF.chr${x}.txt

#Missense in general, to be used below
# first only take missense variants with at least one of the five in-silico algorithms annotation
  awk '/missense_variant/ && /IMPACT=MODERATE/ && (/SIFT_pred/ || /Polyphen2_HVAR_pred/ || /Polyphen2_HDIV_pred/ || /MutationTaster_pred/ || /LRT_pred/ || /am_class/)' $pathAnnot > "${pathTmp}"tmp/Missense.annot.chr${x}.txt

#sift
  grep -o -P "(SIFT_pred=[,\.DT]*D[,\.DT]*;)" "${pathTmp}"tmp/Missense.annot.chr${x}.txt | grep -o -P "(=[,\.DT]*D[,\.DT]*;)" > "${pathTmp}"tmp/tmp.sift.count.chr${x}.txt
  grep "SIFT_pred=[,\.DT]*D[,\.DT]*;" "${pathTmp}"tmp/Missense.annot.chr${x}.txt | awk '{ print $1, $4 }' > "${pathTmp}"tmp/tmp.siftVariant.chr${x}.txt
  paste "${pathTmp}"tmp/tmp.siftVariant.chr${x}.txt "${pathTmp}"tmp/tmp.sift.count.chr${x}.txt > "${pathTmp}"tmp/sift.chr${x}.txt
  awk '{ if (substr($3,1) ~ /D/) { print $1, $2, "SIFT_D", $1$2} }' "${pathTmp}"tmp/sift.chr${x}.txt | sort -u -k 1.5,1 -k 2.6,2 > "${pathTmp}"tmp/sift.chr${x}.del.prelim.txt

#LRT
  awk '(/LRT_pred=D;/)' "${pathTmp}"tmp/Missense.annot.chr${x}.txt | awk '{ print $1, $4, "LDR_D", $1$2}' | sort -u -k 1.5,1 -k 2.6,2 > "${pathTmp}"tmp/LRT.chr${x}.del.prelim.txt 

#mutation taster
  awk '(/MutationTaster_pred=[,ADNP\.]*[AD][,ADNP\.]*;/)' "${pathTmp}"tmp/Missense.annot.chr${x}.txt | awk '{ print $1, $4, "mutation_taster_AD", $1$2}' | sort -u -k 1.5,1 -k 2.6,2  > "${pathTmp}"tmp/mutationTaster.chr${x}.del.prelim.txt

#polyphen HVAR  
  grep -o -P "(Polyphen2_HVAR_pred=[,DPB\.]*D[,DPB\.]*;)" "${pathTmp}"tmp/Missense.annot.chr${x}.txt | grep -o -P "(=[,DPB\.]*D[,DPB\.]*;)"  > "${pathTmp}"tmp/tmp.pphenHVAR.count.chr${x}.txt
  grep "Polyphen2_HVAR_pred=[,DPB\.]*D[,DPB\.]*;" "${pathTmp}"tmp/Missense.annot.chr${x}.txt | awk '{ print $1, $4 }' > "${pathTmp}"tmp/tmp.pphenHVARVariant.chr${x}.txt
  paste "${pathTmp}"tmp/tmp.pphenHVARVariant.chr${x}.txt "${pathTmp}"tmp/tmp.pphenHVAR.count.chr${x}.txt > "${pathTmp}"tmp/pphenHVAR.chr${x}.txt
  awk '{ if (substr($3,1) ~ /D/) { print $1, $2, "HVAR_D", $1$2} }' "${pathTmp}"tmp/pphenHVAR.chr${x}.txt | sort -u -k 1.5,1 -k 2.6,2 > "${pathTmp}"tmp/pphenHVAR.chr${x}.del.prelim.txt

#polyphen HDIV
  grep -o -P "(Polyphen2_HDIV_pred=[,DPB\.]*D[,DPB\.]*;)" "${pathTmp}"tmp/Missense.annot.chr${x}.txt | grep -o -P "(=[,DPB\.]*D[,DPB\.]*;)"  > "${pathTmp}"tmp/tmp.pphenHDIV.count.chr${x}.txt
  grep "Polyphen2_HDIV_pred=[,DPB\.]*D[,DPB\.]*;" "${pathTmp}"tmp/Missense.annot.chr${x}.txt | awk '{ print $1, $4 }' > "${pathTmp}"tmp/tmp.pphenHDIVVariant.chr${x}.txt
  paste "${pathTmp}"tmp/tmp.pphenHDIVVariant.chr${x}.txt "${pathTmp}"tmp/tmp.pphenHDIV.count.chr${x}.txt > "${pathTmp}"tmp/pphenHDIV.chr${x}.txt
  awk '{ if (substr($3,1) ~ /D/) { print $1, $2, "HDIV_D", $1$2} }' "${pathTmp}"tmp/pphenHDIV.chr${x}.txt | sort -u -k 1.5,1 -k 2.6,2 > "${pathTmp}"tmp/pphenHDIV.chr${x}.del.prelim.txt

#AlphaMissense
  grep -o -P "am_class=[^;]*" "${pathTmp}"tmp/Missense.annot.chr${x}.txt > "${pathTmp}"tmp/tmp.alphamiss.count.chr${x}.txt
  grep "am_class=[^;]*" "${pathTmp}"tmp/Missense.annot.chr${x}.txt | awk '{ print $1, $4 }' > "${pathTmp}"tmp/tmp.alphamiss.chr${x}.txt
  paste "${pathTmp}"tmp/tmp.alphamiss.chr${x}.txt "${pathTmp}"tmp/tmp.alphamiss.count.chr${x}.txt > "${pathTmp}"tmp/alpha_miss.chr${x}.txt
  awk '{ if (substr($3,1) ~ /likely_pathogenic/) { print $1, $2, "likely_pathogenic", $1$2}}' "${pathTmp}"tmp/alpha_miss.chr${x}.txt | sort -u -k 1.5,1 -k 2.6,2 > "${pathTmp}"tmp/alphamiss.chr${x}.del.prelim.txt
done

###combine set file
cat "${pathTmp}"regenie.set.list.chr* > "${pathReg}"regenie.set.list.txt
