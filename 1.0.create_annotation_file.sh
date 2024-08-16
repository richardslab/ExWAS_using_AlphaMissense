#path to output folder
pathAnnot=Annotation
pathPlugin=Plugins
pathCache=home/.vep

pathSplit=ukb_merged_1-22_indiv_${chr}.vcf

vep -i $pathSplit \
  --plugin dbNSFP,dbNSFP4.3a_${chr}_grch38.gz,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
  --plugin AlphaMissense,"file=/.vep/Plugins/AlphaMissense_hg38.tsv.gz","cols=all" \
  --buffer_size 100000 \
  --force_overwrite \
  --offline \
  --fork 10 \
  --dir_cache $pathCache \
  --everything \
  --cache -o "${pathAnnot}"finalAnnot_${chr}_dbNSFP_and_AlphaMissense.txt;
