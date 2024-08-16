genotype_pruned_step1="eur_autosomes_ld_pruned"
WES_genotype_prefix="ukb_merged_1-22"
phenotype_file="UKB_eur_binary_phenotypes.txt"
EUR_sample_list="ukb.eurIDsPCA.plink.txt"
output_dir="results"
path_to_annotation="Annotation"


### binary trait example
regenie \
  --step 1 \
  --bed $genotype_pruned_step1 \
  --covarFile $phenotype_file \
  --phenoFile $phenotype_file \
  --keep ${EUR_sample_list} \
  --phenoCol TRAIT \
  --covarColList Age,Age_square,sex_age_interaction,sex_age_sq_interaction,PC{1:10},rarePC{1:20} \
  --catCovarList biological_sex \
  --bt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmp_rg_TRAIT \
  --out fit_step1_TRAIT
