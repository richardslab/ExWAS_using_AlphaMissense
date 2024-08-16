genotype_pruned_step1="eur_autosomes_ld_pruned"
WES_genotype_prefix="ukb_merged_1-22"
phenotype_file="UKB_eur_binary_phenotypes.txt"
EUR_sample_list="ukb.eurIDsPCA.plink.txt"
output_dir="results"
path_to_annotation="annotation_file"

mask_name="pLOF_missense"
annotation_file_name="regenie.anno.file.pLOF.missense.txt"

#### step 2
regenie \
    --step 2 \
    --pred ${output_dir}/fit_step1_TRAIT_pred.list \
    --pgen ${WES_genotype_prefix} \
    --ref-first \
    --keep ${EUR_sample_list} \
    --phenoFile $phenotype_file \
    --covarFile $phenotype_file \
    --bt \
    --firth --approx \
    --phenoCol TRAIT \
    --covarColList Age,Age_square,sex_age_interaction,sex_age_sq_interaction,PC{1:10},rarePC{1:20} \
    --catCovarList biological_sex \
    --set-list "${path_to_annotation}/regenie.set.list.txt" \
    --anno-file "${path_to_annotation}/${annotation_file_name}" \
    --mask-def "${path_to_annotation}/${mask_name}_masks.txt" \
    --vc-tests skato \
    --vc-maxAAF 0.01 \
    --build-mask max \
    --aaf-bins 0.01,0.001\
    --joint acat \
    --bsize 200 \
    --out ${output_dir}/step_2_TRAIT_ExWAS_${mask_name}
