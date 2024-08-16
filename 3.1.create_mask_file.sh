### prepare mask file
path_to_annotation="/annotation_file"

##pLoF+alphamissense
mask_name="pLOF_alpha_missense"
echo "M_pLOF pLoF" > "${path_to_annotation}"/${mask_name}_masks.txt
echo "M_pLOF_and_alpha_missense pLoF,likely_pathogenic" >> "${path_to_annotation}"/${mask_name}_masks.txt

##pLoF+5 methods missense
mask_name="pLOF_missense"
echo "pLOF_and_55missense pLoF,missense.5in5" > "${path_to_annotation}"/${mask_name}_masks.txt
echo "pLOF_and_55_15missense pLoF,missense.5in5,missense.1in5" >> "${path_to_annotation}"/${mask_name}_masks.txt
