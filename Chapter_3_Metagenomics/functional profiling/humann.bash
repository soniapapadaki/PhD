###############################################
# 3) FUNCTIONAL ANNOTATION USING HUMANN
#
# Description:
#   Functional annotation using HUMAnN. Output can then be analysed 
#	further for differential gene abundances between samples using 
#	hclust2 and LefSe (see NEOF tutorial).
#
# Notes:
#   This script is based on NEOF's Microbial Shotgun Genomics course:
#   https://neof-workshops.github.io/Shotgun_ld5ug2/Course/01-Shotgun_metagenomics.html
###############################################

#activate biobakery and create new directory
mamba activate biobakery
mkdir 3-FunctionalProfiling
cd 3-FunctionalProfiling

#create folder for concatenated reads 
mkdir concat_reads
cd concat_reads

# Navigate to the 2-Trimmed directory relative to current path
trimmed_dir="../2-Trimmed"
output_dir="."  # Current directory: concat_reads

# Loop through all R1 files in 2-Trimmed
for r1 in "$trimmed_dir"/*R1.fq.gz; do
  # Extract sample name by removing the path and suffix
  filename=$(basename "$r1")
  sample=${filename%R1.fq.gz}

  r2="${trimmed_dir}/${sample}R2.fq.gz"
  output_file="${output_dir}/${sample}.fq.gz"

  # Concatenate if both R1 and R2 exist
  if [[ -f "$r2" ]]; then
    echo "Concatenating $sample..."
    cat "$r1" "$r2" > "$output_file"
  else
    echo "Warning: Missing R2 file for sample $sample, skipping."
  fi
done

#run humann on concatenated reads
cd ..
for file in concat_reads/*.fq.gz; do
    sample=$(basename "$file" .fq.gz)
    humann --input "$file" --output "${sample}_output" --threads 10
done