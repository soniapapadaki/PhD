# NEOF’s Microbial Shotgun Genomics 2025 course



#first set up directory "1-Raw" with symbolic links to raw data and change to that directory
#activate shotgun_meta environment
mamba activate shotgun_meta

#R1 fastqc
#Make an output directory
mkdir R1_fastqc
#Run fastqc on all the R1.fastq.gz files 
#* matches any pattern
#*R1.fastq.gz matches any file that ends R1.fastq.gz in the current directory
#-t 3 indicates to use 9 threads
fastqc -t 9 -o R1_fastqc *R1.fastq.gz

#R2 fastqc
#Make output directory
mkdir R2_fastqc
#Run fastqc
fastqc -t 9 -o R2_fastqc *R2.fastq.gz

#R1 multiqc fastqc report
#Create output directory
mkdir R1_fastqc/multiqc
#Create multiqc output
multiqc -o R1_fastqc/multiqc R1_fastqc

#R2 multiqc fastqc report
#Create output directory
mkdir R2_fastqc/multiqc
#Create multiqc report
multiqc -o R2_fastqc/multiqc R2_fastqc

#go back to project directory and create directory for trimmed sequences
cd ..
mkdir 2-Trimmed
cd 2-Trimmed

#trim reads with trim galore
for r1 in ../1-Raw/*R1.fastq.gz; do
    r2=${r1/R1.fastq.gz/R2.fastq.gz}
    sample=$(basename "$r1" R1.fastq.gz)

    trim_galore --paired --quality 20 --stringency 4 "$r1" "$r2"
done

#rename files
for file in *_val_?.fq.gz; do
  newname="${file/_val_[12]/}"
  mv "$file" "$newname"
done

#create fastqc and multiqc reports with trimmed reads
#R1
mkdir R1_fastqc
fastqc -t 9 -o R1_fastqc *R1.fq.gz
mkdir R1_fastqc/multiqc
multiqc -o R1_fastqc/multiqc R1_fastqc
#R2
mkdir R2_fastqc
fastqc -t 9 -o R2_fastqc *R2.fq.gz
mkdir R2_fastqc/multiqc
multiqc -o R2_fastqc/multiqc R2_fastqc