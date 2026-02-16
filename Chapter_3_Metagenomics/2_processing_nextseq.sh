#!/bin/bash



# Trimming the 151st base (NextSeq)

for forward in *_1.fastq; do reverse=${forward/_1/_2}; base="$(basename -- $forward | sed 's/_.*//')"; bbduk.sh -Xmx10g in1=$forward in2=$reverse out1="$base"_lbt_1.fastq out2="$base"_lbt_2.fastq ftr=149 t=12; done


# Trimming adapters

for forward in *lbt_1.fastq; do reverse=${forward/_1/_2}; base="$(basename -- $forward | sed 's/_.*//')"; bbduk.sh -Xmx10g in1=$forward in2=$reverse out1="$base"_atrim_1.fastq out2="$base"_atrim_2.fastq ref=~/glacier-mtg/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats="$base"_atrim_stats.txt tpe tbo t=12; done


# Removing PhiX

for forward in *_atrim_1.fastq; do reverse=${forward/_1/_2}; base="$(basename -- $forward | sed 's/_.*//')"; bbduk.sh -Xmx10g in1=$forward in2=$reverse out1="$base"_phix_1.fastq out2="$base"_phix_2.fastq ref=~/glacier-mtg/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats="$base"_phix_stats.txt t=12; done


# Quality trimming (Reads where the quality scores fall towards the end are shortened)

for forward in *_phix_1.fastq; do reverse=${forward/_1/_2}; base="$(basename -- $forward | sed 's/_.*//')"; bbduk.sh -Xmx10g in1=$forward in2=$reverse out1="$base"_processed_1.fastq out2="$base"_processed_2.fastq qtrim=r trimq=15 minlength=50 t=12 stats="$base"_qtrim_stat; done


# Delete intermediate files

rm *lbt*.fastq; rm *atrim*.fastq; rm *phix*.fastq


#------------


# Create fastqc

mkdir ./fastqc_processed; fastqc *processed*.fastq -o ./fastqc_processed -t 5


# And summarise with MultiQC

cd ./fastqc_processed; conda activate multiqc; multiqc *_fastqc.zip


#------------

cd ..

# Add a list of sample IDs

for file in *processed_1.fastq; do sample="$(basename -- $file | sed 's/_.*//')"; echo $sample >> temp.txt; done


# Number of raw pairs per sample (this is total reads, will be divided by 2 later)

for file in *atrim_stats.txt; do grep 'Total' $file | cut -f 2 >> temp.txt; done


# Number of read pairs remaining after preprocessing

seqkit stats *processed_1.fastq -T -j 12 | cut -f 4 | tail -n +2 >> temp.txt


# Shift this long data into columns
# pr will split into 'pages' of 66 lines each, -l 1000 changes that to 1000 # (i.e. unless >1000 samples, should be the 4 columns you want)

pr -ts --columns 3 -l 1000 temp.txt > temp2.txt


# Divide those raw read pairs by 2

awk '{print $1 "\t" $2/2 "\t" $3}' temp2.txt > temp3.txt # Add a column for the percentage of reads retained
awk '{print $1 "\t" $2 "\t" $3 "\t" ($3/$2)*100}' temp3.txt > temp4.txt # Add headers to the final file
echo -e "Sample\tRaw.pairs\tPairs.remain\tPercent.remain" | cat - temp4.txt > read_summary; rm temp*.txt

