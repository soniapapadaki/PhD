###############################################
# QIIME2 (18S amplicon sequencing)
#
# Description:
#   Workflow for sequence quality control, ASV identification, and taxonomic 
#	classification with 18S rRNA amplicon sequencing data.
#
# Notes:
#   This script is adapted from the QIIME Moving Pictures Tutorial:
#   https://docs.qiime2.org/2024.10/tutorials/moving-pictures/
###############################################

#!/bin/bash

## IMPORTING DATA

#Create new directory
mkdir endoliths_sc_18S
cd endoliths_sc_18S

#Import sequencing data, metadata and manifest file
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest-file-18S-endoliths-sc.tsv \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

#Create summary file of sequence quality
qiime demux summarize \
--i-data paired-end-demux.qza \
--o-visualization paired-end-demux.qzv

## SEQUENCE QUALITY CONTROL AND FEATURE TABLE CONSTRUCTION

##DADA2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trunc-len-f 130 \
--p-trunc-len-r 130 \
--o-representative-sequences rep-seqs.qza \
--o-table table.qza \
--o-denoising-stats stats-dada2.qza \
--verbose

qiime metadata tabulate \
--m-input-file stats-dada2.qza \
--o-visualization stats-dada2.qzv

#Create feature table summaries
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file endolith_metadata.tsv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

## TAXONOMIC CLASSIFICATION

#extract reference reads from SILVA 138 database
qiime feature-classifier extract-reads \
--i-sequences silva-138-99-seqs.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--o-reads ref-seqs.qza

#train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier 18S-classifier.qza

#use trained classifier to create taxonomy file for ASVs
qiime feature-classifier classify-sklearn \
--i-classifier 18S-classifier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

#Create taxa bar plots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file endolith_metadata.tsv \
--o-visualization taxa-bar-plots.qzv

## FILTER OUT BACTERIA, ARCHAEA, METAZOA, AND UNASSIGNED FOR 18S

qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Unassigned,d__Bacteria,d__Archaea,p__Rotifera,p__Nematozoa,p__Arthropoda,p__Vertebrata,p__Tardigrada,p__Annelida,p__Gastrotricha,p__Platyhelminthes,c__Embryophyta \
  --p-include p__ \
  --o-filtered-table table-filtered.qza
 
qiime taxa barplot \
  --i-table table-filtered.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered.qzv

qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Unassigned,d__Bacteria,d__Archaea,p__Rotifera,p__Nematozoa,p__Arthropoda,p__Vertebrata,p__Tardigrada,p__Annelida,p__Gastrotricha,p__Platyhelminthes,c__Embryophyta \
  --p-include p__ \
  --o-filtered-sequences rep-seqs-filtered.qza

qiime feature-table summarize \
  --i-table table-filtered.qza \
  --o-visualization table-filtered.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-filtered.qza \
  --o-visualization rep-seqs-filtered.qzv
  
## PRODUCE RELATIVE ABUNDANCE PLOTS FOR TOP 4 PHYLA

#ascomycota
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Bacteria,d__Archaea \
  --p-include p__Ascomycota \
  --o-filtered-table table-filtered-only-asc.qza
  
qiime taxa barplot \
  --i-table table-filtered-only-asc.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-only-asc.qzv
  
qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Bacteria,d__Archaea \
  --p-include p__Ascomycota \
  --o-filtered-sequences rep-seqs-filtered-only-asc.qza
  
#chlorophyta
 qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Bacteria,d__Archaea \
  --p-include p__Chlorophyta \
  --o-filtered-table table-filtered-only-cpt.qza
  
qiime taxa barplot \
  --i-table table-filtered-only-cpt.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-only-cpt.qzv
  
qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Bacteria,d__Archaea \
  --p-include p__Chlorophyta \
  --o-filtered-sequences rep-seqs-filtered-only-cpt.qza
 
#ciliophora
 qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Bacteria,d__Archaea \
  --p-include p__Ciliophora \
  --o-filtered-table table-filtered-only-cil.qza
  
qiime taxa barplot \
  --i-table table-filtered-only-cil.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-only-cil.qzv
  
qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Bacteria,d__Archaea \
  --p-include p__Ciliophora \
  --o-filtered-sequences rep-seqs-filtered-only-cil.qza

#basidiomycota
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Bacteria,d__Archaea,c__Malasseziomycetes \
  --p-include p__Basidiomycota \
  --o-filtered-table table-filtered-only-bas.qza
  
qiime taxa barplot \
  --i-table table-filtered-only-bas.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-only-bas.qzv
  
qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Bacteria,d__Archaea,c__Malasseziomycetes \
  --p-include p__Basidiomycota \
  --o-filtered-sequences rep-seqs-filtered-only-bas.qza

## GENERATE PHYLOGENETIC TREE
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs-filtered.qza \
--o-alignment aligned-rep-seqs-filtered.qza \
--o-masked-alignment masked-aligned-rep-seqs-filtered.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

## ALPHA RAREFACTION - choosing max depth based on median in table-filtered.qzv

qiime diversity alpha-rarefaction \
--i-table table-filtered.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 56943 \
--p-steps 20 \
--m-metadata-file endolith_metadata.tsv \
--o-visualization alpha-rarefaction.qzv

qiime diversity alpha-rarefaction \
--i-table table-filtered.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 56943 \
--p-steps 20 \
--o-visualization alpha-rarefaction-eachsample.qzv

## ALPHA AND BETA DIVERSITY ANALYSIS - using filtered data

#Generate core diversity metrics
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table-filtered.qza \
--p-sampling-depth 4935 \
--m-metadata-file endolith_metadata.tsv \
--output-dir core-metrics-results

#Test for associations between metadata categories and alpha diversity data
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
--m-metadata-file endolith_metadata.tsv \
--o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/evenness_vector.qza \
--m-metadata-file endolith_metadata.tsv \
--o-visualization core-metrics-results/evenness-group-significance.qzv

#For beta diversity, compare between sample types and sites
qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file endolith_metadata.tsv \
--m-metadata-column sample_type \
--o-visualization core-metrics-results/unweighted-unifrac-sample-type-significance.qzv \
--p-pairwise

qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file endolith_metadata.tsv \
--m-metadata-column sample_site \
--o-visualization core-metrics-results/unweighted-unifrac-sample-site-group-significance.qzv \
--p-pairwise

## DIFFERENTIAL ABUNDANCE TESTING WITH ANCOM-BC

#Perform this test at phylum level 
qiime taxa collapse \
--i-table table-filtered.qza \
--i-taxonomy taxonomy.qza \
--p-level 2 \
--o-collapsed-table table-l2.qza

qiime composition ancombc \
--i-table table-l2.qza \
--m-metadata-file endolith_metadata.tsv \
--p-formula 'sample_type' \
--o-differentials l2-ancombc-sample-type.qza

qiime composition da-barplot \
--i-data l2-ancombc-sample-type.qza \
--p-significance-threshold 0.001 \
--p-level-delimiter ';' \
--o-visualization l2-da-barplot-sample-type.qzv

#Perform this test at genus level 
qiime taxa collapse \
--i-table table-filtered.qza \
--i-taxonomy taxonomy.qza \
--p-level 6 \
--o-collapsed-table table-l6.qza

qiime composition ancombc \
--i-table table-l6.qza \
--m-metadata-file endolith_metadata.tsv \
--p-formula 'sample_type' \
--o-differentials l6-ancombc-sample-type.qza

qiime composition da-barplot \
--i-data l6-ancombc-sample-type.qza \
--p-significance-threshold 0.001 \
--p-level-delimiter ';' \
--o-visualization l6-da-barplot-sample-type.qzv

## EXPORTING DATA FOR R (IF USING PHYLOSEQ OR ANCOMBC) - microeco can convert qiime2 
## files and doesn't require this step

#export the ASV table into biom format, then convert to .tsv
qiime tools export \
--input-path table-filtered.qza \
--output-path export

biom convert \
-i export/feature-table.biom \
-o export/feature-table.tsv --to-tsv

#export representative sequences
qiime tools export \
--input-path rep-seqs-filtered.qza \
--output-path export

#export taxonomy table
qiime tools export \
--input-path taxonomy.qza \
--output-path export

#export phylogenetic tree
qiime tools export \
--input-path rooted-tree.qza \
--output-path export



