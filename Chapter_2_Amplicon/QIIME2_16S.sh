
#!/bin/bash

# Workflow for sequence quality control, ASV identification, and taxonomic classification

##IMPORTING DATA

#Create new directory
mkdir endoliths_sc_16S
cd endoliths_sc_16S

#Import sequencing data, metadata and manifest file
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest-file-16S-endoliths-sc.tsv \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

#Create summary file of sequence quality
qiime demux summarize \
--i-data paired-end-demux.qza \
--o-visualization paired-end-demux.qzv

##SEQUENCE QUALITY CONTROL AND FEATURE TABLE CONSTRUCTION

#DADA2 
qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trunc-len-f 190 \
--p-trunc-len-r 140 \
--p-trim-left-r 10 \
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

##TAXONOMIC CLASSIFICATION

#extract reference reads from SILVA 138 database
qiime feature-classifier extract-reads \
--i-sequences silva-138-99-seqs.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--o-reads ref-seqs.qza

#train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier 16S-classifier.qza

#use trained classifier to create taxonomy file for ASVs
qiime feature-classifier classify-sklearn \
--i-classifier 16S-classifier.qza \
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

##FILTER OUT MITOCHONDRIA, CHLOROPLASTS, EUKARYOTES, AND UNASSIGNED FOR 16S
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota \
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
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota \
  --p-include p__ \
  --o-filtered-sequences rep-seqs-filtered.qza

qiime feature-table summarize \
  --i-table table-filtered.qza \
  --o-visualization table-filtered.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-filtered.qza \
  --o-visualization rep-seqs-filtered.qzv
  
##PRODUCE RELATIVE ABUNDANCE PLOTS FOR TOP 4 PHYLA

#cyanobacteria
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota,d__Archaea \
  --p-include p__Cyanobacteria \
  --o-filtered-table table-filtered-only-cyano.qza
  
  qiime taxa barplot \
  --i-table table-filtered-only-cyano.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-only-cyano.qzv
  
qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota,d__Archaea \
  --p-include p__Cyanobacteria \
  --o-filtered-sequences rep-seqs-filtered-only-cyano.qza
  
 #proteobacteria
 qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota,d__Archaea \
  --p-include p__Proteobacteria \
  --o-filtered-table table-filtered-only-pro.qza
  
  qiime taxa barplot \
  --i-table table-filtered-only-pro.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-only-pro.qzv
  
 qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota,d__Archaea \
  --p-include p__Proteobacteria \
  --o-filtered-sequences rep-seqs-filtered-only-pro.qza
 
 #actinobacteriota
  qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota,d__Archaea \
  --p-include p__Actinobacteriota \
  --o-filtered-table table-filtered-only-act.qza
  
  qiime taxa barplot \
  --i-table table-filtered-only-act.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-only-act.qzv
  
 qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota,d__Archaea \
  --p-include p__Actinobacteriota \
  --o-filtered-sequences rep-seqs-filtered-only-act.qza
 
 #bacteroidota
  qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota,d__Archaea \
  --p-include p__Bacteroidota \
  --o-filtered-table table-filtered-only-bda.qza
  
  qiime taxa barplot \
  --i-table table-filtered-only-bda.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file endolith_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-only-bda.qzv
  
 qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,Unassigned,d__Eukaryota,d__Archaea \
  --p-include p__Bacteroidota \
  --o-filtered-sequences rep-seqs-filtered-only-bda.qza

##GENERATE PHYLOGENETIC TREE
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs-filtered.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

##ALPHA RAREFACTION - choosing max depth based on median in table-filtered.qzv
qiime diversity alpha-rarefaction \
--i-table table-filtered.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 108766 \
--p-steps 20 \
--m-metadata-file endolith_metadata.tsv \
--o-visualization alpha-rarefaction.qzv

qiime diversity alpha-rarefaction \
--i-table table-filtered.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 108766 \
--p-steps 20 \
--o-visualization alpha-rarefaction-eachsample.qzv

##ALPHA AND BETA DIVERSITY ANALYSIS

#Generate core diversity metrics, using a sampling depth based on alpha rarefaction curve
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table-filtered.qza \
--p-sampling-depth 10148 \
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

##DIFFERENTIAL ABUNDANCE TESTING WITH ANCOM-BC

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

##EXPORTING DATA FOR R (IF USING PHYLOSEQ) - microeco can convert qiime2 files and doesn't require this step

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
  