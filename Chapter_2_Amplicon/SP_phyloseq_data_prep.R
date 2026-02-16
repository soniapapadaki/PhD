### CLEAR WORKSPACE AND LOAD PACKAGES

rm(list=ls())
graphics.off()

library(tidyverse)
library(phyloseq)
library(readr)
library(seqinr)
library(vegan)
library(magrittr)
library(data.table)
library(ape) #appe::Ntip
library(dplyr)
library(RColorBrewer) #for plotting colours
library(tidyr) #wide to long format

### IMPORT DATA

metadata <- read.csv(file="endolith_metadata.csv", sep=",") #Metadata
metadata <- column_to_rownames(metadata, var = colnames(metadata)[1]) #make first col in our metadata (sample IDs) into a header

#16S rRNA Gene Amplicon Data
pro_count_tab <- read_tsv(file="C:/Users/PC/Documents/QIIME2/endoliths_sc_16S/export/feature-table.tsv", skip=1) #Count table
pro_count_tab <- column_to_rownames(pro_count_tab, var = colnames(pro_count_tab)[1]) #make first col in our count tab (ASV IDs) into a header
pro_count_tab <- pro_count_tab %>% relocate(SV103:SV126, .after = last_col()) #rearrange columns
pro_tax <- read_tsv(file="C:/Users/PC/Documents/QIIME2/endoliths_sc_16S/export/taxonomy.tsv") #Taxonomy table
pro_fasta <- read.fasta(file="C:/Users/PC/Documents/QIIME2/endoliths_sc_16S/export/dna-sequences.fasta") #Sequences
pro_tree <- read_tree("C:/Users/PC/Documents/QIIME2/endoliths_sc_16S/export/tree.nwk") #Phylogenetic tree

#18S rRNA Gene Amplicon Data
euk_count_tab <- read_tsv(file="C:/Users/PC/Documents/QIIME2/endoliths_sc_18S/export/feature-table.tsv", skip=1) #Count table
euk_count_tab <- column_to_rownames(euk_count_tab, var = colnames(euk_count_tab)[1]) #make first col in our count tab (ASV IDs) into a header
euk_count_tab <- euk_count_tab %>% relocate(SV103:SV126, .after = last_col()) #rearrange columns
euk_tax <- read_tsv(file="C:/Users/PC/Documents/QIIME2/endoliths_sc_18S/export/taxonomy.tsv") #Taxonomy table
euk_fasta <- read.fasta(file="C:/Users/PC/Documents/QIIME2/endoliths_sc_18S/export/dna-sequences.fasta") #Sequences
euk_tree <- read_tree("C:/Users/PC/Documents/QIIME2/endoliths_sc_18S/export/tree.nwk") #Phylogenetic tree

### FORMAT TAXONOMY TABLE (edit Qiime taxonomy output into easier format for building a phyloseq object)

#16S
pro_taxa_tab <- pro_tax %>%
  mutate(pro_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=";$", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=" ", replacement="")) %>%
  separate(pro_tax, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

#18S
euk_taxa_tab <- euk_tax %>%
  mutate(euk_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=";$", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=" ", replacement="")) %>%
  separate(euk_tax, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

###OUTPUT NUMBER OF SEQUENCES PER SAMPLE IN RAW READS

#16S
pro_count_tab_df = as.data.frame(colSums(pro_count_tab)) #make into dataframe
pro_count_tab_df$Sample = rownames(pro_count_tab_df) #make sample names a column
names(pro_count_tab_df) = c("Reads", "Sample")

pro_p<-ggplot(data=pro_count_tab_df, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pro_p

#18S
euk_count_tab_df = as.data.frame(colSums(euk_count_tab)) #make into df
euk_count_tab_df$Sample = rownames(euk_count_tab_df) #make sample names a column
names(euk_count_tab_df) = c("Reads", "Sample")

euk_p<-ggplot(data=euk_count_tab_df, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

euk_p

### RAW DATA BASIC STATS

#18S
euk_tot = sum(euk_count_tab_df$Reads) #total reads in the dataset
euk_asv = nrow(euk_count_tab) #number of unique ASVs

#16S
pro_tot = sum(pro_count_tab_df$Reads) #total reads in the dataset
pro_asv = nrow(pro_count_tab) #number of unique ASVs

### CONVERT DATA INTO PHYLOSEQ OBJECTS

#16S
pro_ASV = otu_table(as.matrix(pro_count_tab), taxa_are_rows = TRUE)
pro_TAX = tax_table(as.matrix(pro_taxa_tab))
pro_META = sample_data(data.frame(metadata, row.names = rownames(metadata)))
pro_TREE = pro_tree
pro_ps <- phyloseq(pro_ASV, pro_TAX, pro_META, pro_TREE) #make into phyloseq object

#18S
euk_ASV = otu_table(as.matrix(euk_count_tab), taxa_are_rows = TRUE) 
euk_TAX = tax_table(as.matrix(euk_taxa_tab))
euk_META = sample_data(data.frame(metadata, row.names = rownames(metadata)))
euk_TREE = euk_tree
euk_ps <- phyloseq(euk_ASV, euk_TAX, euk_META, euk_TREE) #make into phyloseq object

### CREATE VERSION OF PHYLOSEQ OBJECT WITH REMOVED SINGLETONS

#16S
pro_ps.filtered <- prune_taxa(taxa_sums(pro_ps) >= 10, pro_ps)

#18S
euk_ps.filtered <- prune_taxa(taxa_sums(euk_ps) >= 10, euk_ps)

### RE-ROOT TREE: Use this function > (https://john-quensen.com/r/unifrac-and-tree-roots/) < picks longest branch

#First define the function from link above to find furthest outgroup
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  #tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  #Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }

#Run on my phyloseq tree

#16S
pro_tree_new <- phy_tree(pro_ps.filtered)
pro_new_outgroup <- pick_new_outgroup(pro_tree_new)
pro_new_tree_root <- ape::root(pro_tree_new, outgroup=pro_new_outgroup, resolve.root=TRUE) # Re-root tree
pro_new_tree_dich <- ape::multi2di(pro_new_tree_root) # Convert to dichotomy tree
phy_tree(pro_ps.filtered) <- pro_new_tree_dich

#18S
euk_tree_new <- phy_tree(euk_ps.filtered)
euk_new_outgroup <- pick_new_outgroup(euk_tree_new)
euk_new_tree_root <- ape::root(euk_tree_new, outgroup=euk_new_outgroup, resolve.root=TRUE) # Re-root tree
euk_new_tree_dich <- ape::multi2di(euk_new_tree_root) # Convert to dichotomy tree
phy_tree(euk_ps.filtered) <- euk_new_tree_dich

### FILTERED DATA BASIC STATS

#16S
filt_pro_tot = sum(colSums(data.frame(otu_table(pro_ps.filtered))))
filt_pro_asv = nrow(data.frame(otu_table(pro_ps.filtered)))

#18S
filt_euk_tot = sum(colSums(data.frame(otu_table(euk_ps.filtered))))
filt_euk_asv = nrow(data.frame(otu_table(euk_ps.filtered)))

### SAVE PHYLOSEQ OBJECTS

saveRDS(pro_ps.filtered, "16S-phylo-object.rds")
saveRDS(euk_ps.filtered, "18S-phylo-object.rds")
