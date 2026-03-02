###############################################
# MICROECO 
#
# Description:
#   Microbiome data analysis in microeco R package using outputs from QIIME. 
#   Produces relative abundance plots (phylum and family level), alpha diversity 
#   box plots, and beta diversity plots (PCoA).
#
# Notes:
#   This script is adapted from sections 4.1, 5.1, and 5.2 in the microeco 
#   R package tutorial:
#   https://chiliubio.github.io/microeco_tutorial/
###############################################

### LOAD PACKAGES

library(file2meco)
library(microeco)
library(ggplot2)
library(patchwork) #for combining plots together in grid
library(ggh4x) #for plot facets
library(ggpubr)
library(paletteer) #for ggthemes plot colour palette
library(GUniFrac) #if using UniFrac for beta diversity analysis
library(RColorBrewer)


#### RELATIVE ABUNDANCE PLOTS ####

## PROKARYOTES

# import data
prok_abund_file_path <- "../16S/table-filtered.qza" #Count table
metadata_file_path <- "endolith_metadata.tsv" #Metadata file
prok_taxonomy_file_path <- "../16S/taxonomy-v1.qza" #Taxonomy table
prok_tree_data <- "../16S/rooted-tree.qza" #Phylogenetic tree
prok_rep_data <- "../16S/rep-seqs-filtered.qza" #Sequences

# construct microtable object
qiime2meco(prok_abund_file_path)
prok_microdata <- qiime2meco(prok_abund_file_path, 
                             sample_table = metadata_file_path, 
                             taxonomy_table = prok_taxonomy_file_path, 
                             phylo_tree = prok_tree_data, 
                             rep_fasta = prok_rep_data)
prok_microdata

# make sample info consistent across all files in the dataset object
prok_microdata$tidy_dataset()
print(prok_microdata)

# save raw taxa abundance at each rank as .tsv file
prok_microdata$cal_abund(rel = FALSE)
prok_microdata$save_abund(dirpath = 'prok_taxa_abund_raw')
prok_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

# save relative taxa abundance at each rank as .tsv file
prok_microdata$cal_abund()
prok_microdata$save_abund(dirpath = 'taxa_abund')
prok_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

# relative abundance of top 8 phyla between sample types and sites
p1 <- trans_abund$new(dataset = prok_microdata, taxrank = 'Phylum', ntaxa = 8)
prok_plot <- p1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
            others_color = 'grey70', facet = c('sample_type','sample_site'), 
            xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Phylum')) & 
  theme(legend.key.size = unit(10, "pt"),
        legend.justification = "right",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "pt")) 


## EUKARYOTES

euk_abund_file_path <- "../18S/table-filtered.qza" #Count table
euk_taxonomy_file_path <- "../18S/taxonomy.qza" #Taxonomy table
euk_tree_data <- "../18S/rooted-tree.qza" #Phylogenetic tree
euk_rep_data <- "../18S/rep-seqs-filtered.qza" #Sequences

qiime2meco(euk_abund_file_path)
euk_microdata <- qiime2meco(euk_abund_file_path, 
                            sample_table = metadata_file_path, 
                            taxonomy_table = euk_taxonomy_file_path, 
                            phylo_tree = euk_tree_data, 
                            rep_fasta = euk_rep_data)
euk_microdata

euk_microdata$tidy_dataset()
print(euk_microdata)

euk_microdata$cal_abund(rel = FALSE)
euk_microdata$save_abund(dirpath = 'euk_taxa_abund_raw')
euk_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

euk_microdata$cal_abund()
euk_microdata$save_abund(dirpath = 'euk_taxa_abund')
euk_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

e1 <- trans_abund$new(dataset = euk_microdata, taxrank = 'Phylum', ntaxa = 8)
euk_plot <- e1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"),
            others_color = 'grey70', facet = c('sample_type','sample_site'), 
            xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Phylum')) & 
  theme(legend.key.size = unit(10, "pt"),
        legend.justification = "right",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "pt"))

## CYANOBACTERIA ONLY

cy_abund_file_path <- "../16S/table-filtered-only-cyano.qza" 
cy_rep_data <- "../16S/rep-seqs-filtered-only-cyano.qza" 

qiime2meco(cy_abund_file_path)
cy_microdata <- qiime2meco(cy_abund_file_path, sample_table = metadata_file_path, 
                           taxonomy_table = prok_taxonomy_file_path, 
                           rep_fasta = cy_rep_data)
cy_microdata

cy_microdata$tidy_dataset()
print(cy_microdata)

cy_microdata$cal_abund()
cy_microdata$save_abund(dirpath = 'cy_taxa_abund')
cy_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

cy1 <- trans_abund$new(dataset = cy_microdata, taxrank = 'Family', ntaxa = 8)
cy_plot <- cy1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
                        others_color = 'grey70', 
                        facet = c('sample_type','sample_site'), 
                        xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Cyanobacteria')) 


## PROTEOBACTERIA ONLY

prot_abund_file_path <- "../16S/table-filtered-only-pro.qza" 
prot_rep_data <- "../16S/rep-seqs-filtered-only-pro.qza"

qiime2meco(prot_abund_file_path)
prot_microdata <- qiime2meco(prot_abund_file_path, sample_table = metadata_file_path, taxonomy_table = prok_taxonomy_file_path, rep_fasta = prot_rep_data)
prot_microdata

prot_microdata$tidy_dataset()
print(prot_microdata)

prot_microdata$cal_abund()
prot_microdata$save_abund(dirpath = 'proteo_taxa_abund')
prot_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

prot1 <- trans_abund$new(dataset = prot_microdata, taxrank = 'Family', ntaxa = 8)
prot_plot <- prot1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
                            others_color = 'grey70', 
                            facet = c('sample_type','sample_site'), 
                            xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Proteobacteria')) 


## ACTINOBACTERIOTA ONLY

act_abund_file_path <- "../16S/table-filtered-only-act.qza" 
act_rep_data <- "../16S/rep-seqs-filtered-only-act.qza" 

qiime2meco(act_abund_file_path)
act_microdata <- qiime2meco(act_abund_file_path, sample_table = metadata_file_path, 
                            taxonomy_table = prok_taxonomy_file_path, 
                            rep_fasta = act_rep_data)
act_microdata

act_microdata$tidy_dataset()
print(act_microdata)

act_microdata$cal_abund()
act_microdata$save_abund(dirpath = 'actino_taxa_abund')
act_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

act1 <- trans_abund$new(dataset = act_microdata, taxrank = 'Family', ntaxa = 8)
actino_plot <- act1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
                             others_color = 'grey70', 
                             facet = c('sample_type','sample_site'), 
                             xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Actinobacteriota')) 


## BACTERIODOTA ONLY

bda_abund_file_path <- "../16S/table-filtered-only-bda.qza" 
bda_rep_data <- "../16S/rep-seqs-filtered-only-bda.qza" 

qiime2meco(bda_abund_file_path)
bda_microdata <- qiime2meco(bda_abund_file_path, sample_table = metadata_file_path, 
                            taxonomy_table = prok25_taxonomy_file_path, 
                            rep_fasta = bda_rep_data)
bda_microdata

bda_microdata$tidy_dataset()
print(bda_microdata)

bda_microdata$cal_abund()
bda_microdata$save_abund(dirpath = 'bacteroidota_taxa_abund')
bda_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

bda1 <- trans_abund$new(dataset = bda_microdata, taxrank = 'Family', ntaxa = 8)
bda_plot <- bda1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
                         others_color = 'grey70', 
                         facet = c('sample_type','sample_site'), 
                         xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Bacteroidota'))


### ASCOMYCOTA ONLY

asc_abund_file_path <- "../18S/table-filtered-only-asc.qza" 
asc_rep_data <- "../18S/rep-seqs-filtered-only-asc.qza" 

qiime2meco(asc_abund_file_path)
asc_microdata <- qiime2meco(asc_abund_file_path, 
                            sample_table = metadata_file_path, 
                            taxonomy_table = euk_taxonomy_file_path, 
                            rep_fasta = asc_rep_data)
asc_microdata

asc_microdata$tidy_dataset()
print(asc_microdata)

asc_microdata$cal_abund()
asc_microdata$save_abund(dirpath = 'asco_taxa_abund')
asc_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

asc1 <- trans_abund$new(dataset = asc_microdata, taxrank = 'Family', ntaxa = 8)
asco_plot <- asc1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
                           others_color = 'grey70', 
                           facet = c('sample_type','sample_site'), 
                           xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Ascomycota')) 


## CHLOROPHYTA ONLY

cpt_abund_file_path <- "../18S/table-filtered-only-cpt.qza" 
cpt_rep_data <- "../18S/rep-seqs-filtered-only-cpt.qza" 

qiime2meco(cpt_abund_file_path)
cpt_microdata <- qiime2meco(cpt_abund_file_path, 
                            sample_table = metadata_file_path, 
                            taxonomy_table = euk_taxonomy_file_path, 
                            rep_fasta = cpt_rep_data)
cpt_microdata

cpt_microdata$tidy_dataset()
print(cpt_microdata)

cpt_microdata$cal_abund()
cpt_microdata$save_abund(dirpath = 'chlorophyt_taxa_abund')
cpt_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

cpt1 <- trans_abund$new(dataset = cpt_microdata, taxrank = 'Family')
cpt_plot <- cpt1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
                          others_color = 'grey70', 
                          facet = c('sample_type','sample_site'), 
                          xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Chlorophyta')) 


### CILIOPHORA ONLY

cil_abund_file_path <- "../18S/table-filtered-only-cil.qza" 
cil_rep_data <- "../18S/rep-seqs-filtered-only-cil.qza" 

qiime2meco(cil_abund_file_path)
cil_microdata <- qiime2meco(cil_abund_file_path, 
                            sample_table = metadata_file_path, 
                            taxonomy_table = euk_taxonomy_file_path, 
                            rep_fasta = cil_rep_data)
cil_microdata

cil_microdata$tidy_dataset()
print(cil_microdata)

cil_microdata$cal_abund()
cil_microdata$save_abund(dirpath = 'cilio_taxa_abund')
cil_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

cil1 <- trans_abund$new(dataset = cil_microdata, taxrank = 'Family', ntaxa = 8)
cil_plot <- cil1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
                          others_color = 'grey70', 
                          facet = c('sample_type','sample_site'), 
                          xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Ciliophora')) 


### BASIDIOMYCOTA ONLY

bas_abund_file_path <- "../18S/table-filtered-only-bas.qza" 
bas_rep_data <- "../18S/rep-seqs-filtered-only-bas.qza" 

qiime2meco(bas_abund_file_path)
bas_microdata <- qiime2meco(bas_abund_file_path, 
                            sample_table = metadata_file_path, 
                            taxonomy_table = euk25_taxonomy_file_path, 
                            rep_fasta = bas_rep_data)
bas_microdata

bas_microdata$tidy_dataset()
print(bas_microdata)

bas_microdata$cal_abund()
bas_microdata$save_abund(dirpath = 'basidio_taxa_abund')
bas_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

# Create bar plot comparing top 8 phyla between sample types and sites
bas1 <- trans_abund$new(dataset = bas_microdata, taxrank = 'Family', ntaxa = 8)
bas_plot <- bas1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
                          others_color = 'grey70', 
                          facet = c('sample_type','sample_site'), 
                          xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Basidiomycota')) 


## COMPILING FAMILY PLOTS TOGETHER

# 16S
prok_plot / prot_plot2 / actino_plot / cy_plot / bda_plot + 
  plot_layout(axes = 'collect_y') & 
  theme(legend.key.size = unit(10, "pt"),
        legend.justification = "right",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "pt"))

# 18S
euk_plot / asco_plot / cpt_plot / cil_plot / bas_plot + 
  plot_layout(axes = 'collect_y') & 
  theme(legend.key.size = unit(10, "pt"),
        legend.justification = "right",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "pt"))



#### ALPHA DIVERSITY PLOTS ####

# Calculate alpha diversity (NULL indicating to calculate all alpha diversity 
# measures) and save to a .csv file

prok_microdata$cal_alphadiv(measures = NULL, PD = FALSE)
prok_microdata$save_alphadiv(dirpath = 'prok_alpha_diversity')
euk_microdata$cal_alphadiv(measures = NULL, PD = FALSE)
euk_microdata$save_alphadiv(dirpath = 'euk_alpha_diversity')

# Create boxplots comparing soil crust and endolith alpha diversity

p2 <- trans_alpha$new(dataset = prok_microdata, group = "sample_type")
p2$data_stat[which(p2$data_stat$Mean == min(p2$data_stat$Mean[p2$data_stat$Measure=='Shannon'])),]
p2$data_stat[which(p2$data_stat$Mean == max(p2$data_stat$Mean[p2$data_stat$Measure=='Shannon'])),]
p2$data_stat[which(p2$data_stat$Mean == min(p2$data_stat$Mean[p2$data_stat$Measure=='InvSimpson'])),]
p2$data_stat[which(p2$data_stat$Mean == max(p2$data_stat$Mean[p2$data_stat$Measure=='InvSimpson'])),]
p2$cal_diff(method = 'wilcox')

e2 <- trans_alpha$new(dataset = euk_microdata, group = "sample_type")
e2$data_stat[which(e2$data_stat$Mean == min(e2$data_stat$Mean[e2$data_stat$Measure=='Shannon'])),]
e2$data_stat[which(e2$data_stat$Mean == max(e2$data_stat$Mean[e2$data_stat$Measure=='Shannon'])),]
e2$data_stat[which(e2$data_stat$Mean == min(e2$data_stat$Mean[e2$data_stat$Measure=='InvSimpson'])),]
e2$data_stat[which(e2$data_stat$Mean == max(e2$data_stat$Mean[e2$data_stat$Measure=='InvSimpson'])),]
e2$cal_diff(method = 'wilcox')

p_Obs <- p2$plot_alpha(measure = 'Observed', add = "jitter", xtext_angle = 0, 
                       add_sig = FALSE) + 
  facet_wrap2(~Measure) + 
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 10, ),axis.text.x = element_blank(), 
        axis.ticks.x.bottom = element_blank(), 
        panel.border = element_rect(colour = 1, fill = NA, linewidth = 1)) +
  ylab("Alpha Diversity Measure") +
  guides(color = guide_legend(title = "Habitat"))

e_Obs <- e2$plot_alpha(measure = 'Observed', add = "jitter", xtext_angle = 0, 
                       add_sig_text_size = 5, y_start = 0.05) + 
  facet_wrap2(~Measure) + ylim (0, 220) + 
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 10, ),axis.text.x = element_blank(), 
        axis.ticks.x.bottom = element_blank(), 
        panel.border = element_rect(colour = 1, fill = NA, linewidth = 1)) +
  ylab("Alpha Diversity Measure") +
  guides(color = guide_legend(title = "Habitat"))

p_Sh <- p2$plot_alpha(measure = 'Shannon', add = "jitter", xtext_angle = 0, 
                      add_sig = FALSE) + 
  facet_wrap2(~Measure) + 
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 10, ),axis.text.x = element_blank(), 
        axis.ticks.x.bottom = element_blank(), 
        panel.border = element_rect(colour = 1, fill = NA, linewidth = 1)) +
  ylab("Alpha Diversity Measure") +
  guides(color = guide_legend(title = "Habitat"))

e_Sh <- e2$plot_alpha(measure = 'Shannon', add = "jitter", xtext_angle = 0, 
                      add_sig = FALSE) + 
  facet_wrap2(~Measure) + 
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 10, ),axis.text.x = element_blank(), 
        axis.ticks.x.bottom = element_blank(), 
        panel.border = element_rect(colour = 1, fill = NA, linewidth = 1)) +
  ylab("Alpha Diversity Measure") +
  guides(color = guide_legend(title = "Habitat"))

p_Inv <- p2$plot_alpha(measure = 'InvSimpson', add = "jitter", xtext_angle = 0, 
                       add_sig_text_size = 5, y_start = 0.05) + 
  facet_wrap2(~Measure) + ylim (0, 1250) +
  theme(axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 10), axis.text.x = element_blank(), 
        axis.ticks.x.bottom = element_blank(), 
        panel.border = element_rect(colour = 1, fill = NA, linewidth = 1)) +
  ylab("Alpha Diversity Measure") +
  guides(color = guide_legend(title = "Habitat"))

e_Inv <- e2$plot_alpha(measure = 'InvSimpson', add = "jitter", xtext_angle = 0, 
                       add_sig_text_size = 5, y_start = 0.05) + 
  facet_wrap2(~Measure) + ylim (0, 38) +
  theme(axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 10), axis.text.x = element_blank(), 
        axis.ticks.x.bottom = element_blank(), 
        panel.border = element_rect(colour = 1, fill = NA, linewidth = 1)) +
  ylab("Alpha Diversity Measure") +
  guides(color = guide_legend(title = "Habitat"))

# Combine boxplots into single image with shared legend
# export image as 875 x 437 pixels

combined_alpha_p <- p_Obs + p_Sh + p_Inv & 
  theme(legend.position = "right")

final_alpha_p <- combined_alpha_p + 
  plot_layout(axis_titles = "collect", guides = "collect")

ggsave("combined_alpha_prokaryotes.png",
       final_alpha_p,
       width = 875/100, height = 437/100, dpi = 100)


combined_alpha_e <- e_Obs + e_Sh + e_Inv & 
  theme(legend.position = "right")

final_alpha_e <- combined_alpha_e + 
  plot_layout(axis_titles = "collect", guides = "collect")

ggsave("combined_alpha_eukaryotes.png",
       final_alpha_e,
       width = 875/100, height = 437/100, dpi = 100)

# Save alpha diversity stats as .csv

write.csv(p2[["data_stat"]], "../16S/alpha_diversity/alpha_div_stats.csv")
write.csv(p2[["res_diff"]], "../16S/alpha_diversity/alpha_div_p.csv")
write.csv(e2[["data_stat"]], "../18S/alpha_diversity/alpha_div_stats.csv")
write.csv(e2[["res_diff"]], "../18S/alpha_diversity/alpha_div_p.csv")


### BETA DIVERSITY PLOTS

# Calculate beta diversity with unifrac and save to a .csv file

prok_microdata$cal_betadiv(unifrac = TRUE)
prok_microdata$save_betadiv(dirpath = "16S/beta_diversity")

euk_microdata$cal_betadiv(unifrac = TRUE)
euk_microdata$save_betadiv(dirpath = "18S/beta_diversity")

# Create 16S PCoA plot using UniFrac

p3 <- trans_beta$new(dataset = prok_microdata, group = 'sample_type', 
                     measure = 'wei_unifrac')
p3$cal_ordination(method = 'PCoA')

p3$cal_manova(manova_all = FALSE) #run permanova for comparison of two groups
p3$res_manova
p3$cal_betadisper() # run permdisp to assess group variances
p3$res_betadisper

p_pcoa_base <- p3$plot_ordination(plot_color = 'sample_site', 
                                  plot_shape = 'sample_type',
                                  plot_type = c('point'), point_size = 3, 
                                  shape_values = c(19, 17), 
                                  color_values = paletteer_d("ggthemes::gdoc"))

# Add ellipses, annotation, and theme
p_pcoa_final <- p_pcoa_base +
  stat_ellipse(aes(group = sample_type), linetype = "dashed", color = "black", 
               linewidth = 0.5) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        plot.subtitle = element_text(hjust = 0, 
        size = 10, face = "italic", margin = margin(b = 2))) +
  labs(subtitle = paste0("PERMANOVA R² Stat = ", round(p3$res_manova$R2, 2), 
                         ", p-value = ", p3$res_manova$p.value))

p_pcoa_final

# Create 18S PCoA plot using UniFrac
e3 <- trans_beta$new(dataset = euk_microdata, group = 'sample_type', 
                     measure = 'wei_unifrac')
e3$cal_ordination(method = 'PCoA')
e3$cal_manova(manova_all = FALSE) #run permanova for comparison of two groups
e3$res_manova
e3$cal_betadisper() # run permdisp to assess group variances
e3$res_betadisper

e_pcoa_base <- e3$plot_ordination(plot_color = 'sample_site', 
                                  plot_shape = 'sample_type',
                                  plot_type = c('point'), point_size = 3, 
                                  shape_values = c(19, 17), 
                                  color_values = paletteer_d("ggthemes::gdoc"))

# Add ellipses, annotation, and theme
e_pcoa_final <- e_pcoa_base +
  stat_ellipse(aes(group = sample_type), linetype = "dashed", color = "black", 
               linewidth = 0.5) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        plot.subtitle = element_text(hjust = 0, 
        size = 10, face = "italic", margin = margin(b = 2))) +
  labs(subtitle = paste0("PERMANOVA R² Stat = ", round(e3$res_manova$R2, 2), 
                         ", p-value = ", e3$res_manova$p.value))
e_pcoa_final

# combine PCoA plots and save as 1100 x 500 pixels
# convert pixels to inches (dpi = 100)

width_in <- 1100 / 100
height_in <- 500 / 100

combined_pcoa <- p_pcoa_final + e_pcoa_final & 
  theme(legend.position = "right")

final_pcoa <- combined_pcoa + 
  plot_layout(axis_titles = "collect", guides = "collect")

ggsave("combined_pcoa.png",
       final_pcoa,
       width = width_in,
       height = height_in,
       dpi = 100)

