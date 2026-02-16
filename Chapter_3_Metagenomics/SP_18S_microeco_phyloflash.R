#### SET WORKING DIRECTORY, CLEAR WORKSPACE AND LOAD PACKAGES ####

library(file2meco)
library(microeco)
library(ggplot2)
library(patchwork) #for combining plots together in grid
library(ggh4x) #for plot facets
library(ggpubr)
#library(magrittr)
#library(paletteer) #for using different colour palettes in plots
library(GUniFrac) #if using UniFrac for beta diversity analysis
library(phyloseq)
#library(RColorBrewer)



#### RELATIVE ABUNDANCE PLOT (PHYLUM) ####

physeq <- readRDS("18S_phyloflash.rds")

euk_microdata <- phyloseq2meco(physeq)

euk_microdata$tidy_dataset()
print(euk_microdata)

euk_microdata$cal_abund(rel = TRUE)
euk_microdata$save_abund(dirpath = 'euk_taxa_abund_raw')
euk_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

euk_microdata$cal_abund()
euk_microdata$save_abund(dirpath = 'euk_taxa_abund')
euk_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

e1 <- trans_abund$new(dataset = euk_microdata, taxrank = 'Phylum', ntaxa = 10)
euk_plot <- e1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"),
            others_color = 'grey70', facet = 'sample_type', 
            xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Phylum')) & 
  theme(legend.key.size = unit(10, "pt"),
        legend.justification = "right",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "pt"))

euk_plot


#### RELATIVE ABUNDANCE PLOT (FAMILY) ####

e2 <- trans_abund$new(dataset = euk_microdata, taxrank = 'Family', ntaxa = 15)
euk_plot_2 <- e2$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"),
                        others_color = 'grey70', facet = 'sample_type', 
                        xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Family')) & 
  theme(legend.key.size = unit(10, "pt"),
        legend.justification = "right",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "pt"))

euk_plot_2


