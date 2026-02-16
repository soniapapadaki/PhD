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

physeq <- readRDS("16S_phyloflash.rds")

prok_microdata <- phyloseq2meco(physeq)

prok_microdata$tidy_dataset()
print(prok_microdata)

prok_microdata$cal_abund(rel = TRUE)
prok_microdata$save_abund(dirpath = 'prok_taxa_abund_raw')
prok_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

prok_microdata$cal_abund()
prok_microdata$save_abund(dirpath = 'prok_taxa_abund')
prok_microdata$save_abund(merge_all = TRUE, sep = '\t', quote = FALSE)

p1 <- trans_abund$new(dataset = prok_microdata, taxrank = 'Phylum', ntaxa = 10)
prok_plot <- p1$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"),
            others_color = 'grey70', facet = 'sample_type', 
            xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Phylum')) & 
  theme(legend.key.size = unit(10, "pt"),
        legend.justification = "right",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "pt"))

prok_plot


#### RELATIVE ABUNDANCE PLOT (FAMILY) ####

p2 <- trans_abund$new(dataset = prok_microdata, taxrank = 'Family', ntaxa = 15)
prok_plot_2 <- p2$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Dark2"),
                        others_color = 'grey70', facet = 'sample_type', 
                        xtext_keep = FALSE, legend_text_italic = FALSE) +
  guides(fill = guide_legend(ncol = 1, title = 'Family')) & 
  theme(legend.key.size = unit(10, "pt"),
        legend.justification = "right",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "pt"))

prok_plot_2


