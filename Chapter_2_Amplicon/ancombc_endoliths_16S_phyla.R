### SET WORKING DIRECTORY, CLEAR WORKSPACE AND LOAD PACKAGES

setwd("C:/Users/PC/Documents/R/ENDOLITHS_SC/ancombc/16S")
rm(list=ls())
graphics.off()

library(dplyr)
library(ANCOMBC)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(ggpubr)


### IMPORT PHYLOSEQ OBJECT

filename <- file.choose()
prok <- readRDS(filename)

# change "sample type" column into a factor object with levels "Endolith" and "Soil crust"
prok@sam_data[["sample_type"]] <- as.factor(prok@sam_data[["sample_type"]])

# change order of levels (putting "Soil crust" first so that this is used as the intercept
# instead of "Endolith")
prok@sam_data[["sample_type"]] = factor(prok@sam_data[["sample_type"]], levels = c("Soil crust", "Endolith"))

# verify
levels(prok@sam_data[["sample_type"]])


### FILTER TAXA FOR ANCOMBC

# Initial run of ancombc came up with error requesting removal of two taxa, see 
# "structural zeros" in ancombc doc
taxa_to_remove <- c("eea54e724bb7511c47c5b7b56111878b", "31e8143a1cf9735a24bd5f6921752174", "5574b2c299eb1c01eb825a17d96aca57") 
prok <- prune_taxa(!taxa_names(prok) %in% taxa_to_remove, prok)
any(taxa_names(prok) %in% taxa_to_remove)  # Should return FALSE


### PERFORM ANCOM-BC

set.seed(123)
output = ancombc2(data = prok, tax_level = "Phylum",
                  fix_formula = "sample_type",
                  rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.1, lib_cut = 1000, s0_perc = 0.05,
                  group = "sample_type", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2),
                                       solver = "ECOS",
                                       B = 100))

###

res_prim = output$res

# from the output data, select only the data relevant for your plot, i.e. you don't want any
# of the intercept results, and you only want the measures for a particular variable like
# sample site or sample type, etc.

df_type = res_prim %>%
  dplyr::select(taxon, contains("sample_type")) 

# from your selected data, filter it down to only the taxa that pass the significant difference
# threshold, in your data table this is the "diff" column with TRUE (1) or FALSE (0)
df_fig_type = df_type %>%
  dplyr::filter(`diff_sample_typeEndolith` == 1) %>% 
  dplyr::arrange(desc(`lfc_sample_typeEndolith`)) %>%

# label >0 lfc values as "Enriched" and if not >0 then "Depleted"
# and change colour of taxon text if it passed the sensitivity analysis
  dplyr::mutate(direct = ifelse(`lfc_sample_typeEndolith` > 0, "Enriched", "Depleted"),
                color = ifelse(`passed_ss_sample_typeEndolith` == 1, "aquamarine3", "black"))

df_fig_type$taxon = factor(df_fig_type$taxon, levels = df_fig_type$taxon)
df_fig_type$direct = factor(df_fig_type$direct, 
                           levels = c("Enriched", "Depleted"))

# Remove taxa not identified at chosen level, and taxa that did not pass sensitivity analysis
#df_fig_type_filtered = df_fig_type %>%
  #dplyr::filter(!grepl('Bacteria_', taxon) & !grepl('Unknown', taxon) & !grepl('uncultured', taxon)) #%>%
  #dplyr::filter(`passed_ss_sample_typeEndolith` == 1)

# plot
fig_type = df_fig_type %>%
  ggplot(aes(x = `lfc_sample_typeEndolith`, y = taxon, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = `lfc_sample_typeEndolith` - `se_sample_typeEndolith`, xmax = `lfc_sample_typeEndolith` + `se_sample_typeEndolith`), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(y = "Prokaryotic Phyla", x = "Log fold change", 
       title = NULL) + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1,
                                   color = df_fig_type$color, size = 8))
# or, if you want to colour code taxon name by whether it passed sensitivity analysis:
#theme(plot.title = element_text(hjust = 0.5),
      #panel.grid.minor.y = element_blank(),
      #axis.text.x = element_text(angle = 60, hjust = 1,
                                 #color = df_fig_type_filtered$color))

#formatting p values and adding to plot
p_formatted_c <- sapply(df_fig_type$p_sample_typeEndolith, function(x) {
  if (x < 0.001) {
    return(formatC(x, format = "e", digits = 3))  # Scientific notation for small values
  } else {
    return(formatC(x, format = "f", digits = 3))  # 3 decimal places for larger values
  }
})

fig_type +
  geom_text(
    aes(
      x = ifelse(seq_along(taxon) %in% 13:16, 0.2, 0),  # Shift only row 13 to the right
      y = taxon, 
      label = paste("p =", p_formatted_c),  # Add "p =" to each label
      hjust = ifelse(seq_along(taxon) %in% 13:16, 0, 1.2),  # Adjust hjust for row 13
    ),
    vjust = 0.5,  # 0.5 to keep vertical centering
    size = 3,
    color = "gray36"
  ) +
  xlim(-2, 3)
