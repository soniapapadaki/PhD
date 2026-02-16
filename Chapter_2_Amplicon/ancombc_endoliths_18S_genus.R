### SET WORKING DIRECTORY, CLEAR WORKSPACE AND LOAD PACKAGES

setwd("C:/Users/PC/Documents/R/ENDOLITHS_SC/ancombc/18S")
rm(list=ls())
graphics.off()

library(dplyr)
library(ANCOMBC)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(grid)
library(gtable)


### IMPORT PHYLOSEQ OBJECT

filename <- file.choose()
euk <- readRDS(filename)

# change "sample type" column into a factor object with levels "Endolith" and "Soil crust"
euk@sam_data[["sample_type"]] <- as.factor(euk@sam_data[["sample_type"]])

# change order of levels (putting "Soil crust" first so that this is used as the intercept
# instead of "Endolith")
euk@sam_data[["sample_type"]] = factor(euk@sam_data[["sample_type"]], levels = c("Soil crust", "Endolith"))

# verify
levels(euk@sam_data[["sample_type"]])


### PERFORM ANCOM-BC

set.seed(123)
output = ancombc2(data = euk, tax_level = "Genus",
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
df_fig_type_filtered = df_fig_type %>%
  dplyr::filter(`passed_ss_sample_typeEndolith` == 1) %>%
  dplyr::filter(!grepl('Eukaryota_', taxon) & !grepl('Unknown', taxon)) %>%
  dplyr::filter(!grepl('Helotiaceae', taxon) &
                  !grepl('MAST-12C', taxon) &
                  !grepl('Dactylopodida', taxon))


# plot 
fig_type = df_fig_type_filtered %>%
  ggplot(aes(x = `lfc_sample_typeEndolith`, y = taxon, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  scale_fill_manual(name = NULL, values = c("Enriched" = "#1b9e77", "Depleted" = "#d95f02")) +
  geom_errorbar(aes(xmin = `lfc_sample_typeEndolith` - `se_sample_typeEndolith`, 
                    xmax = `lfc_sample_typeEndolith` + `se_sample_typeEndolith`), 
                width = 0.2, position = position_dodge(0.05)) + 
  labs(y = NULL, x = "Log fold change", title = NULL) + 
  coord_cartesian(xlim = c(-3.3, 2.5)) +
  theme_bw(base_size = 12) + 
  theme(axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.key.size = unit(1.5, "lines"))

fig_type

# Set desired panel size
pixels_per_unit <- 65
x_range <- 2.5 - (-3.5)  # based on xlim setting
panel_width_inches <- (x_range * pixels_per_unit) / 100  # dpi = 100

# Set height per bar
dpi <- 100
height_per_bar <- 33
n_taxa <- length(unique(df_fig_type_filtered$taxon))
plot_height_inches <- (n_taxa * height_per_bar) / dpi

# Convert ggplot to grob (table layout)
g <- ggplotGrob(fig_type)

# Set panel width manually
panel_index <- g$layout[g$layout$name == "panel", ]
g$widths[panel_index$l] <- unit(panel_width_inches, "in")

# Save with ggsave using modified grob
ggsave("18S_genus_fixed_panel.png",
       g,
       width = sum(convertWidth(g$widths, "in", valueOnly = TRUE)),
       height = plot_height_inches,
       dpi = dpi)


########## OPTIONS #########

# If you are including taxa that failed the sensitivity analysis and want
# to colour code which ones failed/passed, add this to 'plot' section:
#theme(plot.title = element_text(hjust = 0.5),
      #panel.grid.minor.y = element_blank(),
      #axis.text.x = element_text(angle = 60, hjust = 1,
                                 #color = df_fig_type_filtered$color))

#To add p values to plot ("fig_type")
#p_formatted_c <- sapply(df_fig_type_filtered$p_sample_typeEndolith, function(x) {
#  if (x < 0.001) {
#    return(formatC(x, format = "e", digits = 3))  # Scientific notation for small values
#  } else {
#    return(formatC(x, format = "f", digits = 3))  # 3 decimal places for larger values
#  }
#})
#fig_type +
#  geom_text(
#    aes(x = ifelse(seq_along(taxon) %in% 9:11, 0.2, 0),  # Shift only specified row(s) to the right
#      y = taxon, 
#      label = paste("p =", p_formatted_c),  # Add "p =" to each label
#      hjust = ifelse(seq_along(taxon) %in% 9:11, 0, 1.2),  # Adjust hjust for specified row(s)
#    ),
#    vjust = 0.5,  # 0.5 to keep vertical centering
#    size = 2.5,   color = "gray36") 
