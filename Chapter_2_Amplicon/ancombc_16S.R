###############################################
# ANCOM-BC (16S Family + Genus)
#
# Description:
#   Runs ANCOM-BC at family and genus levels, and generates custom plots for 
#   comparison of enriched taxa in endolith vs soil crust microbiomes.
#
# Notes:
#   This script is adapted from:
#   https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
###############################################

### LOAD PACKAGES

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
prok <- readRDS(filename)

# change "sample type" column into a factor object with levels "Endolith" and 
# "Soil crust"
prok@sam_data[["sample_type"]] <- as.factor(prok@sam_data[["sample_type"]])

# change order of levels (putting "Soil crust" first so that this is used as the 
# intercept instead of "Endolith")
prok@sam_data[["sample_type"]] = factor(prok@sam_data[["sample_type"]], 
                                        levels = c("Soil crust", "Endolith"))

# verify
levels(prok@sam_data[["sample_type"]])

### FILTER TAXA FOR ANCOMBC

# Initial run of ancombc came up with error requesting removal of two taxa, see 
# "structural zeros" in ancombc doc
taxa_to_remove <- c("eea54e724bb7511c47c5b7b56111878b", 
                    "31e8143a1cf9735a24bd5f6921752174", 
                    "5574b2c299eb1c01eb825a17d96aca57") 
prok <- prune_taxa(!taxa_names(prok) %in% taxa_to_remove, prok)
any(taxa_names(prok) %in% taxa_to_remove)  # Should return FALSE

### PERFORM ANCOM-BC (FAMILY LEVEL)

set.seed(123)
output_f = ancombc2(data = prok, tax_level = "Family",
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

res_prim_f = output_f$res

# from the output data, select only the data relevant for your plot, i.e. you 
# don't want any of the intercept results, and you only want the measures for a 
# particular variable like sample site or sample type, etc.

df_type_f = res_prim_f %>%
  dplyr::select(taxon, contains("sample_type")) 

# from your selected data, filter it down to only the taxa that pass the 
# significant difference threshold, in your data table this is the "diff" column 
# with TRUE (1) or FALSE (0)

df_fig_type_f = df_type_f %>%
  dplyr::filter(`diff_sample_typeEndolith` == 1) %>% 
  dplyr::arrange(desc(`lfc_sample_typeEndolith`)) %>%

# label >0 lfc values as "Enriched" and if not >0 then "Depleted"
# and change colour of taxon text if it passed the sensitivity analysis
  
  dplyr::mutate(direct = ifelse(`lfc_sample_typeEndolith` > 0, "Enriched", 
                                "Depleted"),
                color = ifelse(`passed_ss_sample_typeEndolith` == 1, 
                               "aquamarine3", "black"))

df_fig_type_f$taxon = factor(df_fig_type_f$taxon, levels = df_fig_type_f$taxon)
df_fig_type_f$direct = factor(df_fig_type_f$direct, 
                           levels = c("Enriched", "Depleted"))

# Remove taxa not identified at chosen level (in this case, Family), 
# and taxa that did not pass sensitivity analysis

df_fig_type_filtered_f = df_fig_type_f %>%
  dplyr::filter(`passed_ss_sample_typeEndolith` == 1) %>%
  dplyr::filter(!grepl('Bacteria_', taxon) & !grepl('Unknown', taxon)) %>%
  dplyr::filter(!grepl('PLTA13', taxon) &
                !grepl('DS-100', taxon) &
                !grepl('MB-A2-108', taxon) &
                !grepl('Rokubacteriales', taxon) &
                !grepl('Saccharimonadales', taxon) &
                !grepl('SAR324', taxon) &
                !grepl('C0119', taxon) &
                !grepl('AKAU4049', taxon) &
                !grepl('Gitt-GS-136', taxon)) 

# initial plot 
fig_type_f = df_fig_type_filtered_f %>%
  ggplot(aes(x = `lfc_sample_typeEndolith`, y = taxon, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  scale_fill_manual(name = NULL, values = c("Enriched" = "#1b9e77", 
                                            "Depleted" = "#d95f02")) +
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

fig_type_f

### PERFORM ANCOM-BC (GENUS LEVEL)

set.seed(123)
output_g = ancombc2(data = prok, tax_level = "Genus",
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

res_prim_g = output_g$res

df_type_g = res_prim_g %>%
  dplyr::select(taxon, contains("sample_type")) 

df_fig_type_g = df_type_g %>%
  dplyr::filter(`diff_sample_typeEndolith` == 1) %>% 
  dplyr::arrange(desc(`lfc_sample_typeEndolith`)) %>%
  dplyr::mutate(direct = ifelse(`lfc_sample_typeEndolith` > 0, "Enriched", 
                                "Depleted"),
                color = ifelse(`passed_ss_sample_typeEndolith` == 1, 
                               "aquamarine3", "black"))

df_fig_type_g$taxon = factor(df_fig_type_g$taxon, levels = df_fig_type_g$taxon)
df_fig_type_g$direct = factor(df_fig_type_g$direct, 
                            levels = c("Enriched", "Depleted"))

df_fig_type_filtered_g = df_fig_type_g %>%
  dplyr::filter(`passed_ss_sample_typeEndolith` == 1) %>%
  dplyr::filter(!grepl('Bacteria_', taxon) & !grepl('Unknown', taxon)) %>%
  dplyr::filter(!grepl('BIrii41', taxon) &
                  !grepl('PLTA13', taxon) &
                  !grepl('DS-100', taxon) &
                  !grepl('MB-A2-108', taxon) &
                  !grepl('Fimbriimonadaceae', taxon) &
                  !grepl('Rokubacteriales', taxon) &
                  !grepl('Saccharimonadales', taxon) &
                  !grepl('SAR324', taxon) &
                  !grepl('SM2D12', taxon) &
                  !grepl('C0119', taxon) &
                  !grepl('Obscuribacteraceae', taxon) &
                  !grepl('AKAU4049', taxon) &
                  !grepl('AKIW781', taxon) &
                  !grepl('Gitt-GS-136', taxon)) 

fig_type_g = df_fig_type_filtered_g %>%
  ggplot(aes(x = `lfc_sample_typeEndolith`, y = taxon, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  scale_fill_manual(name = NULL, values = c("Enriched" = "#1b9e77", 
                                            "Depleted" = "#d95f02")) +
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

fig_type_g


### FINAL PLOTS 

# shared plot settings:
dpi <- 100
x_range <- 2.5 - (-3.5)        
panel_width_inches <- x_range * 0.65

# Family-level plot:
n_taxa_f <- length(unique(df_fig_type_filtered_f$taxon))
plot_height_inches_f <- n_taxa_f * 0.32   # 32 px per bar at 100 dpi

g1 <- ggplotGrob(fig_type_f)
panel_index_f <- g1$layout[g1$layout$name == "panel", ]
g1$widths[panel_index_f$l] <- unit(panel_width_inches, "in")

ggsave(
  "ancombc_16S_family.png",
  g1,
  width = sum(convertWidth(g1$widths, "in", valueOnly = TRUE)),
  height = plot_height_inches_f,
  dpi = dpi
)

# Genus-level plot:
n_taxa_g <- length(unique(df_fig_type_filtered_g$taxon))
plot_height_inches_g <- n_taxa_g * 0.24   # 24 px per bar at 100 dpi

g2 <- ggplotGrob(fig_type_g)
panel_index_g <- g2$layout[g2$layout$name == "panel", ]
g2$widths[panel_index_g$l] <- unit(panel_width_inches, "in")

ggsave(
  "ancombc_16S_genus.png",
  g2,
  width = sum(convertWidth(g2$widths, "in", valueOnly = TRUE)),
  height = plot_height_inches_g,
  dpi = dpi
)

# OPTIONAL: Color-code taxa that failed/passed sensitivity analysis (if you 
# choose to not filter out taxa that failed). Add the following lines inside the 
# initial ggplot (fig_type) theme section:

# theme(
#   plot.title = element_text(hjust = 0.5),
#   panel.grid.minor.y = element_blank(),
#   axis.text.x = element_text(
#     angle = 60, hjust = 1,
#     color = df_fig_type_filtered$color
#   )
# )

# OPTIONAL: Add p-values to the plot (fig_type)
#
# p_formatted_c <- sapply(
#   df_fig_type_filtered$p_sample_typeEndolith,
#   function(x) {
#     if (x < 0.001) {
#       formatC(x, format = "e", digits = 3)  # scientific notation
#     } else {
#       formatC(x, format = "f", digits = 3)  # 3 decimal places
#     }
#   }
# )
#
# fig_type +
#   geom_text(
#     aes(
#       x     = ifelse(seq_along(taxon) %in% 9:11, 0.2, 0),  # shift rows 9–11 
                                                             # slightly right
#       y     = taxon,
#       label = paste("p =", p_formatted_c),
#       hjust = ifelse(seq_along(taxon) %in% 9:11, 0, 1.2)
#     ),
#     vjust = 0.5,
#     size  = 2.5,
#     color = "gray36"
#   )
