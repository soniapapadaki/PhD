###############################################
# Principal Coordinates Analysis (with elemental concentration data)
#
# Description:
#   PCoA analysis of 18S rRNA gene diversity with elemental concentrations as
#   loading arrows, using microeco and envfit function from vegan R package.
###############################################

### LOAD PACKAGES
library(file2meco)
library(microeco)
library(ggplot2)
library(ggpubr)
library(GUniFrac)
library(paletteer)
library(vegan)
library(ggrepel)
library(dplyr)
library(compositions)
library(tidyr)
library(patchwork)
library(tibble)

#### MICROBIOME DATA PREP ####
metadata_file_path   <- "endolith_metadata.tsv"
euk_abund_file_path <- "../18S/table-filtered.qza"
euk_taxonomy_file_path <- "../18S/taxonomy.qza"
euk_tree_data <- "../18S/rooted-tree.qza"
euk_rep_data <- "../18S/rep-seqs-filtered.qza"

# Load raw microeco object
euk_raw <- qiime2meco(
  euk_abund_file_path,
  sample_table = metadata_file_path,
  taxonomy_table = euk_taxonomy_file_path,
  phylo_tree = euk_tree_data,
  rep_fasta = euk_rep_data
)
euk_raw$tidy_dataset()
abund_raw <- euk_raw$otu_table
meta_raw  <- euk_raw$sample_table

################ PCoA WITH ALL SAMPLES ################
abund_t <- as.data.frame(t(abund_raw))
abund_t$sample_type_site <- meta_raw$sample_type_site

abund_avg <- abund_t %>%
  group_by(sample_type_site) %>%
  summarise(across(where(is.numeric), mean))

abund_avg_final <- as.data.frame(t(abund_avg[ , -1]))
colnames(abund_avg_final) <- abund_avg$sample_type_site

metadata_avg <- meta_raw %>%
  group_by(sample_type_site) %>%
  summarise(
    sample_type = first(sample_type),
    sample_site = first(sample_site),
    across(where(is.numeric), mean, na.rm = TRUE)
  ) %>%
  as.data.frame()
rownames(metadata_avg) <- metadata_avg$sample_type_site

euk_microdata_full <- microtable$new(
  otu_table = abund_avg_final,
  sample_table = metadata_avg,
  tax_table = euk_raw$tax_table,
  phylo_tree = euk_raw$phylo_tree,
  rep_fasta = euk_raw$rep_fasta
)
euk_microdata_full$tidy_dataset()

# Save averaged abundances
euk_microdata_full$cal_abund(rel = FALSE)
euk_microdata_full$save_abund(dirpath = 'euk_taxa_abund_full_raw_AVERAGED',
                              merge_all = TRUE, sep = '\t', quote = FALSE)
euk_microdata_full$cal_abund()
euk_microdata_full$save_abund(dirpath = 'euk_taxa_abund_full_AVERAGED',
                              merge_all = TRUE, sep = '\t', quote = FALSE)

# Beta diversity and PCoA
euk_microdata_full$cal_betadiv(unifrac = TRUE)
p_full <- trans_beta$new(dataset = euk_microdata_full, group = 'sample_type', 
                         measure = 'wei_unifrac')
p_full$cal_ordination(method = 'PCoA')
p_pcoa_full <- p_full$plot_ordination(
  plot_color = 'sample_site',
  plot_shape = 'sample_type',
  plot_type = 'point',
  point_size = 3,
  shape_values = c(19,17),
  color_values = paletteer_d("ggthemes::gdoc")
) +
  stat_ellipse(aes(group = sample_type), linetype = "dashed", color = "black", 
               linewidth = 0.5) +
  theme_bw() +
  theme(axis.title = element_text(size=12), legend.text = element_text(size=12))
p_pcoa_full

################ PCoA WITH ONLY ENDOLITH SAMPLES ################
endolith_samples <- meta_raw$sample_type == "Endolith"
meta_endo  <- meta_raw[endolith_samples, ]
abund_endo <- abund_raw[, rownames(meta_endo)]

abund_t <- as.data.frame(t(abund_endo))
abund_t$sample_type_site <- meta_endo$sample_type_site

abund_avg <- abund_t %>%
  group_by(sample_type_site) %>%
  summarise(across(where(is.numeric), mean))
abund_avg_final <- as.data.frame(t(abund_avg[,-1]))
colnames(abund_avg_final) <- abund_avg$sample_type_site

metadata_avg <- meta_endo %>%
  group_by(sample_type_site) %>%
  summarise(
    sample_type = first(sample_type),
    sample_site = first(sample_site),
    across(where(is.numeric), mean, na.rm = TRUE)
  ) %>%
  as.data.frame()
rownames(metadata_avg) <- metadata_avg$sample_type_site

euk_microdata_endo <- microtable$new(
  otu_table = abund_avg_final,
  sample_table = metadata_avg,
  tax_table = euk_raw$tax_table,
  phylo_tree = euk_raw$phylo_tree,
  rep_fasta = euk_raw$rep_fasta
)
euk_microdata_endo$tidy_dataset()

# Save averaged abundances
euk_microdata_endo$cal_abund(rel = FALSE)
euk_microdata_endo$save_abund(dirpath = 'ENDO_euk_taxa_abund_raw_AVERAGED',
                              merge_all = TRUE, sep = '\t', quote = FALSE)
euk_microdata_endo$cal_abund()
euk_microdata_endo$save_abund(dirpath = 'ENDO_euk_taxa_abund_AVERAGED',
                              merge_all = TRUE, sep = '\t', quote = FALSE)

# Beta diversity and PCoA
euk_microdata_endo$cal_betadiv(unifrac = TRUE)
p_endo <- trans_beta$new(dataset = euk_microdata_endo, group = 'sample_site', 
                         measure = 'wei_unifrac')
p_endo$cal_ordination(method = 'PCoA')
p_pcoa_endo <- p_endo$plot_ordination(
  plot_color = 'sample_site',
  plot_shape = 'sample_type',
  plot_type = 'point',
  point_size = 3,
  shape_values = c(19,17),
  color_values = paletteer_d("ggthemes::gdoc")
) +
  stat_ellipse(aes(group = sample_type), linetype = "dashed", color = "black", 
               linewidth = 0.5) +
  theme_bw() +
  theme(axis.title = element_text(size=12), legend.text = element_text(size=12))
p_pcoa_endo

################ ADD ELEMENT CONCENTRATION DATA ################

#### ELEMENTAL DATA PREP (combined workflow) ####

# Load elemental data
ea_raw <- read.delim("EA_PCOA_concentrations.csv", stringsAsFactors = FALSE)
ea_loq <- read.delim("EA_PCOA_LOQ.csv", stringsAsFactors = FALSE)
ea_meta <- read.delim("EA_PCOA_metadata.csv", stringsAsFactors = FALSE)

# Convert all values to mg/kg
ea_raw <- ea_raw %>%
  mutate(value_mgkg = case_when(
    unit == "wt%"  ~ value_raw * 10000,   # wt% to mg/kg
    unit == "mg/kg" ~ value_raw,          # already mg/kg
    TRUE ~ NA_real_
  ))

# Replace censored values (BDL) with LOQ/sqrt(2)
ea_raw <- ea_raw %>%
  mutate(value_clean = case_when(
    BDL == "<" ~ value_mgkg / sqrt(2),
    TRUE       ~ value_mgkg
  ))

# Merge LOQ + metadata
ea_raw2 <- ea_raw %>%
  left_join(ea_loq, by = c("element" = "Element")) %>%
  left_join(ea_meta, by = "sampleid")

# Pivot table to wide format
ea_wide <- ea_raw2 %>%
  select(sampleid, element, value_clean) %>%
  pivot_wider(names_from = element, values_from = value_clean)

# Prepare metadata
ea_meta2 <- ea_meta %>%
  mutate(
    sample_type_abbrev = case_when(
      sample_type == "Endolith"   ~ "E",
      sample_type == "Soil crust" ~ "S",
      TRUE ~ NA_character_
    ),
    sample_type_site = paste0(sample_site, "-", sample_type_abbrev)
  )

# Merge wide elemental data with metadata
ea_wide2 <- ea_wide %>%
  left_join(ea_meta2 %>% select(sampleid, sample_type_site, sample_type),
            by = "sampleid")


### Function to run envfit for a given subset of samples
run_envfit <- function(pcoa_scores, ea_wide2_subset, output_prefix) {
  
  # Make matrix of PCoA site scores
  pc_mat <- as.matrix(pcoa_scores)
  
  # Prepare elemental matrix
  ea_mat <- ea_wide2_subset %>%
    column_to_rownames("sample_type_site") %>%
    select(-sampleid)
  
  # Ensure matching row order
  ea_mat <- ea_mat[rownames(pc_mat), ]
  
  # Convert first 3 columns to numeric
  pc_num <- as.data.frame(apply(pc_mat[,1:3], 2, function(x) as.numeric(trimws(x))))
  
  # Fit envfit
  fit <- envfit(pc_num ~ ., data = ea_mat)
  
  # Save table
  scores  <- as.data.frame(fit$vectors$arrows)
  r_vals  <- as.numeric(fit$vectors$r)
  p_vals  <- as.numeric(fit$vectors$pvals)
  
  envfit_table <- data.frame(
    Variable = rownames(scores),
    PCo1 = scores[,1],
    PCo2 = scores[,2],
    r2 = r_vals,
    p_value = p_vals,
    stringsAsFactors = FALSE
  )
  
  envfit_table[, -1] <- lapply(envfit_table[, -1], function(x) as.numeric(as.character(x)))
  envfit_table_pp <- envfit_table
  envfit_table_pp[, -1] <- lapply(envfit_table_pp[, -1], function(x) round(x, 3))
  
  write.csv(envfit_table_pp, paste0(output_prefix, "_envfit_results.csv"), row.names = FALSE)
  
  # Plot arrows
  fit_df <- data.frame(element = rownames(fit$vectors$arrows), r2 = fit$vectors$r)
  
  sig_vars <- names(fit$vectors$pvals)[fit$vectors$pvals < 0.05]
  top10_vars <- fit_df %>% arrange(desc(r2)) %>% slice(1:10) %>% pull(element)
  vars_to_keep <- intersect(sig_vars, top10_vars)
  
  arrow_df <- as.data.frame(scores(fit, "vectors"))
  arrow_df$element <- rownames(arrow_df)
  colnames(arrow_df)[1:2] <- c("PCo1","PCo2")
  arrow_df <- arrow_df %>% filter(element %in% vars_to_keep)
  
  max_pcoa <- max(abs(c(pcoa_scores[,1], pcoa_scores[,2])))
  max_arrow <- max(abs(c(arrow_df$PCo1, arrow_df$PCo2)))
  arrow_scaler <- (0.2 * max_pcoa) / max_arrow
  
  p_pcoa_with_arrows <- p_pcoa_final +  # Use whichever p_pcoa object corresponds
    geom_segment(
      data = arrow_df,
      inherit.aes = FALSE,
      aes(x=0, y=0, xend=PCo1*arrow_scaler, yend=PCo2*arrow_scaler),
      arrow = arrow(length = unit(0.25, "cm")),
      color = "red",
      linewidth = 0.8
    ) +
    geom_text_repel(
      data = arrow_df,
      inherit.aes = FALSE,
      aes(x=PCo1*arrow_scaler, y=PCo2*arrow_scaler, label=element),
      color="red",
      size=3,
      max.overlaps = Inf
    )
  
  return(p_pcoa_with_arrows)
}


### Run envfit on full dataset (all samples)
p_pcoa_full_with_arrows <- run_envfit(pcoa_scores = p3$res_ordination$scores,
                                      ea_wide2_subset = ea_wide2,
                                      output_prefix = "euk")

### Run envfit on only endolith samples
ea_wide2_endo <- ea_wide2 %>% filter(sample_type == "Endolith")
p_pcoa_endo_with_arrows <- run_envfit(pcoa_scores = p3$res_ordination$scores,
                                      ea_wide2_subset = ea_wide2_endo,
                                      output_prefix = "endo_euk")