#### LOAD LIBRARIES ####
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

## --- PROKARYOTES --- ##
prok_abund_file_path <- "C:/Users/PC/Documents/QIIME2/ENDOLITHS_SC_2024/endoliths_sc_16S/table-filtered.qza"
metadata_file_path   <- "C:/Users/PC/Documents/R/ENDOLITHS_SC/microeco/endolith_metadata.tsv"
prok_taxonomy_file_path <- "C:/Users/PC/Documents/QIIME2/ENDOLITHS_SC_2024/endoliths_sc_16S/taxonomy-v1.qza"
prok_tree_data       <- "C:/Users/PC/Documents/QIIME2/ENDOLITHS_SC_2024/endoliths_sc_16S/rooted-tree.qza"
prok_rep_data        <- "C:/Users/PC/Documents/QIIME2/ENDOLITHS_SC_2024/endoliths_sc_16S/rep-seqs-filtered.qza"

# Load raw microeco object
prok_raw <- qiime2meco(
  prok_abund_file_path,
  sample_table = metadata_file_path,
  taxonomy_table = prok_taxonomy_file_path,
  phylo_tree = prok_tree_data,
  rep_fasta = prok_rep_data
)

prok_raw$tidy_dataset()

# Extract abundance and metadata
abund_raw <- prok_raw$otu_table    # ASV x sample
meta_raw  <- prok_raw$sample_table # includes "sample_type_site"

### --------------------------
### AVERAGE BY sample_type_site
### --------------------------

## 1. Average ASV abundances
abund_t <- as.data.frame(t(abund_raw))  # samples in rows
abund_t$sample_type_site <- meta_raw$sample_type_site

abund_avg <- abund_t %>%
  group_by(sample_type_site) %>%
  summarise(across(where(is.numeric), mean))

# Convert back to microeco-required format: ASV rows × samples
abund_avg_final <- as.data.frame(t(abund_avg[ , -1]))
colnames(abund_avg_final) <- abund_avg$sample_type_site

## 2. Average metadata (and keep sample_type or other descriptors)
metadata_avg <- meta_raw %>%
  group_by(sample_type_site) %>%
  summarise(
    sample_type = first(sample_type),
    sample_site = first(sample_site),      # <-- ADD THIS LINE
    across(where(is.numeric), mean, na.rm = TRUE)
  ) %>%
  as.data.frame()

rownames(metadata_avg) <- metadata_avg$sample_type_site


### --------------------------
### Build NEW microeco object
### --------------------------

prok_microdata <- microtable$new(
  otu_table = abund_avg_final,
  sample_table = metadata_avg,
  tax_table = prok_raw$tax_table,
  phylo_tree = prok_raw$phylo_tree,
  rep_fasta = prok_raw$rep_fasta
)

prok_microdata$tidy_dataset()

# You can still save the new averaged abundances for inspection
prok_microdata$cal_abund(rel = FALSE)
prok_microdata$save_abund(dirpath = 'prok_taxa_abund_raw_AVERAGED',
                          merge_all = TRUE, sep = '\t', quote = FALSE)

prok_microdata$cal_abund()
prok_microdata$save_abund(dirpath = 'prok_taxa_abund_AVERAGED',
                          merge_all = TRUE, sep = '\t', quote = FALSE)


#### PCOA PLOTS, AVERAGED BY SAMPLE SITE, (NO ELEMENTAL DATA) ####

# Calculate beta diversity on averaged data
prok_microdata$cal_betadiv(unifrac = TRUE)

# 16S PCoA using averaged microeco object
p3 <- trans_beta$new(
  dataset = prok_microdata,
  group = 'sample_type',        # your grouping variable
  measure = 'wei_unifrac'
)

p3$cal_ordination(method = 'PCoA')
p3$cal_manova(manova_all = FALSE)
p3$cal_betadisper()
p3$cal_anosim(group = "sample_type")

p_pcoa_final <- p3$plot_ordination(
  plot_color = 'sample_site',       # now valid
  plot_shape = 'sample_type',
  plot_type = 'point',
  point_size = 3,
  shape_values = c(19, 17),
  color_values = paletteer_d("ggthemes::gdoc")
) +
  stat_ellipse(
    aes(group = sample_type),
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.subtitle = element_text(
      hjust = 0,
      size = 10,
      face = "italic",
      margin = margin(b = 2)
    )
  ) +
  labs(
    subtitle = paste0(
      "ANOSIM R-Stat = ",
      round(p3$res_anosim$statistic, 3),
      ", p-value = ",
      p3$res_anosim$p.value
    )
  )
p_pcoa_final

#### ELEMENTAL DATA PREP ####

ea_raw <- read.delim("EA_PCOA_concentrations.csv", stringsAsFactors = FALSE)
ea_loq <- read.delim("EA_PCOA_LOQ.csv", stringsAsFactors = FALSE)
ea_meta <- read.delim("EA_PCOA_metadata.csv", stringsAsFactors = FALSE)

# ─────────────────────────────────────────
# 2. Convert all values to mg/kg
# ─────────────────────────────────────────
ea_raw <- ea_raw %>%
  mutate(value_mgkg = case_when(
    unit == "wt%"  ~ value_raw * 10000,   # convert wt% → mg/kg
    unit == "mg/kg" ~ value_raw,          # already mg/kg
    TRUE ~ NA_real_
  ))

# ─────────────────────────────────────────
# 3. Replace censored values (BDL) with LOQ/sqrt(2)
# ─────────────────────────────────────────
ea_raw <- ea_raw %>%
  mutate(value_clean = case_when(
    BDL == "<" ~ value_mgkg / sqrt(2),
    TRUE       ~ value_mgkg
  ))




# Merge LOQ + metadata
ea_raw2 <- ea_raw %>%
  left_join(ea_loq, by = c("element" = "Element")) %>%
  left_join(ea_meta, by = "sampleid")

# pivot table to wide format

ea_wide <- ea_raw2 %>%
  select(sampleid, element, value_clean) %>%
  pivot_wider(names_from = element, values_from = value_clean)

# match elemental samples to microbiome samples

ea_meta2 <- ea_meta %>%
  mutate(
    sample_type_abbrev = case_when(
      sample_type == "Endolith"    ~ "E",
      sample_type == "Soil crust"  ~ "S",
      TRUE ~ NA_character_
    ),
    sample_type_site = paste0(sample_site, "-", sample_type_abbrev)
  )

ea_wide2 <- ea_wide %>%
  left_join(ea_meta2 %>% select(sampleid, sample_type_site),
            by = "sampleid")

# Extract pcoa coordinates from earlier

pcoa_scores <- p3$res_ordination$scores

# Make a matrix of PCoA site scores
pc_mat <- as.matrix(pcoa_scores)

# prep ea_wide2 data into format for envfit

ea_mat <- ea_wide2 %>%
  column_to_rownames("sample_type_site") %>%
  select(-sampleid)

# put rownames in same order in both ea_mat and pc_mat
ea_mat <- ea_mat[rownames(pc_mat), ]

# prep pc_mat into format for envfit

# convert first 3 columns to numeric
pc_num <- as.data.frame(apply(pc_mat[, 1:3], 2, function(x) as.numeric(trimws(x))))

# check
str(pc_num)

# Fit environmental vectors with "envfit" in vegan
fit <- envfit(pc_num ~ ., data = ea_mat)
fit

# save fit scores as a table

scores  <- as.data.frame(fit$vectors$arrows)
r_vals  <- as.numeric(fit$vectors$r)
p_vals  <- as.numeric(fit$vectors$pvals)


envfit_table <- data.frame(
  Variable = rownames(scores),
  PCo1 = scores[, 1],
  PCo2 = scores[, 2],
  r2 = r_vals,
  p_value = p_vals,
  stringsAsFactors = FALSE
)

envfit_table[, -1] <- lapply(envfit_table[, -1], function(x) as.numeric(as.character(x)))

envfit_table_pp <- envfit_table
envfit_table_pp[, -1] <- lapply(envfit_table_pp[, -1], function(x) round(x, 3))

envfit_table_pp

write.csv(envfit_table_pp, "prok_envfit_results.csv", row.names = FALSE)

###### PLOTTING ####
# Extract r2 values from envfit
fit_df <- data.frame(
  element = rownames(fit$vectors$arrows),
  r2 = fit$vectors$r
)

# identify only variables with p < 0.05
sig_vars <- names(fit$vectors$pvals)[fit$vectors$pvals < 0.05]

# identify top 10 variables by r2
top10_vars <- fit_df %>% 
  arrange(desc(r2)) %>% 
  slice(1:10) %>% 
  pull(element)

# only keep variables which are present in both 
vars_to_keep <- intersect(sig_vars, top10_vars)
vars_to_keep

# Extract vector scores + p-values
arrow_df <- as.data.frame(scores(fit, "vectors"))
arrow_df$element <- rownames(arrow_df)

colnames(arrow_df)[1:2] <- c("PCo1", "PCo2")

# Filter to strongest variables
arrow_df <- arrow_df %>% 
  filter(element %in% vars_to_keep)

# Determine a safe scaling factor
max_pcoa <- max(abs(c(p3$res_ordination$scores[,1],
                      p3$res_ordination$scores[,2])))

max_arrow <- max(abs(c(arrow_df$PCo1, arrow_df$PCo2)))

arrow_scaler <- (0.2 * max_pcoa) / max_arrow   # arrows limited to 20% of axis range

p_pcoa_final_with_arrows <-
  p_pcoa_final +
  geom_segment(
    data = arrow_df,
    inherit.aes = FALSE,
    aes(x = 0, y = 0,
        xend = PCo1 * arrow_scaler,
        yend = PCo2 * arrow_scaler),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "red",
    linewidth = 0.8
  ) +
  geom_text_repel(
    data = arrow_df,
    inherit.aes = FALSE,
    aes(x = PCo1 * arrow_scaler,
        y = PCo2 * arrow_scaler,
        label = element),
    color = "red",
    size = 3,
    max.overlaps = Inf
  )

p_pcoa_final_with_arrows

