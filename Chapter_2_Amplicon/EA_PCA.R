library(tidyverse)
library(readr)
library(compositions)
library(ggrepel)

# ─────────────────────────────────────────
# 1. Import data
# ─────────────────────────────────────────
chem     <- read.delim("EA_concentrations.csv", stringsAsFactors = FALSE)
loq      <- read.delim("EA_LOQ.csv", stringsAsFactors = FALSE)
metadata <- read.delim("EA_metadata.csv", stringsAsFactors = FALSE)

# ─────────────────────────────────────────
# 2. Convert all values to mg/kg
# ─────────────────────────────────────────
chem <- chem %>%
  mutate(value_mgkg = case_when(
    unit == "wt%"  ~ value_raw * 10000,   # convert wt% → mg/kg
    unit == "mg/kg" ~ value_raw,          # already mg/kg
    TRUE ~ NA_real_
  ))

# ─────────────────────────────────────────
# 3. Replace censored values (BDL) with LOQ/sqrt(2)
# ─────────────────────────────────────────
chem <- chem %>%
  mutate(value_clean = case_when(
    BDL == "<" ~ value_mgkg / sqrt(2),
    TRUE       ~ value_mgkg
  ))

# ─────────────────────────────────────────
# 4. Pivot into wide form
# ─────────────────────────────────────────
chem_wide <- chem %>%
  select(sampleid, element, value_clean) %>%
  pivot_wider(names_from = element,
              values_from = value_clean)

# ─────────────────────────────────────────
# 5. Add metadata
# ─────────────────────────────────────────
chem_wide <- chem_wide %>%
  left_join(metadata, by = "sampleid")

# Check for missing sample metadata:
missing_meta <- chem_wide %>% filter(is.na(sample_site) | is.na(sample_type) | is.na(method))
missing_meta   # should be empty

# ─────────────────────────────────────────
# 6. Build compositional numeric matrix
# (INCLUDE ONLY ELEMENT COLUMNS)
# ─────────────────────────────────────────
elem_matrix <- chem_wide %>%
  select(-sampleid, -sample_site, -sample_type, -method) %>%
  as.matrix()

# Make sure rows stay aligned with sampleid
rownames(elem_matrix) <- chem_wide$sampleid

# ─────────────────────────────────────────
# 7. CLR transformation
# ─────────────────────────────────────────
elem_clr <- clr(elem_matrix)

# ─────────────────────────────────────────
# 8. PCA on CLR-transformed data
# ─────────────────────────────────────────
pca <- prcomp(elem_clr, scale = FALSE)

# ─────────────────────────────────────────
# 9. Build PCA dataframe with metadata
# ─────────────────────────────────────────
pca_df <- as.data.frame(pca$x) %>%
  mutate(sampleid   = rownames(elem_clr)) %>%
  left_join(metadata, by = "sampleid")

# ─────────────────────────────────────────
# 10. PCA PLOT 
# ─────────────────────────────────────────
ggplot(pca_df, aes(PC1, PC2,
                   color = sample_site,
                   shape = sample_type)) +
  
  geom_point(size = 3, alpha = 0.9) +
  
  # Sample labels (excluded from legend)
  geom_text_repel(aes(label = sampleid),
                  size = 3,
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  
  # Filled, print-friendly shapes
  scale_shape_manual(values = c(16, 17, 15)) +
  
  labs(
    title = "PCA of Elemental Composition",
    color = "Sampling site",
    shape = "Sample type",
    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )



#######################################

# ─────────────────────────────────────────
# FILTER: Only colonised + uncolonised gypsum
# ─────────────────────────────────────────
gypsum_only <- chem_wide %>%
  filter(sample_type %in% c("Colonised Gypsum", "Uncolonised Gypsum"))

# ─────────────────────────────────────────
# Build element-only matrix for gypsum
# ─────────────────────────────────────────
elem_matrix_gypsum <- gypsum_only %>%
  select(-sampleid, -sample_site, -sample_type, -method) %>%
  as.matrix()

rownames(elem_matrix_gypsum) <- gypsum_only$sampleid

# ─────────────────────────────────────────
# CLR transform (compositional normalisation)
# ─────────────────────────────────────────
elem_clr_gypsum <- clr(elem_matrix_gypsum)

# ─────────────────────────────────────────
# PCA on gypsum-only CLR values
# ─────────────────────────────────────────
pca_gypsum <- prcomp(elem_clr_gypsum, scale = FALSE)

# ─────────────────────────────────────────
# Build PCA dataframe with metadata
# ─────────────────────────────────────────
pca_gypsum_df <- as.data.frame(pca_gypsum$x) %>%
  mutate(sampleid = rownames(elem_clr_gypsum)) %>%
  left_join(metadata, by = "sampleid")

# ─────────────────────────────────────────
# Gypsum-only PCA plot
# ─────────────────────────────────────────
ggplot(pca_gypsum_df, aes(PC1, PC2,
                          color = sample_site,
                          shape = sample_type)) +
  geom_point(size = 3, alpha = 0.9) +
  
  # Sample labels (not shown in legend)
  geom_text_repel(aes(label = sampleid),
                  size = 3,
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  
  # Manual shapes (print-safe)
  scale_shape_manual(values = c(16, 17, 15)) +
  
  labs(
    title = "PCA of Gypsum Samples Only",
    color = "Sampling site",
    shape = "Sample type",
    x = paste0("PC1 (", round(summary(pca_gypsum)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_gypsum)$importance[2,2] * 100, 1), "%)")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )


