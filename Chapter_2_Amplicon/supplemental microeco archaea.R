# Load libraries
library(tidyverse)

# -----------------------------
# 1. Import data
# -----------------------------
# Phylum-level raw counts table: first column = "Taxa", rest = samples
phylum <- read.csv("phylum.csv", check.names = FALSE)

# Metadata table: must contain columns "sampleid", "sample_type", "sample_site"
metadata <- read.csv("metadata.csv")

# -----------------------------
# 2. Convert raw counts to relative abundances (%)
# -----------------------------
# Separate taxonomy and abundance data
taxa <- phylum$Taxa
abund <- phylum[, -1]

# Normalize each sample column to sum to 100
rel_abund <- sweep(abund, 2, colSums(abund), "/") * 100

# Reattach taxonomy labels
phylum_rel <- data.frame(Taxa = taxa, rel_abund)

# -----------------------------
# 3. Keep only Archaea
# -----------------------------
archaea_rel <- phylum_rel %>%
  filter(str_detect(Taxa, "k__Archaea"))

# -----------------------------
# 4. Reshape to long format
# -----------------------------
archaea_long <- archaea_rel %>%
  pivot_longer(
    cols = -Taxa,
    names_to = "sampleid",
    values_to = "rel_abundance"
  )

# -----------------------------
# 5. Join with metadata
# -----------------------------
archaea_long <- archaea_long %>%
  left_join(metadata, by = "sampleid")

# -----------------------------
# 6. Extract phylum names
# -----------------------------
archaea_long <- archaea_long %>%
  mutate(Phylum = str_extract(Taxa, "p__[^|]+") %>% str_replace("p__", ""))

# -----------------------------
# 7. Summarise mean relative abundance by sample_type
# -----------------------------
archaea_by_type <- archaea_long %>%
  group_by(Phylum, sample_type) %>%
  summarise(
    mean_rel_abundance = mean(rel_abundance, na.rm = TRUE),
    detected_samples = sum(rel_abundance > 0),
    .groups = "drop"
  ) %>%
  mutate(mean_rel_abundance = round(mean_rel_abundance, 3))

# -----------------------------
# 8. Pivot wider for easy comparison
# -----------------------------
archaea_compare <- archaea_by_type %>%
  select(Phylum, sample_type, mean_rel_abundance) %>%
  pivot_wider(
    names_from = sample_type,
    values_from = mean_rel_abundance
  )

# -----------------------------
# 9. Save output
# -----------------------------
write.csv(archaea_compare, "archaea_phyla_by_sampletype.csv", row.names = FALSE)

# -----------------------------
# 10. (Optional) Check total archaeal contribution to all prokaryotes
# -----------------------------
total_archaea <- archaea_long %>%
  group_by(sampleid, sample_type) %>%
  summarise(total_archaea = sum(rel_abundance), .groups = "drop") %>%
  group_by(sample_type) %>%
  summarise(mean_total_archaea = mean(total_archaea))
print(total_archaea)
