### LOAD PACKAGES
library(tidyverse)

### IMPORT DATA

# Phylum raw counts table (first column = "Taxa", rest of the columns = samples)
phylum <- read.csv("phylum.csv", check.names = FALSE)

# Metadata table (needs columns "sampleid", "sample_type", "sample_site")
metadata <- read.csv("metadata.csv")

### CONVERT RAW COUNTS TO RELATIVE ABUNDANCE

# Separate taxonomy and abundance data
taxa <- phylum$Taxa
abund <- phylum[, -1]

# Normalize each sample column to sum to 100
rel_abund <- sweep(abund, 2, colSums(abund), "/") * 100

# Reattach taxonomy labels
phylum_rel <- data.frame(Taxa = taxa, rel_abund)

### FILTER FOR ONLY ARCHAEA

archaea_rel <- phylum_rel %>%
  filter(str_detect(Taxa, "k__Archaea"))

### PUT IN LONG FORMAT

archaea_long <- archaea_rel %>%
  pivot_longer(
    cols = -Taxa,
    names_to = "sampleid",
    values_to = "rel_abundance"
  )

### ADD METADATA

archaea_long <- archaea_long %>%
  left_join(metadata, by = "sampleid")

### EXTRACT PHYLUM NAMES

archaea_long <- archaea_long %>%
  mutate(Phylum = str_extract(Taxa, "p__[^|]+") %>% str_replace("p__", ""))

### SUMMARY OF MEAN RELATIVE ABUNDANCE BY SAMPLE TYPE

archaea_by_type <- archaea_long %>%
  group_by(Phylum, sample_type) %>%
  summarise(
    mean_rel_abundance = mean(rel_abundance, na.rm = TRUE),
    detected_samples = sum(rel_abundance > 0),
    .groups = "drop"
  ) %>%
  mutate(mean_rel_abundance = round(mean_rel_abundance, 3))

archaea_compare <- archaea_by_type %>%
  select(Phylum, sample_type, mean_rel_abundance) %>%
  pivot_wider(
    names_from = sample_type,
    values_from = mean_rel_abundance
  )

### SAVE OUTPUT AS CSV

write.csv(archaea_compare, "archaea_phyla_by_sampletype.csv", row.names = FALSE)
