###############################################
# SULFUR CYCLING GENE ABUNDANCES 
#
# Description:
#   Generates bar plots showing mean abundances of sulfur-cycling 
#   genes (SCycDB) in gypsum endoliths, soil crusts, and beach 
#   endoliths. Abundances are summed by metabolic pathway and 
#   visualized with mean ± SD.
#
# Notes:
#	Requires gene abundance table produced via DIAMOND BLASTx 
#   functional annotation. Sulfur metabolism gene sequences 
#	used for annotation are from SCycDB (Yu et al., 2021):
#   https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13306
###############################################

### LOAD PACKAGES

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

### IMPORT AND CLEAN UP DATA

# read CSV
df <- read.csv("scyc_gene_abundances.csv", stringsAsFactors = FALSE)

# make sure sample columns are numeric
sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Metabolic.pathway"))]

df[sample_cols] <- lapply(df[sample_cols], function(x) as.numeric(as.character(x)))

# replace NAs with 0 
df[sample_cols] <- lapply(df[sample_cols], function(x) ifelse(is.na(x), 0, x))

# reshape into long format
df_long <- df %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance"
  )

# Trim whitespace in pathway names
df_long <- df_long %>%
  mutate(metabolic_pathway = trimws(Metabolic.pathway))
  
### SUM ABUNDANCES PER PATHWAY
  
df_pathway <- df_long %>%
  group_by(metabolic_pathway, sample) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

### ADD SAMPLE INFO

sample_type_lookup <- data.frame(
  sample = sample_cols,
  sample_type = rep(c("Arid endoliths", "Soil crusts", "Beach endoliths"), each = 3)
)

df_pathway <- df_pathway %>%
  left_join(sample_type_lookup, by = "sample")

### CALCULATE MEAN AND SD PER PATHWAY PER SAMPLE TYPE

df_pathway_summary <- df_pathway %>%
  group_by(metabolic_pathway, sample_type) %>%
  summarise(
    mean_abundance = mean(total_abundance),
    sd_abundance = sd(total_abundance),
    .groups = "drop"
  )

### FILTER OUT UNWANTED METABOLIC PATHWAYS

df_plot <- df_pathway_summary %>%
  filter(!metabolic_pathway %in% "Assimilatory sulfate reduction/Dissimilatory sulfur reduction/oxidation")

### PLOT

df_plot <- df_plot %>%
  mutate(metabolic_pathway = str_wrap(metabolic_pathway, width = 30))

ggplot(df_plot, 
       aes(x = sample_type, y = mean_abundance, fill = sample_type)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.7), 
           width = 0.5) +  # thinner bars
  geom_errorbar(aes(ymin = mean_abundance - sd_abundance, 
                    ymax = mean_abundance + sd_abundance),
                width = 0.2, 
                position = position_dodge(width = 0.7)) +
  facet_wrap(~ metabolic_pathway, scales = "free_y") +
  scale_fill_manual(values = c(
    "Arid endoliths" = "#1b9e77",
    "Beach endoliths" = "#d95f02",
    "Soil crusts" = "#7570b3"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  
  labs(
    x = NULL,  
    y = "Average gene copies per organism",
    fill = "Sample type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "gray90"), 
    strip.text = element_text(face = "bold", size = 8),
    panel.grid = element_blank(),
    axis.text.y = element_text(color = "black")
  )