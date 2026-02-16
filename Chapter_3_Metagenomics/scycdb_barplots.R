

setwd("~/Writing Up [WD]/Chapter 3 - Metagenomics/metagenomics figs/raw images/scycdb barplot")

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# read your CSV
df <- read.csv("scyc_gene_abundances.csv", stringsAsFactors = FALSE)

# Make sure sample columns are numeric
sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Metabolic.pathway"))]

df[sample_cols] <- lapply(df[sample_cols], function(x) as.numeric(as.character(x)))

# Replace NAs with 0 (if NA means no abundance)
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

###### sum abundances per pathway

df_pathway <- df_long %>%
  group_by(metabolic_pathway, sample) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")


#### add sample info ###

sample_type_lookup <- data.frame(
  sample = sample_cols,
  sample_type = rep(c("Arid endoliths", "Soil crusts", "Beach endoliths"), each = 3)
)

df_pathway <- df_pathway %>%
  left_join(sample_type_lookup, by = "sample")

####calculate mean and SD per pathway per sample type

df_pathway_summary <- df_pathway %>%
  group_by(metabolic_pathway, sample_type) %>%
  summarise(
    mean_abundance = mean(total_abundance),
    sd_abundance = sd(total_abundance),
    .groups = "drop"
  )


######plot

# Filter out unwanted metabolic pathways
df_plot <- df_pathway_summary %>%
  filter(!metabolic_pathway %in% "Assimilatory sulfate reduction/Dissimilatory sulfur reduction/oxidation")


df_plot <- df_plot %>%
  mutate(metabolic_pathway = str_wrap(metabolic_pathway, width = 30))

# Plot
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # 0% below, 5% above
  labs(
    x = NULL,  # remove x-axis labels
    y = "Average gene copies per organism",
    fill = "Sample type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "gray90"),  # optional: nicer facet header
    strip.text = element_text(face = "bold", size = 8),
    panel.grid = element_blank(),
    axis.text.y = element_text(color = "black")
  )