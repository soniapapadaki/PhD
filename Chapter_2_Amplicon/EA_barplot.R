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
# NOTE: you already put LOQ into value_raw, so use that.
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



### SAMPLING SITES AVERAGED
library(tidyverse)
library(ggpubr)

# 1. Filter only Colonised and Uncolonised Gypsum
df <- chem_wide %>%
  filter(sample_type %in% c("Colonised Gypsum", "Uncolonised Gypsum"))

# 2. Pivot element concentrations long
df_long <- df %>%
  select(sampleid, sample_type, sample_site, where(is.numeric)) %>%
  pivot_longer(cols = where(is.numeric),
               names_to = "element",
               values_to = "value") %>%
  # Ensure consistent order for plotting
  mutate(sample_type = factor(sample_type, levels = c("Uncolonised Gypsum", "Colonised Gypsum")))

# 3. Paired t-tests per element
ttest_results <- df_long %>%
  group_by(element) %>%
  summarise(
    p_value = t.test(
      x = value[sample_type == "Colonised Gypsum"],
      y = value[sample_type == "Uncolonised Gypsum"],
      paired = TRUE
    )$p.value
  ) %>%
  mutate(
    sig = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# 4. Merge significance into long dataframe
df_plot <- df_long %>%
  left_join(ttest_results, by = "element")

# 5. Plot: Initial plots with SE and asterisks for all elements
ggplot(df_plot, aes(x = sample_type, y = value, fill = sample_type)) +
  
  # Bars = mean
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  
  # Error bars = standard error
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  
  # Significance asterisks
  geom_text(
    aes(label = sig),
    stat = "summary",
    fun = median,
    vjust = -0.5,
    size = 5
  ) +
  
  # Facet by element
  facet_wrap(~ element, scales = "free_y") +
  
  # Labels and theme
  labs(title = "Elemental Concentration: Colonised vs Uncolonised Gypsum",
       x = NULL,
       y = "Concentration (mg/kg)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

### TO VIEW MEAN AND SE VALUES ###

# Compute mean, SE, and p-value for each element
df_summary_all <- df_plot %>%
  group_by(element) %>%
  summarise(
    mean_colonised   = mean(value[sample_type == "Colonised Gypsum"]),
    se_colonised     = sd(value[sample_type == "Colonised Gypsum"]) / sqrt(sum(sample_type == "Colonised Gypsum")),
    mean_uncolonised = mean(value[sample_type == "Uncolonised Gypsum"]),
    se_uncolonised   = sd(value[sample_type == "Uncolonised Gypsum"]) / sqrt(sum(sample_type == "Uncolonised Gypsum")),
    p_value          = t.test(
      value[sample_type == "Colonised Gypsum"],
      value[sample_type == "Uncolonised Gypsum"],
      paired = TRUE
    )$p.value,
    .groups = "drop"
  ) %>%
  # Add significance column based on p-value
  mutate(
    sig = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )

# View the summary table
df_summary_all


# 6. Plot: Clean barplots with SE and asterisks for selected elements

#BATCH 1
elements_to_plot <- c("Al", "As", "Ba", "C", "Ca", "Co", "Cr", "Cu", "Fe")
df_plot_sub <- df_plot %>%
  filter(element %in% elements_to_plot)

ggplot(df_plot_sub, aes(x = sample_type, y = value, fill = sample_type)) +
  
  # Bars = mean
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  
  # Error bars = standard error
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  
  # Significance asterisks
  geom_text(aes(label = sig),
            stat = "summary",
            fun = median,
            vjust = -0.5,
            size = 5) +
  
  # Facet by element
  facet_wrap(~ element, scales = "free_y") +
  
  # Labels and theme
  labs(title = "Colonised vs Uncolonised Gypsum",
       x = NULL,
       y = "Concentration (mg/kg)",
       fill = "Sample type") +   
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

#BATCH 2
elements_to_plot_b2 <- c("H", "K", "Mg", "Mn", "Na", "Ni", "P", "S", "Si")
df_plot_sub_b2 <- df_plot %>%
  filter(element %in% elements_to_plot_b2)

ggplot(df_plot_sub_b2, aes(x = sample_type, y = value, fill = sample_type)) +
  
  # Bars = mean
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  
  # Error bars = standard error
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  
  # Significance asterisks
  geom_text(aes(label = sig),
            stat = "summary",
            fun = median,
            vjust = -0.2,
            size = 5) +
  
  # Facet by element
  facet_wrap(~ element, scales = "free_y") +
  
  # Labels and theme
  labs(title = "Colonised vs Uncolonised Gypsum",
       x = NULL,
       y = "Concentration (mg/kg)",
       fill = "Sample type") +   
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

#BATCH 3
elements_to_plot_b3 <- c("Sr", "Ti", "V", "Zn")
df_plot_sub_b3 <- df_plot %>%
  filter(element %in% elements_to_plot_b3)

ggplot(df_plot_sub_b3, aes(x = sample_type, y = value, fill = sample_type)) +
  
  # Bars = mean
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  
  # Error bars = standard error
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  
  # Significance asterisks
  geom_text(aes(label = sig),
            stat = "summary",
            fun = median,
            vjust = -0.2,
            size = 5) +
  
  # Facet by element
  facet_wrap(~ element, scales = "free_y") +
  
  # Labels and theme
  labs(title = "Colonised vs Uncolonised Gypsum",
       x = NULL,
       y = "Concentration (mg/kg)",
       fill = "Sample type") +   
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )
