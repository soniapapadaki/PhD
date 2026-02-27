###############################################
# Svalbard map
#
# Description:
#   Creates general map of Svalbard with option to add markers for sampling sites
#
###############################################

## LOAD PACKAGES

library(ggplot2)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(tidyverse)

## GET BASE MAP OF SVALBARD

# Load subnational regions and filter for Svalbard
svalbard <- ne_states(country = "Norway", returnclass = "sf") %>%
  filter(name == "Svalbard")

# Ensure the coordinate reference system is correct
svalbard <- st_transform(svalbard, crs = 4326)  # Standard lat/lon


## ADD SAMPLING SITE COORDINATES

#sampling_site <- data.frame(
  #name = "My Sampling Site",
  #lon = 15.6,   # Longitude
  #lat = 78.2    # Latitude
#)

# Convert to sf object
#sampling_site_sf <- st_as_sf(sampling_site, coords = c("lon", "lat"), crs = 4326)


## CREATE MAP

ggplot() +
  geom_sf(data = svalbard, fill = "gray90", color = "black") +  # Svalbard Mainland
  #geom_sf(data = sampling_site_sf, color = "red", size = 3) +   # Sampling Site
  annotation_scale(location = "bl", width_hint = 0.3) + 
  coord_sf(xlim = c(9, 30), ylim = c(76, 80.5)) +  # Zoomed-in extent
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue", color = NA),
        panel.grid.major = element_line(color = "grey60", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey80", linetype = "dashed")
)