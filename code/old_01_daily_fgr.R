### Fast fires impacts on biodiversity

### This script plots daily fire perimeters of the five largest fires
### Matt Bitters
### matthew.bitters@colorado.edu


# Install and load required packages
packages <- c("here", "sf", "dplyr", "ggplot2", "ggspatial", "RColorBrewer")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(RColorBrewer)



### Setup

# Path to geodatabase
fired_daily_path <- here("data", "raw", "fired", "fired_australia_2000_to_2024_daily.gpkg")

# Check layers
st_layers(fired_daily_path)

# Read combined fire dataset
firedailydat <- st_read(dsn = fired_daily_path, layer = "fired_australia_2000_to_2024_daily")

# Check crs
crs(firedailydat)

# Check fields
names(firedailydat)

# Convert km² to hectares (1 km² = 100 ha)
firedailydat$mx_grw_ha <- firedailydat$mx_grw_km2 * 100

# Filter to just top five fastest fires
firedailydat <- firedailydat %>%
  filter(id %in% c("779368", "72475", "57181", "780632", "780645"))

# Check nrow
nrow(firedailydat) #191 rows coming from five largest fires





### Single fire test




# Filter to one fire by ID 
single_fire <- firedailydat %>%
  filter(id == "780645") %>%
  arrange(event_day)  

# Assign a date-based factor for coloring
single_fire <- single_fire %>%
  mutate(date_factor = as.factor(event_day))  

# Create sf point from ignition coordinates
sinusoidal_crs <- "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"

ignition_point <- single_fire %>%
  slice(1) %>%  # One row per fire is enough
  st_as_sf(coords = c("ig_utm_x", "ig_utm_y"), crs = sinusoidal_crs)

# Get ignition coordinates as a data frame
ignition_df <- ignition_point %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(x = X, y = Y)

# Lookup vector
fire_names_vec <- c(
  "779368" = "Dunns Road Fire (779368)",
  "72475" = "Green Valley – Talmalmo (72475)",
  "57181" = "Kangaroo Island Complex (57181)",
  "780632" = "Shark Creek / Myall Creek (780632)",
  "780645" = "Yarrowitch Complex (780645)"
)

# Get current fire id as string
current_id <- as.character(unique(single_fire$id))

# Get label for current fire
subtitle_text <- fire_names_vec[current_id]

# Plot
ggplot(single_fire) +
  geom_sf(aes(fill = as.integer(date_factor)), color = NA, alpha = 0.7) +
  geom_point(data = ignition_df, aes(x = x, y = y), 
             shape = 4, color = "black", size = 4) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(9, "Reds")),
    name = "Days since ignition",
    breaks = seq(1, length(levels(single_fire$date_factor)), by = 4),
    labels = seq(1, length(levels(single_fire$date_factor)), by = 4),
    guide = guide_colorbar(reverse = TRUE)
  ) +
  labs(
    title = "Daily fire perimeter growth",
    subtitle = subtitle_text,
    x = "", y = ""
  ) +
  annotation_scale(location = "br", width_hint = 0.3) +        # scale bar bottom left
  annotation_north_arrow(location = "br", which_north = "true", # compass rose bottom left
                         pad_x = unit(0.1, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.subtitle = element_text(hjust = 0)
  )









### Loop through all five fastest fires




# Vector of fire IDs and names
fire_names_vec <- c(
  "779368" = "Dunns Road Fire (779368)",
  "72475" = "Green Valley – Talmalmo (72475)",
  "57181" = "Kangaroo Island Complex (57181)",
  "780632" = "Shark Creek / Myall Creek (780632)",
  "780645" = "Yarrowitch Complex (780645)"
)

# Make sure the output directory exists
if(!dir.exists("figures")) dir.create("figures")

sinusoidal_crs <- "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"

for (fire_id in names(fire_names_vec)) {
  
  # Filter to one fire
  single_fire <- firedailydat %>%
    filter(id == fire_id) %>%
    arrange(event_day) %>%
    mutate(date_factor = as.factor(event_day))
  
  # Get centroid of the first day's fire perimeter
  ignition_point <- single_fire %>%
    filter(event_day == min(event_day)) %>%
    st_union() %>%
    st_centroid()
  
  # Convert to data frame for plotting
  ignition_df <- st_coordinates(ignition_point) %>%
    as.data.frame() %>%
    rename(x = X, y = Y)
  
  print(ignition_df)
  
  subtitle_text <- fire_names_vec[fire_id]
  
  # Determine max fire duration for the current fire
  max_day <- max(as.integer(levels(single_fire$date_factor)))
  
  # Create nicely spaced breaks, always starting with 1
  breaks_vec <- unique(c(1, pretty(2:max_day, n = 6)))
  
  p <- ggplot(single_fire) +
    geom_sf(aes(fill = as.integer(date_factor)), color = NA, alpha = 0.7) +
    geom_point(data = ignition_df, aes(x = x, y = y),
               shape = 4, color = "black", size = 2, stroke = 2) +
    # Use this in the fill scale
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(9, "Reds")),
      name = "Days since ignition",
      breaks = breaks_vec,
      labels = breaks_vec,
      guide = guide_colorbar(reverse = TRUE)
    ) +
    labs(
      title = "Daily fire perimeter growth",
      subtitle = subtitle_text,
      x = "", y = ""
    ) +
    annotation_scale(location = "br", width_hint = 0.2) +
    annotation_north_arrow(
      location = "br",
      which_north = "true",
      pad_x = unit(0.1, "in"),
      pad_y = unit(0.3, "in"),
      style = north_arrow_fancy_orienteering(
        line_width = 0.5,
        text_size = 6
      ),
      height = unit(0.4, "in"),  
      width = unit(0.4, "in")    
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      plot.subtitle = element_text(hjust = 0),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Save plot
  ggsave(
    filename = paste0("figures/fire_", fire_id, ".png"),
    plot = p,
    width = 6,
    height = 6,
    dpi = 300
  )
}














