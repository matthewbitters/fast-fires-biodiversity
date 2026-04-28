### Fast fires impacts on biodiversity

### This script creates necessary folders and plots fire perimeters
### Matt Bitters
### matthew.bitters@colorado.edu


# Install and load required packages
packages <- c("here", "sf", "dplyr", "ggplot2", "terra", "rnaturalearth", "rnaturalearthdata", "ggspatial", "scales", "ggrepel")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(scales)
library(ggrepel)


# Set up necessary directories
dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Create subfolders inside 'data'
dir.create("data/raw", showWarnings = FALSE)
dir.create("data/raw/fired", showWarnings = FALSE)
dir.create("data/derived", showWarnings = FALSE)
dir.create("data/output", showWarnings = FALSE)



# ==== READ AND FILTER DATA ====

# Path to geodatabase
fired_path <- here("data", "raw", "fired", "fired_australia_2000_to_2024_events.gpkg")

# Check layers
st_layers(fired_path)

# Read combined fire dataset
firedat <- st_read(dsn = fired_path, layer = "fired_australia_2000_to_2024_events")

# Check crs
crs(firedat)

# Check fields
names(firedat)

# Check and convert date
str(firedat$ig_date)
firedat$ig_date <- as.Date(firedat$ig_date)

# Filter between July 2019 and May 2020
firedat <- firedat %>%
  filter(ig_date >= as.Date("2019-06-01") & ig_date <= as.Date("2020-05-31"))

# Double check dates and number of rows
summary(firedat$ig_date)
nrow(firedat) # 26,126 total files








### Make max fgr histogram (all in ha/day)

# Single fastest fire
max(firedat$mx_grw_km2*100) # 261,454.3 ha/day

# Raw
ggplot(firedat, aes(x = mx_grw_km2 * 100)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  theme_classic() +
  labs(x = "max FRG (ha/day)", y = "Count")

# Logged
ggplot(firedat, aes(x = mx_grw_km2 * 100)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  scale_x_log10() +
  theme_classic() +
  labs(x = "log max FRG (ha/day)", y = "Count")

# Faceted by land cover type
ggplot(firedat, aes(x = mx_grw_km2 * 100)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  scale_x_log10() +
  facet_wrap(~lc_name, nrow = 3) +
  theme_minimal() +
  labs(x = "log max FRG (ha/day)", y = "Count") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )






### Plot FGR on map (each point represents a fire; color represents fgr)

# Convert km² to hectares (1 km² = 100 ha)
firedat$mx_grw_ha <- firedat$mx_grw_km2 * 100

# Remove inherited polygon geometry if present
firedat <- st_drop_geometry(firedat)

# Assign custom Sinusoidal projection
sinusoidal_crs <- "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"

# Create points from coordinate columns directly (overwrite previous sf object)
fire_points <- st_as_sf(
  firedat,
  coords = c("ig_utm_x", "ig_utm_y"),
  crs = sinusoidal_crs,
  remove = FALSE
)

fire_points_wgs84 <- st_transform(fire_points, crs = 4326)

# Get Australia map
australia <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(admin == "Australia")

ggplot() +
  geom_sf(data = australia, fill = "gray95", color = "black") +
  geom_sf(data = fire_points_wgs84, aes(color = mx_grw_km2 * 100), size = 1, alpha = 0.7) +
  scale_color_viridis_c(
    name = "Max Growth Rate (ha)",
    trans = "sqrt",
    labels = comma
  ) +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10), expand = FALSE) +
  labs(
    title = "Australian Bushfires (2019–2020)",
    subtitle = "Each point = one fire, colored by max growth rate (ha)",
    x = "", y = ""
  ) +
  theme_minimal() +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true",
                         style = north_arrow_fancy_orienteering)





### Add in fire size to scale points on map

# Read in data
fire_perims <- read_sf(here("data", "raw", "fired", "fired_australia_2000_to_2024_events.shp"))

# Ensure fire_perims is in a projected CRS (meters) for area calculation
fire_perims_proj <- st_transform(fire_perims, crs = "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")

# Calculate area in hectares
fire_perims_proj$fire_size_ha <- as.numeric(st_area(fire_perims_proj)) / 10000  # m² → ha

# Summarize total fire size per id
fire_area_summary <- fire_perims_proj |>
  st_drop_geometry() |>
  group_by(id) |>
  summarise(fire_size_ha = sum(fire_size_ha, na.rm = TRUE))

# Join to fire_points by id
fire_points <- fire_points |>
  left_join(fire_area_summary, by = "id")

# Reproject fire points to WGS84 for plotting
fire_points_wgs84 <- st_transform(fire_points, crs = 4326)

# Load Australia base map
australia <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  dplyr::filter(admin == "Australia")

# Plot: Each point = one fire
ggplot() +
  geom_sf(data = australia, fill = "gray95", color = "black") +
  geom_sf(data = fire_points_wgs84,
          aes(color = mx_grw_ha, size = fire_size_ha),
          alpha = 0.7) +
  scale_color_viridis_c(
    name = "Max Growth Rate (ha)",
    trans = "sqrt",
    labels = comma
  ) +
  scale_size_continuous(
    name = "Total Fire Size (ha)",
    range = c(0.5, 6),         # adjust for visual balance
    trans = "sqrt",            # useful for skewed size data
    labels = comma
  ) +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10), expand = FALSE) +
  labs(
    title = "Australian Bushfires (2000–2024)",
    subtitle = "Each point = one fire\nColor = max daily growth (ha), Size = total area burned (ha)",
    x = "", y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         style = north_arrow_fancy_orienteering)




### Only show top 100 fastest fires

# Add a flag for the top 100 fires by max growth rate
fire_points_wgs84 <- fire_points_wgs84 |>
  arrange(desc(mx_grw_ha)) |>
  mutate(top_100 = row_number() <= 100)

# Arrange points
fire_points_wgs84 <- fire_points_wgs84 %>%
  arrange(mx_grw_ha)  # small to large, so large ones are last and drawn on top

# Plot
ggplot() +
  geom_sf(data = australia, fill = "gray95", color = "black") +
  
  # Background: all other fires
  geom_sf(data = filter(fire_points_wgs84, !top_100),
          color = "grey70", size = 0.4, alpha = 0.5) +
  
  # Highlighted: top 100 fastest-growing fires
  geom_sf(data = filter(fire_points_wgs84, top_100),
          aes(color = mx_grw_ha, size = fire_size_ha),
          alpha = 0.8) +
  scale_color_gradientn(
    name = "Max fire growth rate (ha/day)",
    colours = c("#FEE08B", "#F46D43", "#B2182B", "#440154"), 
    trans = "sqrt",
    labels = comma
  ) +
  scale_size_continuous(
    name = "Fire size (ha)",
    range = c(1, 4),
    trans = "sqrt",
    labels = comma
  ) +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10), expand = FALSE) +
  labs(
    title = "Top 100 fastest growing fires during the 2019-2020 Black Summer bushfires",
    x = "", y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering)






### Identify and label fastest fires
top5_fires <- fire_points_wgs84 |>
  arrange(desc(mx_grw_ha)) |>
  slice_head(n = 5)

# Label top 5
top5_fires$fire_name <- c(
  "Dunns Road Fire",
  "Green Valley – Talmalmo",
  "Kangaroo Island Complex",
  "Shark Creek / Myall Creek",
  "Yarrowitch Complex"
)


# Plot
fast_fires_aus <- ggplot() +
  geom_sf(data = australia, fill = "gray95", color = "black") +
  
  # Background: all other fires
  geom_sf(data = filter(fire_points_wgs84, !top_100),
          color = "grey70", size = 0.4, alpha = 0.5) +
  
  # Highlighted: top 100 fastest-growing fires
  geom_sf(data = filter(fire_points_wgs84, top_100),
          aes(color = mx_grw_ha, size = fire_size_ha),
          alpha = 0.8) +
  scale_color_gradientn(
    name = "Max fire growth rate (ha/day)",
    colours = c("#FEE08B", "#F46D43", "#B2182B", "#440154"), 
    trans = "sqrt",
    labels = comma,
    guide = guide_colorbar(order = 1)
  ) +
  scale_size_continuous(
    name = "Fire size (ha)",
    range = c(1, 4),
    trans = "sqrt",
    labels = comma,
    guide = guide_legend(order = 2)
  ) +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10), expand = FALSE) +
  labs(
    title = "Top 100 fastest growing fires during the 2019-2020 Black Summer bushfires",
    x = "", y = ""
  ) +
  geom_text_repel(
    data = top5_fires,
    aes(geometry = geometry, label = fire_name),
    stat = "sf_coordinates",
    size = 3,                  
    fontface = "bold",        
    segment.size = 0.6,       
    segment.color = "black",  
    box.padding = 0.6,
    point.padding = 0.4,
    min.segment.length = 0
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(legend.position = "right") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering)



ggsave(
  filename = "figures/fast_fires_aus.png",  # Change this to your desired path and filename
  plot = fast_fires_aus,          # The ggplot object
  width = 7,              # Width in inches
  height = 5,              # Height in inches
  dpi = 300                # Resolution (dots per inch), 300 is good for print quality
)



###-----------------------------------------------------------------------------


### Max FGR vs. final fire size


# Convert total area burned to hectares
firedat <- firedat %>%
  mutate(fire_size_ha = tot_ar_km2 * 100) %>%
  filter(mx_grw_ha > 0, fire_size_ha > 0)

# Fit log-log linear model
lm_log <- lm(log(fire_size_ha) ~ log(mx_grw_ha), data = firedat)

# Extract coefficients
coef_intercept <- coef(lm_log)[1]
coef_slope <- coef(lm_log)[2]
adj_r2 <- summary(lm_log)$adj.r.squared

# Format equation text (rounded)
eq_text <- paste0(
  "log(Final size) = ", round(coef_intercept, 2),
  " + ", round(coef_slope, 2), " * log(FGR)\n",
  "Adj. R² = ", round(adj_r2, 3)
)

# Create plot with annotation
fgr_size <- ggplot(firedat, aes(x = mx_grw_ha, y = fire_size_ha)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = TRUE) +
  scale_x_log10(labels = comma_format()) +
  scale_y_log10(labels = comma_format()) +
  labs(
    x = "Max fire growth rate (ha/day, log scale)",
    y = "Final fire size (ha, log scale)",
    ) +
  annotate(
    "text",
    x = min(firedat$mx_grw_ha, na.rm = TRUE) * 1.2,
    y = max(firedat$fire_size_ha, na.rm = TRUE) / 2,
    label = eq_text,
    color = "black",
    size = 4,
    hjust = 0
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = "figures/fgr_size.png",  # Change this to your desired path and filename
  plot = fgr_size,          # The ggplot object
  width = 5,              # Width in inches
  height = 5,              # Height in inches
  dpi = 300                # Resolution (dots per inch), 300 is good for print quality
)






###----------------------------------------------------------------------------




### Calculate mean FGR per ecosystem type



# Group and summarize
lc_summary <- firedat %>%
  group_by(lc_name) %>%
  summarise(
    n = n(),
    mean_grw = mean(mx_grw_ha, na.rm = TRUE),
    sd_grw = sd(mx_grw_ha, na.rm = TRUE),
    se_grw = sd_grw / sqrt(n),
    ci95 = qt(0.975, df = n - 1) * se_grw,
    lower_ci = mean_grw - ci95,
    upper_ci = mean_grw + ci95
  ) %>%
  arrange(desc(mean_grw))

# Plot
lc_means <- ggplot(lc_summary, aes(x = mean_grw, y = reorder(lc_name, mean_grw))) +
  geom_pointrange(
    aes(xmin = lower_ci, xmax = upper_ci),
    color = "black", size = 0.4
  ) +
  geom_text(
    aes(label = paste0("n = ", n)),
    hjust = -0.2, vjust = 1.2, size = 3.5
  ) +
  labs(
    x = "Mean max fire growth rate (ha/day)",
    y = "Land cover type"
    ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  xlim(-250, max(lc_summary$upper_ci) * 1.1)  



ggsave(
  filename = "figures/lc_means.png",  # Change this to your desired path and filename
  plot = lc_means,          # The ggplot object
  width = 7,              # Width in inches
  height = 5,              # Height in inches
  dpi = 300                # Resolution (dots per inch), 300 is good for print quality
)





###---------------------------------------------------------------------------



### Calculate overall mean, minimum, max, and threshold for FGR

# Remove missing values
grw_values <- na.omit(firedat$mx_grw_ha)

# Sample size
n <- length(grw_values)

# Mean and standard error
mean_grw <- mean(grw_values)
se_grw <- sd(grw_values) / sqrt(n)

# 95% CI
ci_95 <- qt(0.975, df = n - 1) * se_grw
lower_ci <- mean_grw - ci_95
upper_ci <- mean_grw + ci_95

# Print
cat("Mean max fire growth rate:", round(mean_grw, 1), "ha/day\n")
cat("95% Confidence Interval: [", round(lower_ci, 1), ",", round(upper_ci, 1), "]\n")

# Min FGR
cat("Minimum max fire growth rate:", round(min(firedat$mx_grw_ha), 1), "ha/day\n")

# Max FGR
cat("Maximum max fire growth rate:", round(max(firedat$mx_grw_ha), 1), "ha/day\n")

# 97th percentile
quantile_97 <- quantile(firedat$mx_grw_ha, probs = 0.97, na.rm = TRUE)
cat("97th percentile of max fire growth rate:", round(quantile_97, 1), "ha/day\n")

# Count how many values exceed 97th percentile
n_total <- sum(!is.na(firedat$mx_grw_ha))
n_above <- sum(firedat$mx_grw_ha > quantile_97, na.rm = TRUE)

# Number above
cat("Number of fires above the 97th percentile:", n_above, "\n")

# Calculate percentage
percent_above <- (n_above / n_total) * 100

cat("Percentage of fires above the 97th percentile:", round(percent_above, 2), "%\n")



# Plot overall histogram
overall_hist <- ggplot(firedat, aes(x = mx_grw_ha)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = quantile_97, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = quantile_97, y = 9000, 
           label = paste0("97th percentile\n(", round(quantile_97, 1), " ha/day)"),
           vjust = 1.3, hjust = -0.1, color = "red", size = 4) +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Max fire growth rate (ha/day, log scale)", y = "Count") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) 

ggsave(
  filename = "figures/overall_hist.png",  # Change this to your desired path and filename
  plot = overall_hist,          # The ggplot object
  width = 5,              # Width in inches
  height = 3,              # Height in inches
  dpi = 300                # Resolution (dots per inch), 300 is good for print quality
)


# Plot histogram faceted by land cover

# Compute 97th percentiles by land cover
p97_df <- firedat %>%
  group_by(lc_name) %>%
  summarise(p97 = quantile(mx_grw_ha, 0.97, na.rm = TRUE))

# Estimate y-position for labels (max bin count)
hist_max_y <- firedat %>%
  mutate(bin = cut(mx_grw_ha, breaks = 50)) %>%
  group_by(lc_name, bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(lc_name) %>%
  summarise(y = max(n), .groups = "drop")

# Join to get label data
label_df <- p97_df %>%
  left_join(hist_max_y, by = "lc_name") %>%
  mutate(label = paste0(round(p97, 1), " ha/day"))

# Main plot (no text in geom_text yet!)
firedat_p97 <- firedat %>%
  left_join(p97_df, by = "lc_name")

# Plot
lc_hist <- ggplot(firedat_p97, aes(x = mx_grw_ha)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(data = p97_df, aes(xintercept = p97), color = "red", linetype = "dashed", size = 0.8) +
  geom_text(
    data = label_df,
    aes(x = p97, y = y * 1.1, label = label),
    inherit.aes = FALSE,
    color = "red", hjust = -0.1, vjust = 1.2, size = 3
  ) +
  facet_wrap(~lc_name, nrow = 5, scales = "free_y") +
  scale_x_log10(labels = scales::comma) +
  labs(
    x = "Max fire growth rate (ha/day, log scale)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 8)
  )

# Save 
ggsave(
  filename = "figures/lc_hist.png",
  plot = lc_hist,
  width = 7,
  height = 6,
  dpi = 300
)
