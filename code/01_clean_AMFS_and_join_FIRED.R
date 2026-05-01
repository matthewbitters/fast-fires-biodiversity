### Fast fires impacts on biodiversity

### This script creates necessary folders, filters to suitable projects, cleans 
### AMFS data for joining with FIRED FGRs, and joins FIRED FGRs.

### Matt Bitters
### matthew.bitters@colorado.edu




# =========================================================
# Install and load required packages
# =========================================================

packages <- c("here", "dplyr", "tidyr", "sf", "lubridate", "ggplot2", "stringr", "readr")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(dplyr)
library(tidyr)
library(sf)
library(lubridate)
library(ggplot2)
library(stringr)
library(readr)


# =========================================================
# Set up directories
# =========================================================
dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

dir.create("data/raw", showWarnings = FALSE)
dir.create("data/raw/fired", showWarnings = FALSE)
dir.create("data/derived", showWarnings = FALSE)
dir.create("data/output", showWarnings = FALSE)


# =========================================================
# Read and inspect AMFS data
# =========================================================

amfs <- read.csv("data/raw/The_Australian_MegaFire_Synthesis_Dataset.csv")

str(amfs)
head(amfs)
nrow(amfs) # 809,831 total observations


# =========================================================
# Summarize projects
# =========================================================

# Inspect project function
inspect_project <- function(df, pid) {
  
  library(dplyr)
  
  dat <- df %>%
    filter(unique_project_ID == pid)
  
  fix_inf <- function(x) {
    ifelse(is.infinite(x), NA_real_, x)
  }
  
  summary_tbl <- dat %>%
    summarise(
      n_rows = n(),
      n_sites = n_distinct(site_name),
      
      hours_min = fix_inf(suppressWarnings(min(survey_effort_hours, na.rm = TRUE))),
      hours_max = fix_inf(suppressWarnings(max(survey_effort_hours, na.rm = TRUE))),
      traps_min = fix_inf(suppressWarnings(min(survey_effort_traps, na.rm = TRUE))),
      traps_max = fix_inf(suppressWarnings(max(survey_effort_traps, na.rm = TRUE))),
      area_min  = fix_inf(suppressWarnings(min(survey_area, na.rm = TRUE))),
      area_max  = fix_inf(suppressWarnings(max(survey_area, na.rm = TRUE))),
      
      response_min = fix_inf(suppressWarnings(min(response_value, na.rm = TRUE))),
      response_max = fix_inf(suppressWarnings(max(response_value, na.rm = TRUE))),
      response_unique = n_distinct(response_value)
    )
  
  print(summary_tbl, width = 200)
  
  cat("\n--- Diagnostics ---\n")
  
  if (!is.na(summary_tbl$hours_min) && summary_tbl$hours_min == summary_tbl$hours_max) {
    cat("survey_effort_hours is constant\n")
  }
  if (!is.na(summary_tbl$traps_min) && summary_tbl$traps_min == summary_tbl$traps_max) {
    cat("survey_effort_traps is constant\n")
  }
  if (!is.na(summary_tbl$area_min) && summary_tbl$area_min == summary_tbl$area_max) {
    cat("survey_area is constant\n")
  }
  if (!is.na(summary_tbl$response_unique) && summary_tbl$response_unique <= 2) {
    cat("response_value looks binary or near-binary\n")
  }
  
  cat("\n--- Sample rows ---\n")
  
  dat %>%
    select(
      site_name,
      survey_method,
      survey_effort_hours,
      survey_effort_traps,
      survey_area,
      survey_effort_text,
      response_variable,
      response_variable_text,
      response_value,
      response_units,
      broad_taxon
    ) %>%
    slice_head(n = 20) %>%
    tibble::as_tibble() %>%
    print(n = 20, width = 200)
}



# Check one project at a time
# Then fill in summary table below
inspect_project(amfs, "UP_085")




# Make project summary table
project_review <- tribble(
  ~unique_project_ID, ~include_pa, ~include_abundance, ~effort_clean, ~effort_variable, ~effort_use_type, ~response_transformation, ~response_type, ~notes,
  
  "UP_003", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Vegetation mortality metrics (% basal area dead, % stems dead, % topkill); not abundance or presence.",
  "UP_004", TRUE, FALSE, "camera deployment days constant within project", "survey_effort_traps", "none", "convert response_value > 0 to presence/absence", "binary_derived_from_index", "Camera trap relative abundance index; not used as abundance.",
  "UP_005", TRUE, TRUE, "plot area varies", "survey_area", "offset_abundance_covariate_pa", "retain raw counts; convert count > 0 to PA", "count", "Number of live individuals; use area to standardize effort.",
  "UP_006", TRUE, TRUE, "standard 2-ha/20-min bird surveys; effort constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count", "Standardized bird count data.",
  "UP_007", TRUE, TRUE, "survey hours vary; area constant", "survey_effort_hours", "offset_abundance_covariate_pa", "retain raw counts; convert count > 0 to PA", "count", "Pollinator visit counts; standardize by survey duration.",
  "UP_008", TRUE, TRUE, "fixed-area belt transects; effort constant", "none", "none", "filter to 'number of individuals alive'; convert count > 0 to PA", "count", "Exclude total trees and % killed rows.",
  "UP_010", TRUE, TRUE, "non-standard search area and time vary", "survey_area + survey_effort_hours", "covariate", "filter to 'all plants (mature + juvenile)'; convert count > 0 to PA", "count", "Search-based plant counts; treat effort as covariates, not offsets.",
  "UP_011", TRUE, TRUE, "unconstrained flock search; time and area vary", "survey_area + survey_effort_hours", "covariate", "retain raw counts; convert count > 0 to PA", "count", "Relative encounter counts; effort covariates only.",
  "UP_012", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Tree survival/resprouting numerator-denominator data.",
  "UP_013", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Modeled koala density estimates; not raw observations.",
  "UP_014", TRUE, TRUE, "camera trap effort constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count_activity", "Camera trap detection/activity counts.",
  "UP_015", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Mixed methods, extreme area variation, incompatible response types.",
  "UP_017", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Modeled occupancy probabilities, not raw detections.",
  "UP_018", TRUE, FALSE, "visual survey hours vary strongly", "survey_effort_hours", "covariate", "use binary response directly", "binary_raw", "PA only; effort covariate for detectability.",
  "UP_020", TRUE, FALSE, "camera trap nights vary", "survey_effort_traps", "covariate", "use binary response directly", "binary_raw", "PA only; effort covariate.",
  "UP_021", TRUE, FALSE, "audio recording nights vary", "survey_effort_traps", "covariate", "use binary response directly", "binary_raw", "PA only; effort covariate.",
  "UP_024", TRUE, TRUE, "10-min point counts; effort constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count", "Standardized bird point counts.",
  "UP_025", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Detection rate not recoverable to raw detections.",
  "UP_026", TRUE, FALSE, "visual quadrat effort constant", "none", "none", "use binary response directly", "binary_raw", "Clean plant PA data.",
  "UP_027", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Aquatic/semi-aquatic eDNA detection data; detection and exposure are influenced by hydrology and downstream DNA transport, making it poorly comparable to site-level terrestrial biodiversity responses.",
  "UP_028", TRUE, FALSE, "maximum count across 1–3 breeding-season surveys; effort embedded in max-count response", "none", "none", "convert max count > 0 to PA", "binary_derived", "Counts represent maximum across repeated surveys, which biases abundance upward with number of surveys; converted to PA only and effort not modeled separately.",
  "UP_029", TRUE, FALSE, "aural survey effort constant", "none", "none", "use binary response directly", "binary_raw", "Calling male amphibian PA.",
  "UP_030", TRUE, FALSE, "camera trap nights vary", "survey_effort_traps", "covariate", "reconstruct detection nights = round(response_value * survey_effort_traps); convert > 0 to PA", "binary_derived", "Reporting rate recoverable to detection nights; PA only.",
  "UP_031", TRUE, FALSE, "visual survey hours vary slightly", "survey_effort_hours", "covariate", "use binary response directly", "binary_raw", "PA with effort covariate for consistency.",
  "UP_032", TRUE, TRUE, "pitfall trapping effort constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count_activity", "Standardized pitfall counts.",
  "UP_034", TRUE, FALSE, "vegetation plot effort constant", "none", "none", "convert cover > 0 to PA", "binary_derived", "Cover not abundance.",
  "UP_036", TRUE, TRUE, "trap effort constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count_activity", "Standardized trapping counts.",
  "UP_037", TRUE, FALSE, "observer effort varies but is not usable", "none", "none", "aggregate scat/run rows per site; convert summed value > 0 to PA", "binary_derived", "Indirect sign detections; PA only.",
  "UP_039", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Flower grazing/pollination process rates.",
  "UP_040", TRUE, TRUE, "transect length varies; survey_area is meters walked", "survey_area", "offset_abundance_covariate_pa", "count = round(response_value * (survey_area / 1000)); convert count > 0 to PA", "count_reconstructed", "Detections/km back-calculated to counts.",
  "UP_041", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Experimental seeding/restoration treatment.",
  "UP_042", TRUE, FALSE, "spotlight transect effort constant", "none", "none", "convert detection count > 0 to PA", "binary_derived", "Detection counts, not abundance.",
  "UP_043", TRUE, TRUE, "litter sampling effort constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count", "Standardized litter invertebrate counts.",
  "UP_045", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Aquatic fish electrofishing dataset; fire effects likely operate through stream/catchment hydrology, sediment, and riparian pathways not comparable to terrestrial site-level fire-speed exposure.",
  "UP_046", TRUE, FALSE, "PA survey effort constant", "none", "none", "use binary response directly", "binary_raw", "Clean PA data.",
  "UP_048", TRUE, FALSE, "timed search; response is non-integer density/rate", "none", "none", "convert response_value > 0 to PA", "binary_derived", "Not recoverable to counts.",
  "UP_050", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "cover", "Weed cover area, not biodiversity abundance/PA.",
  "UP_051", TRUE, FALSE, "bird transect effort constant", "none", "none", "convert detection count > 0 to PA", "binary_derived", "Detection counts, not abundance.",
  "UP_052", TRUE, FALSE, "visual search effort constant", "none", "none", "use binary response directly", "binary_raw", "Abundance classes ignored.",
  "UP_053", TRUE, FALSE, "pole camera effort roughly constant", "none", "none", "use binary response directly", "binary_raw", "Artificial hollow occupancy; include PA only with caution.",
  "UP_055", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "unknown", "Undefined response variable.",
  "UP_056", TRUE, TRUE, "seed count quadrats (0.2x0.1 m x 7 samples; area constant but mis-entered in raw data)", "none", "none", "set survey_area = 0.000014 for all rows; retain counts; convert count > 0 to PA", "count", "Sampling area was inconsistently recorded (0.002 vs 0.000014) but should be constant; corrected to 0.000014 ha for all observations.",
  "UP_059", TRUE, TRUE, "pitfall trap effort constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count_activity", "Trap captures/activity-density.",
  "UP_060", TRUE, FALSE, "camera trap nights vary", "survey_effort_traps", "covariate", "use binary response directly", "binary_raw", "PA only; effort covariate.",
  "UP_061", TRUE, FALSE, "acoustic survey effort constant", "none", "none", "use binary response directly", "binary_raw", "Clean acoustic PA.",
  "UP_062", TRUE, TRUE, "visual survey time constant; mixed response rows", "none", "none", "filter to 'nest occupancy and reproduction'; retain counts; convert count > 0 to PA", "count", "Drop binary nest-presence rows.",
  "UP_065", TRUE, TRUE, "survey area varies modestly; time constant", "survey_area", "offset_abundance_covariate_pa", "retain raw counts; convert count > 0 to PA", "count", "Small reptile count dataset; area correction included for consistency.",
  "UP_066", TRUE, FALSE, "survey duration varies strongly", "survey_effort_hours", "covariate", "use binary response directly", "binary_raw", "PA only; effort covariate.",
  "UP_067", TRUE, TRUE, "bird point count duration varies; area constant", "survey_effort_hours", "covariate", "retain raw counts; convert count > 0 to PA", "count_detection", "Use survey duration as covariate in both models.",
  "UP_068", TRUE, TRUE, "drone survey area varies", "survey_area", "offset_abundance_covariate_pa", "count = as.numeric(response_variable_text); convert count > 0 to PA", "count_reconstructed", "Use counts, not density response_value; area offset for abundance.",
  "UP_069", TRUE, TRUE, "two repeated surveys; effort constant", "none", "none", "count = response_value * 2; convert count > 0 to PA", "count_reconstructed", "Average count converted back to total count.",
  "UP_070", TRUE, TRUE, "camera trap nights vary", "survey_effort_traps", "offset_abundance_covariate_pa", "count = round(response_value * survey_effort_traps); convert count > 0 to PA", "count_reconstructed_activity", "Detection proportions converted to detection counts.",
  "UP_071", TRUE, FALSE, "call playback effort constant", "none", "none", "use binary response directly", "binary_raw", "Clean PA data.",
  "UP_074", TRUE, FALSE, "camera trap effort constant", "none", "none", "use binary response directly", "binary_raw", "Camera trap PA; counts unavailable.",
  "UP_075", TRUE, TRUE, "500 m transect effort constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count_detection", "Transect detection counts treated as relative abundance.",
  "UP_078", TRUE, TRUE, "time and area effectively constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count", "Standardized invertebrate counts.",
  "UP_079", TRUE, FALSE, "survey duration varies; true area unrecorded", "survey_effort_hours", "covariate", "convert count > 0 to PA", "binary_derived", "Abundance excluded because area varied but was not recorded.",
  "UP_080", FALSE, FALSE, "exclude", NA_character_, "none", "exclude", "other", "Targeted non-random sweep sampling of flowering vegetation; no fixed patch size or effort.",
  "UP_082", TRUE, TRUE, "camera trap nights vary", "survey_effort_traps", "offset_abundance_covariate_pa", "count = as.numeric(response_variable_text); convert count > 0 to PA", "count_detection_frequency", "Use detection-day counts with trap-night effort.",
  "UP_083", TRUE, TRUE, "pitfall trap nights vary", "survey_effort_traps", "offset_abundance_covariate_pa", "retain raw counts; convert count > 0 to PA", "count_activity", "Large pitfall capture dataset; trap-night offset.",
  "UP_084", TRUE, TRUE, "fixed-duration visual walks; time constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count_time_standardized", "Timed butterfly/insect counts.",
  "UP_085", TRUE, TRUE, "fixed-duration visual walks; time constant", "none", "none", "retain raw counts; convert count > 0 to PA", "count_time_standardized", "Treat consistently with UP_084."
)


# =========================================================
# Clean AMFS + assign FIRED speed + build analysis datasets
# Outputs:
#   1) cleaned master dataset with FIRED metrics
#   2) analysis-ready presence/absence dataset
#   3) analysis-ready abundance dataset
# =========================================================

safe_num <- function(x) suppressWarnings(as.numeric(x))

make_pa <- function(x) {
  as.integer(!is.na(x) & x > 0)
}

# =========================================================
# 1. Clean AMFS dates, coordinates, burn status, and severity
# =========================================================

amfs_base <- amfs %>%
  mutate(
    # Dates
    date_survey = parse_date_time(date_survey, orders = c("mdy", "dmy", "ymd")) %>% as.Date(),
    date_previous_record = parse_date_time(date_previous_record, orders = c("mdy", "dmy", "ymd")) %>% as.Date(),
    
    # UP_056 area correction
    survey_area = if_else(
      unique_project_ID == "UP_056",
      0.000014,
      survey_area
    ),
    
    # Clean labels
    before_after_clean = str_squish(as.character(before_after)),
    burnt_unburnt_clean = str_squish(as.character(burnt_unburnt)),
    
    # Coordinate decimal precision
    lat_decimals = if_else(
      is.na(latitude) | !str_detect(as.character(latitude), "\\."),
      0L,
      nchar(str_extract(as.character(latitude), "(?<=\\.)\\d+"))
    ),
    lon_decimals = if_else(
      is.na(longitude) | !str_detect(as.character(longitude), "\\."),
      0L,
      nchar(str_extract(as.character(longitude), "(?<=\\.)\\d+"))
    )
  ) %>%
  filter(
    latlongres == "original",
    !is.na(latitude),
    !is.na(longitude),
    lat_decimals >= 3,
    lon_decimals >= 3
  ) %>%
  mutate(
    is_before = before_after_clean == "Before Fire",
    is_after  = before_after_clean == "After Fire",
    is_burnt_raw = burnt_unburnt_clean == "Burnt",
    is_unburnt_raw = burnt_unburnt_clean == "Unburnt",
    
    # Important:
    # before-fire observations are treated as baseline/control,
    # even if original burnt_unburnt == "Burnt"
    burn_status_clean = case_when(
      is_before ~ "Unburnt",
      is_after & is_burnt_raw ~ "Burnt",
      is_after & is_unburnt_raw ~ "Unburnt",
      TRUE ~ NA_character_
    ),
    
    # Only after-fire burnt observations can receive fire speed
    eligible_fire_speed = burn_status_clean == "Burnt",
    
    fire_severity_clean = case_when(
      burn_status_clean == "Unburnt" ~ "Unburnt/None",
      burn_status_clean == "Burnt" ~ as.character(fire_severity),
      TRUE ~ NA_character_
    ),
    
    lat_r = round(latitude, 5),
    lon_r = round(longitude, 5),
    
    fired_key = paste(
      unique_project_ID,
      site_name,
      burn_status_clean,
      lat_r,
      lon_r,
      sep = "__"
    )
  )

# Coordinate filtering audit
coord_audit <- amfs %>%
  count(latlongres, name = "n_raw")

print(coord_audit)

amfs_base %>%
  summarise(
    n_after_coord_filter = n(),
    min_lat_decimals = min(lat_decimals, na.rm = TRUE),
    min_lon_decimals = min(lon_decimals, na.rm = TRUE)
  ) %>%
  print()


# =========================================================
# 2. Read and filter FIRED daily polygons
# =========================================================

fired_daily <- st_read(
  "data/raw/fired/fired_australia_2000_to_2024_daily.gpkg",
  layer = "fired_australia_2000_to_2024_daily"
)

start_date <- as.POSIXct("2019-06-01", tz = "UTC")
end_date   <- as.POSIXct("2020-05-31 23:59:59", tz = "UTC")

fired_1920_window <- fired_daily %>%
  filter(date >= start_date & date <= end_date) %>%
  mutate(date_day = as.Date(date))

days <- sort(unique(fired_1920_window$date_day))

fired_1920_window %>%
  st_drop_geometry() %>%
  summarise(
    min_date = min(date_day, na.rm = TRUE),
    max_date = max(date_day, na.rm = TRUE),
    n_polygons = n(),
    n_events = n_distinct(id)
  ) %>%
  print()


# =========================================================
# 3. Create eligible burnt locations for FIRED matching
# =========================================================

pts_burnt <- amfs_base %>%
  filter(eligible_fire_speed) %>%
  distinct(
    fired_key,
    unique_project_ID,
    site_name,
    burn_status_clean,
    lat_r,
    lon_r
  )

pts_burnt_sf <- st_as_sf(
  pts_burnt,
  coords = c("lon_r", "lat_r"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(st_crs(fired_1920_window)) %>%
  mutate(
    burn_date  = as.Date(NA),
    fire_id    = NA_real_,
    did        = NA_character_,
    fsr_km2_dy = NA_real_,
    dy_ar_km2  = NA_real_,
    match_type = NA_character_
  )


# =========================================================
# 4. Assign FIRED metrics by exact daily polygon intersection
# =========================================================

for (i in seq_along(days)) {
  d <- days[i]
  
  idx_unassigned <- which(is.na(pts_burnt_sf$burn_date))
  if (length(idx_unassigned) == 0) break
  
  pts_unassigned <- pts_burnt_sf[idx_unassigned, ]
  
  polys_day <- fired_1920_window %>% filter(date_day == d)
  if (nrow(polys_day) == 0) next
  
  hits <- st_intersects(pts_unassigned, polys_day, sparse = TRUE)
  hit_pts_local <- which(lengths(hits) > 0)
  if (length(hit_pts_local) == 0) next
  
  chosen_poly_local <- vapply(hit_pts_local, function(j) {
    poly_idx <- hits[[j]]
    o <- order(
      polys_day$fsr_km2_dy[poly_idx],
      polys_day$dy_ar_km2[poly_idx],
      decreasing = TRUE,
      na.last = TRUE
    )
    poly_idx[o[1]]
  }, integer(1))
  
  hit_pts_global <- idx_unassigned[hit_pts_local]
  
  pts_burnt_sf$burn_date[hit_pts_global]  <- d
  pts_burnt_sf$fire_id[hit_pts_global]    <- polys_day$id[chosen_poly_local]
  pts_burnt_sf$did[hit_pts_global]        <- polys_day$did[chosen_poly_local]
  pts_burnt_sf$fsr_km2_dy[hit_pts_global] <- polys_day$fsr_km2_dy[chosen_poly_local]
  pts_burnt_sf$dy_ar_km2[hit_pts_global]  <- polys_day$dy_ar_km2[chosen_poly_local]
  pts_burnt_sf$match_type[hit_pts_global] <- "intersect"
  
  if (i %% 10 == 0) {
    message(
      "Processed day ", i, "/", length(days),
      " (", d, ") | remaining unmatched: ",
      sum(is.na(pts_burnt_sf$burn_date))
    )
  }
}


# =========================================================
# 5. Recover unmatched eligible burnt locations within 1 km
# =========================================================

pts_unmatched_sf <- pts_burnt_sf %>%
  filter(is.na(burn_date))

if (nrow(pts_unmatched_sf) > 0) {
  
  nearest_idx <- st_nearest_feature(pts_unmatched_sf, fired_1920_window)
  nearest_polys <- fired_1920_window[nearest_idx, ]
  nearest_dist_m <- st_distance(pts_unmatched_sf, nearest_polys, by_element = TRUE)
  
  pts_recovered <- pts_unmatched_sf %>%
    mutate(
      nearest_dist_m = as.numeric(nearest_dist_m),
      burn_date = nearest_polys$date_day,
      fire_id = nearest_polys$id,
      did = nearest_polys$did,
      fsr_km2_dy = nearest_polys$fsr_km2_dy,
      dy_ar_km2 = nearest_polys$dy_ar_km2,
      match_type = "nearest_1km"
    ) %>%
    filter(nearest_dist_m <= 1000)
  
} else {
  pts_recovered <- pts_unmatched_sf %>%
    mutate(nearest_dist_m = numeric())
}


# =========================================================
# 6. Combine FIRED metrics
# =========================================================

pt_metrics_exact <- pts_burnt_sf %>%
  filter(!is.na(burn_date)) %>%
  st_drop_geometry() %>%
  select(
    fired_key,
    burn_date,
    fire_id,
    did,
    fsr_km2_dy,
    dy_ar_km2,
    match_type
  )

pt_metrics_recovered <- pts_recovered %>%
  st_drop_geometry() %>%
  select(
    fired_key,
    burn_date,
    fire_id,
    did,
    fsr_km2_dy,
    dy_ar_km2,
    match_type,
    any_of("nearest_dist_m")
  )

pt_metrics <- bind_rows(
  pt_metrics_exact,
  pt_metrics_recovered
) %>%
  arrange(fired_key, match_type) %>%
  distinct(fired_key, .keep_all = TRUE)


# =========================================================
# 7. Join FIRED metrics back to AMFS master
# =========================================================

amfs_master <- amfs_base %>%
  left_join(pt_metrics, by = "fired_key") %>%
  mutate(
    fire_speed = if_else(
      eligible_fire_speed,
      fsr_km2_dy,
      NA_real_
    ),
    
    fire_daily_area_km2 = if_else(
      eligible_fire_speed,
      dy_ar_km2,
      NA_real_
    ),
    
    days_since_fire = case_when(
      eligible_fire_speed & !is.na(burn_date) & !is.na(date_survey) ~
        as.integer(date_survey - burn_date),
      TRUE ~ NA_integer_
    ),
    
    timing_rel_fire = case_when(
      !eligible_fire_speed ~ "Control_or_before",
      is.na(burn_date) | is.na(date_survey) ~ NA_character_,
      days_since_fire < 0 ~ "Before",
      days_since_fire == 0 ~ "DayOf",
      days_since_fire > 0 ~ "After"
    ),
    
    # Model-safe: controls/before rows are 0 so they are not dropped
    fire_speed_model = if_else(
      eligible_fire_speed & !is.na(fire_speed),
      fire_speed,
      NA_real_
    ),
    
    has_fire_speed = eligible_fire_speed & !is.na(fire_speed)
  ) %>%
  left_join(project_review, by = "unique_project_ID")

# Cap days since fire at 0 (0.3% of values were down to -7)
amfs_master <- amfs_master %>%
  mutate(
    days_since_fire_raw = days_since_fire,
    days_since_fire = if_else(
      eligible_fire_speed & !is.na(days_since_fire) & days_since_fire < 0,
      0L,
      days_since_fire
    ),
    timing_rel_fire = case_when(
      !eligible_fire_speed ~ "Control_or_before",
      is.na(burn_date) | is.na(date_survey) ~ NA_character_,
      days_since_fire_raw < 0 ~ "Near_fire_date",
      days_since_fire == 0 ~ "DayOf",
      days_since_fire > 0 ~ "After"
    )
  )

# Create before/after control/impact column as a control variable
amfs_master <- amfs_master %>%
  mutate(
    study_design = case_when(
      is_before ~ "Before",
      is_after & burn_status_clean == "Burnt" ~ "After_Burnt",
      is_after & burn_status_clean == "Unburnt" ~ "After_Unburnt",
      TRUE ~ NA_character_
    )
  )

# Make it a factor
amfs_master <- amfs_master %>%
  mutate(
    study_design = factor(
      study_design,
      levels = c("Before", "After_Unburnt", "After_Burnt")
    )
  )

# Remove the few NAs
amfs_master <- amfs_master %>%
  filter(!is.na(study_design))

# Create unique site ID (unique to each project)
amfs_master <- amfs_master %>%
  mutate(
    site_id = paste(unique_project_ID, site_name, sep = "__"),
    site_id = factor(site_id)
  )

# =========================================================
# 8. FIRED / burn-status diagnostics
# =========================================================

table(amfs_master$before_after_clean, amfs_master$burnt_unburnt_clean, useNA = "ifany")
table(amfs_master$burn_status_clean, useNA = "ifany")
table(amfs_master$eligible_fire_speed, !is.na(amfs_master$fire_speed), useNA = "ifany")

amfs_master %>%
  filter(burn_status_clean == "Unburnt") %>%
  count(fire_severity_clean) %>%
  print()

match_summary <- amfs_master %>%
  distinct(fired_key, unique_project_ID, eligible_fire_speed, fire_speed) %>%
  filter(eligible_fire_speed) %>%
  group_by(unique_project_ID) %>%
  summarise(
    n_eligible_burnt_locations = n(),
    n_matched = sum(!is.na(fire_speed)),
    match_rate = n_matched / n_eligible_burnt_locations,
    .groups = "drop"
  ) %>%
  arrange(match_rate)

print(match_summary)

match_summary_analysis <- amfs_clean %>%
  distinct(fired_key, unique_project_ID, eligible_fire_speed, fire_speed) %>%
  filter(eligible_fire_speed) %>%
  group_by(unique_project_ID) %>%
  summarise(
    n_eligible_burnt_locations = n(),
    n_matched = sum(!is.na(fire_speed)),
    match_rate = n_matched / n_eligible_burnt_locations,
    .groups = "drop"
  ) %>%
  arrange(match_rate)

print(match_summary_analysis)
# =========================================================
# 9. Check project review decisions
# =========================================================

missing_review <- amfs_master %>%
  filter(is.na(include_pa), is.na(include_abundance)) %>%
  distinct(unique_project_ID, project_name)

if (nrow(missing_review) > 0) {
  message("Some projects are missing from project_review. They will not be included in PA/abundance datasets.")
  print(missing_review, n = Inf)
}


# =========================================================
# 10. Project-specific row filters
# =========================================================

amfs_filtered <- amfs_master %>%
  filter(include_pa | include_abundance) %>%
  filter(
    !(unique_project_ID == "UP_008" &
        !str_detect(str_to_lower(response_variable_text), "alive")),
    
    !(unique_project_ID == "UP_010" &
        !str_detect(str_to_lower(response_variable_text), "all plants|mature|juvenile")),
    
    !(unique_project_ID == "UP_062" &
        response_variable_text != "nest occupancy and reproduction")
  )


# =========================================================
# 11. Create standardized response_count and response_pa
# =========================================================

amfs_clean <- amfs_filtered %>%
  mutate(
    response_count = case_when(
      unique_project_ID == "UP_030" ~ round(response_value * survey_effort_traps),
      unique_project_ID == "UP_040" ~ round(response_value * (survey_area / 1000)),
      unique_project_ID == "UP_068" ~ safe_num(response_variable_text),
      unique_project_ID == "UP_069" ~ response_value * 2,
      unique_project_ID == "UP_070" ~ round(response_value * survey_effort_traps),
      unique_project_ID == "UP_082" ~ safe_num(response_variable_text),
      TRUE ~ response_value
    ),
    
    response_pa = case_when(
      response_type %in% c("binary", "binary_raw") ~ as.integer(response_value > 0),
      TRUE ~ make_pa(response_count)
    )
  )


# =========================================================
# 12. Create modeling effort columns
# =========================================================

amfs_clean <- amfs_clean %>%
  mutate(
    effort_offset_abundance = case_when(
      unique_project_ID == "UP_005" ~ survey_area,
      unique_project_ID == "UP_007" ~ survey_effort_hours,
      unique_project_ID == "UP_040" ~ survey_area / 1000,
      unique_project_ID == "UP_065" ~ survey_area,
      unique_project_ID == "UP_068" ~ survey_area,
      unique_project_ID == "UP_070" ~ survey_effort_traps,
      unique_project_ID == "UP_082" ~ survey_effort_traps,
      unique_project_ID == "UP_083" ~ survey_effort_traps,
      TRUE ~ 1
    ),
    
    effort_covariate_abundance = case_when(
      unique_project_ID == "UP_010" ~ survey_effort_hours,
      unique_project_ID == "UP_011" ~ survey_effort_hours,
      unique_project_ID == "UP_067" ~ survey_effort_hours,
      TRUE ~ 1
    ),
    
    effort_covariate2_abundance = case_when(
      unique_project_ID == "UP_010" ~ survey_area,
      unique_project_ID == "UP_011" ~ survey_area,
      TRUE ~ 1
    ),
    
    effort_covariate_pa = case_when(
      unique_project_ID %in% c("UP_018", "UP_031", "UP_066", "UP_067", "UP_079") ~ survey_effort_hours,
      unique_project_ID == "UP_028" ~ 1,
      unique_project_ID %in% c("UP_020", "UP_021", "UP_030", "UP_060", "UP_070", "UP_082", "UP_083") ~ survey_effort_traps,
      unique_project_ID == "UP_040" ~ survey_area / 1000,
      unique_project_ID %in% c("UP_005", "UP_065", "UP_068") ~ survey_area,
      unique_project_ID %in% c("UP_007", "UP_010", "UP_011") ~ survey_effort_hours,
      TRUE ~ 1
    ),
    
    effort_covariate2_pa = case_when(
      unique_project_ID == "UP_010" ~ survey_area,
      unique_project_ID == "UP_011" ~ survey_area,
      TRUE ~ 1
    )
  ) %>%
  mutate(
    effort_offset_abundance = if_else(is.na(effort_offset_abundance) | effort_offset_abundance <= 0, 1, effort_offset_abundance),
    effort_covariate_abundance = if_else(is.na(effort_covariate_abundance) | effort_covariate_abundance <= 0, 1, effort_covariate_abundance),
    effort_covariate2_abundance = if_else(is.na(effort_covariate2_abundance) | effort_covariate2_abundance <= 0, 1, effort_covariate2_abundance),
    effort_covariate_pa = if_else(is.na(effort_covariate_pa) | effort_covariate_pa <= 0, 1, effort_covariate_pa),
    effort_covariate2_pa = if_else(is.na(effort_covariate2_pa) | effort_covariate2_pa <= 0, 1, effort_covariate2_pa)
  )


# =========================================================
# 13. Remove rows where effort is required but missing
# =========================================================

amfs_clean <- amfs_clean %>%
  filter(
    !(unique_project_ID %in% c("UP_010", "UP_079") &
        include_pa == TRUE &
        is.na(survey_effort_hours))
  )

# =========================================================
# 13b. Remove burnt-after rows without fire speed (analysis only)
# =========================================================

# Audit BEFORE filtering
fire_missing_audit <- amfs_clean %>%
  filter(eligible_fire_speed, !has_fire_speed) %>%
  count(unique_project_ID, name = "n_missing_fire_speed") %>%
  arrange(desc(n_missing_fire_speed))

print(fire_missing_audit)

# Determine how many rows are lost from each project if they don't have fire speed
fire_missing_summary <- amfs_clean %>%
  filter(eligible_fire_speed) %>%
  group_by(unique_project_ID) %>%
  summarise(
    n_total = n(),
    n_missing = sum(!has_fire_speed),
    n_matched = sum(has_fire_speed),
    prop_missing = n_missing / n_total,
    .groups = "drop"
  ) %>%
  arrange(desc(prop_missing))

print(fire_missing_summary, n = Inf)

# Apply filter
amfs_clean <- amfs_clean %>%
  filter(
    !eligible_fire_speed | has_fire_speed
  )


# =========================================================
# 13c. Create shared modeling transformations
# =========================================================

amfs_clean <- amfs_clean %>%
  mutate(
    # Fire speed should only exist for matched burned-after rows
    fire_speed_model = if_else(
      eligible_fire_speed & has_fire_speed,
      fire_speed,
      NA_real_
    ),
    
    # Log fire speed only where biologically meaningful
    fire_speed_log = log(fire_speed_model),
    
    # Log effort variables
    log_effort_pa = log(effort_covariate_pa),
    log_effort2_pa = log(effort_covariate2_pa),
    
    log_effort_offset_abundance = log(effort_offset_abundance),
    log_effort_covariate_abundance = log(effort_covariate_abundance),
    log_effort_covariate2_abundance = log(effort_covariate2_abundance)
  ) %>%
  mutate(
    # Scaled shared predictors
    fire_speed_z = as.numeric(scale(fire_speed_log)),
    
    log_effort_pa_z = as.numeric(scale(log_effort_pa)),
    log_effort2_pa_z = as.numeric(scale(log_effort2_pa)),
    
    log_effort_offset_abundance_z = as.numeric(scale(log_effort_offset_abundance)),
    log_effort_covariate_abundance_z = as.numeric(scale(log_effort_covariate_abundance)),
    log_effort_covariate2_abundance_z = as.numeric(scale(log_effort_covariate2_abundance))
  )


# =========================================================
# 14. Build analysis-ready PA dataset
# =========================================================

pa_dat <- amfs_clean %>%
  filter(include_pa) %>%
  mutate(
    model_response = response_pa
  )


# =========================================================
# 15. Build analysis-ready abundance dataset
# =========================================================

abundance_dat <- amfs_clean %>%
  filter(include_abundance) %>%
  mutate(
    model_response = response_count,
    model_response = if_else(
      model_response %% 1 == 0,
      model_response,
      round(model_response)
    )
  )


# =========================================================
# 16. Final audits
# =========================================================

pa_bad_response <- pa_dat %>%
  filter(!model_response %in% c(0, 1) | is.na(model_response))

if (nrow(pa_bad_response) > 0) {
  print(pa_bad_response %>% count(unique_project_ID, model_response))
  stop("PA dataset contains non-binary or missing responses.")
}

abund_bad_response <- abundance_dat %>%
  filter(is.na(model_response) | model_response < 0 | !is.finite(model_response))

if (nrow(abund_bad_response) > 0) {
  print(abund_bad_response %>% count(unique_project_ID))
  stop("Abundance dataset contains missing, negative, or non-finite responses.")
}

stopifnot(all(pa_dat$effort_covariate_pa > 0))
stopifnot(all(pa_dat$effort_covariate2_pa > 0))
stopifnot(all(abundance_dat$effort_offset_abundance > 0))
stopifnot(all(abundance_dat$effort_covariate_abundance > 0))
stopifnot(all(abundance_dat$effort_covariate2_abundance > 0))

reconstructed_projects <- c("UP_030", "UP_040", "UP_069", "UP_070", "UP_082")

reconstructed_check <- abundance_dat %>%
  filter(unique_project_ID %in% reconstructed_projects) %>%
  mutate(diff_from_integer = abs(model_response - round(model_response))) %>%
  group_by(unique_project_ID) %>%
  summarise(
    n = n(),
    max_diff_from_integer = max(diff_from_integer, na.rm = TRUE),
    .groups = "drop"
  )

print(reconstructed_check)

pa_summary <- pa_dat %>%
  count(unique_project_ID, project_name, response_type, name = "n_pa_rows")

abundance_summary <- abundance_dat %>%
  count(unique_project_ID, project_name, response_type, name = "n_abundance_rows")

options(na.print = "NA")

pa_summary %>% as.data.frame() %>% print(row.names = FALSE)
abundance_summary %>% as.data.frame() %>% print(row.names = FALSE)

abund_not_pa <- project_review %>%
  filter(include_abundance, !include_pa)

if (nrow(abund_not_pa) > 0) {
  print(abund_not_pa)
  stop("Some abundance projects are not included in PA.")
}


# =========================================================
# 17. Save outputs
# =========================================================

dir.create("data/derived/cleaned_amfs", recursive = TRUE, showWarnings = FALSE)

write_csv(amfs_master, "data/derived/cleaned_amfs/amfs_master_clean_fired.csv")
write_csv(pa_dat, "data/derived/cleaned_amfs/amfs_presence_absence_clean.csv")
write_csv(abundance_dat, "data/derived/cleaned_amfs/amfs_abundance_clean.csv")
write_csv(pa_summary, "data/derived/cleaned_amfs/pa_project_summary.csv")
write_csv(abundance_summary, "data/derived/cleaned_amfs/abundance_project_summary.csv")
write_csv(match_summary, "data/derived/cleaned_amfs/fired_match_summary.csv")

message("Cleaning complete.")