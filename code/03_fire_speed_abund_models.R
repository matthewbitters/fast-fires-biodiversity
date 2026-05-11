### Fast fires impacts on biodiversity

### Modeling abundance

### Matt Bitters
### matthew.bitters@colorado.edu

# =========================================================
# Install and load required packages
# =========================================================

# Run this line in the terminal on cyverse
# conda install -c conda-forge r-rstan r-stanheaders r-rstantools r-rcpp r-rcppeigen --yes # run this in the terminal

# Run these in R
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)

install_cmdstan(cores = 8)

cmdstan_path()
cmdstan_version()


# =========================================================
# Install and load required packages
# =========================================================

packages <- c(
  "here", "dplyr", "readr", "stringr", "forcats",
  "brms", "posterior", "bayesplot", "tidybayes", "ggplot2", "rstan")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(dplyr)
library(readr)
library(stringr)
library(forcats)
library(brms)
library(posterior)
library(bayesplot)
library(tidybayes)
library(ggplot2)
library(rstan)

# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# =========================================================
# Set working directory
# =========================================================

setwd("/home/jovyan/data-store/fast-fires-biodiversity")
here::i_am("fast-fires-biodiversity.Rproj")


# ============================================================
# Read cleaned abundance dataset
# ============================================================

abund_dat <- read_csv(
  here("data", "derived", "cleaned_amfs", "amfs_abundance_clean.csv"),
  show_col_types = FALSE
)

glimpse(abund_dat)

summary(abund_dat$model_response)
table(abund_dat$study_design, useNA = "ifany")
table(abund_dat$broad_taxon, useNA = "ifany")


# ============================================================
# Final abundance modeling prep
# ============================================================

abund_mod_dat <- abund_dat %>%
  filter(
    !is.na(model_response),
    model_response >= 0,
    !is.na(study_design),
    !is.na(unique_project_ID),
    !is.na(site_id),
    !is.na(species),
    !is.na(broad_taxon),
    !is.na(ecosystem_group),
    !is.na(log_effort_offset_abundance),
    !is.na(log_effort_covariate_abundance_z),
    !is.na(log_effort_covariate2_abundance_z)
  ) %>%
  mutate(
    model_response = as.integer(round(model_response)),
    
    study_design = factor(
      study_design,
      levels = c("Before", "After_Unburnt", "After_Burnt")
    ),
    
    broad_taxon = fct_lump_min(as.factor(broad_taxon), min = 100),
    broad_taxon = fct_relevel(broad_taxon, "Bird", after = 0),
    
    ecosystem_group = as.factor(ecosystem_group),
    unique_project_ID = as.factor(unique_project_ID),
    site_id = as.factor(site_id),
    species = as.factor(species)
  )


# ============================================================
# Burned-only abundance dataset
# ============================================================

abund_fire_dat <- abund_mod_dat %>%
  filter(
    study_design == "After_Burnt",
    !is.na(fire_speed_z),
    !is.na(days_since_fire_z)
  ) %>%
  mutate(
    fire_severity_model = factor(
      fire_severity_model,
      levels = c(
        "Low (burnt understory, unburnt canopy)",
        "Moderate (partial canopy scorch)",
        "High (full canopy scorch)",
        "Extreme (full canopy consumption)"
      )
    )
  ) %>%
  filter(!is.na(fire_severity_model))


# ============================================================
# Quick checks
# ============================================================

abund_mod_dat %>%
  summarise(
    n = n(),
    n_projects = n_distinct(unique_project_ID),
    n_sites = n_distinct(site_id),
    n_species = n_distinct(species),
    n_taxa = n_distinct(broad_taxon),
    mean_count = mean(model_response),
    prop_zero = mean(model_response == 0),
    max_count = max(model_response)
  )

abund_fire_dat %>%
  summarise(
    n = n(),
    n_projects = n_distinct(unique_project_ID),
    n_sites = n_distinct(site_id),
    n_species = n_distinct(species),
    n_taxa = n_distinct(broad_taxon),
    mean_count = mean(model_response),
    prop_zero = mean(model_response == 0),
    fire_speed_min = min(fire_speed_model, na.rm = TRUE),
    fire_speed_max = max(fire_speed_model, na.rm = TRUE),
    days_since_fire_min = min(days_since_fire_model, na.rm = TRUE),
    days_since_fire_max = max(days_since_fire_model, na.rm = TRUE)
  )

table(abund_mod_dat$study_design, useNA = "ifany")
table(abund_fire_dat$fire_severity_model, useNA = "ifany")
table(abund_fire_dat$broad_taxon, useNA = "ifany")

# ============================================================
# More scaling for sensitivity model
# ============================================================

summary(abund_fire_dat$no_fires_model)
table(abund_fire_dat$no_fires_model)
summary(abund_fire_dat$unburnt5)
summary(abund_fire_dat$drought_spi_final)

# Make no_fires a factor, log and scale unburnt5, scale drought spi
abund_fire_dat <- abund_fire_dat %>%
  mutate(
    # Fire history
    no_fires_model = factor(no_fires_model),
    
    # Refugia
    unburnt5_log = log1p(unburnt5),
    unburnt5_z   = scale(unburnt5_log)[,1],
    
    # Drought
    drought_spi_z = scale(drought_spi_final)[,1]
  )

# ============================================================
# Priors
# ============================================================

abund_priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "shape")
)


# ============================================================
# Abundance Model 1: BACI baseline
# ============================================================

abund_m1_final <- brm(
  model_response ~
    study_design +
    ecosystem_group +
    offset(log_effort_offset_abundance) +
    log_effort_covariate_abundance_z +
    log_effort_covariate2_abundance_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = abund_mod_dat,
  family = negbinomial(link = "log"),
  prior = abund_priors,
  
  chains = 4,
  cores = 4,
  threads = threading(8),
  backend = "cmdstanr",
  iter = 6000,
  warmup = 3000,
  seed = 123,
  refresh = 100,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  
  file = here("models", "abund_m1_final")
)

summary(abund_m1_final)


# ============================================================
# Abundance Model 3: Fire speed + severity
# ============================================================

abund_m3_final <- brm(
  model_response ~
    fire_speed_z +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    offset(log_effort_offset_abundance) +
    log_effort_covariate_abundance_z +
    log_effort_covariate2_abundance_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = abund_fire_dat,
  family = negbinomial(link = "log"),
  prior = abund_priors,
  
  chains = 4,
  cores = 4,
  threads = threading(8),
  backend = "cmdstanr",
  iter = 6000,
  warmup = 3000,
  seed = 123,
  refresh = 100,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  
  file = here("models", "abund_m3_final")
)

summary(abund_m3_final)


# ============================================================
# Abundance Model 4: Fire speed x taxon
# ============================================================

abund_m4_final <- brm(
  model_response ~
    fire_speed_z * broad_taxon +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    offset(log_effort_offset_abundance) +
    log_effort_covariate_abundance_z +
    log_effort_covariate2_abundance_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = abund_fire_dat,
  family = negbinomial(link = "log"),
  prior = abund_priors,
  
  chains = 4,
  cores = 4,
  threads = threading(8),
  backend = "cmdstanr",
  iter = 6000,
  warmup = 3000,
  seed = 123,
  refresh = 100,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  
  file = here("models", "abund_m4_final")
)

summary(abund_m4_final)


# ============================================================
# Abundance Model 7: Fire speed x taxon + time since fire x taxon
# ============================================================

abund_m7_final <- brm(
  model_response ~
    fire_speed_z * broad_taxon +
    days_since_fire_z * broad_taxon +
    fire_severity_model +
    ecosystem_group +
    broad_taxon +
    offset(log_effort_offset_abundance) +
    log_effort_covariate_abundance_z +
    log_effort_covariate2_abundance_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = abund_fire_dat,
  family = negbinomial(link = "log"),
  prior = abund_priors,
  
  chains = 4,
  cores = 4,
  threads = threading(8),
  backend = "cmdstanr",
  iter = 6000,
  warmup = 3000,
  seed = 123,
  refresh = 100,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  
  file = here("models", "abund_m7_final")
)

summary(abund_m7_final)


# ============================================================
# Basic model checks
# ============================================================

pp_check(abund_m3_final)
pp_check(abund_m4_final)

bayes_R2(abund_m3_final)
bayes_R2(abund_m4_final)


################################################################################
################################################################################


### Final sensitivity checks

abund_m3_sens <- brm(
  model_response ~
    fire_speed_z +
    fire_severity_model +
    no_fires_model +
    unburnt5_z +
    drought_spi_z +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    offset(log_effort_offset_abundance) +
    log_effort_covariate_abundance_z +
    log_effort_covariate2_abundance_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = abund_fire_dat,
  family = negbinomial(link = "log"),
  prior = abund_priors,
  
  chains = 4,
  cores = 4,
  threads = threading(8),
  backend = "cmdstanr",
  iter = 6000,
  warmup = 3000,
  seed = 123,
  refresh = 100,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  
  file = here("models", "abund_m3_sens")
)

summary(abund_m3_final)

