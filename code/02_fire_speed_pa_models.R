### Fast fires impacts on biodiversity

### Modeling...

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
# Install and load required packages
# =========================================================

setwd("/home/jovyan/data-store/fast-fires-biodiversity")
here::i_am("fast-fires-biodiversity.Rproj")

# ============================================================
# Read cleaned PA dataset
# ============================================================

pa_dat <- read_csv(
  here("data", "derived", "cleaned_amfs", "amfs_presence_absence_clean.csv"),
  show_col_types = FALSE
)

glimpse(pa_dat)

table(pa_dat$model_response, useNA = "ifany")
table(pa_dat$study_design, useNA = "ifany")
table(pa_dat$broad_taxon, useNA = "ifany")


# ============================================================
# 3. Final modeling prep
# ============================================================

pa_mod_dat <- pa_dat %>%
  filter(
    !is.na(model_response),
    !is.na(study_design),
    !is.na(unique_project_ID),
    !is.na(site_name),
    !is.na(species),
    !is.na(broad_taxon)
  ) %>%
  mutate(
    
    # Factor structure
    study_design = factor(
      study_design,
      levels = c("Before", "After_Unburnt", "After_Burnt")
    ),
    
    broad_taxon = fct_lump_min(as.factor(broad_taxon), min = 100),
    broad_taxon = fct_relevel(broad_taxon, "Bird", after = 0),
    
    unique_project_ID = as.factor(unique_project_ID),
    site_id = as.factor(site_id),
    species = as.factor(species)
    
      ) %>%
  filter(
    !is.na(fire_speed_z),
    !is.na(log_effort_pa_z),
    !is.na(log_effort2_pa_z)
  )

# Quick checks
pa_mod_dat %>%
  summarise(
    n = n(),
    n_projects = n_distinct(unique_project_ID),
    n_sites = n_distinct(site_id),
    n_species = n_distinct(species),
    n_taxa = n_distinct(broad_taxon),
    prevalence = mean(model_response)
  )

table(pa_mod_dat$study_design, useNA = "ifany")
table(pa_mod_dat$broad_taxon, useNA = "ifany")


# ============================================================
# Priors
# ============================================================

pa_priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(exponential(1), class = "sd")
)



# ============================================================
# Model 1: core BACI/fire-speed model
# ============================================================

pa_m1 <- brm(
  model_response ~
    study_design +
    fire_speed_z +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_mod_dat,
  family = bernoulli(link = "logit"),
  prior = pa_priors,
  
  chains = 4,
  cores = 4,
  threads = threading(8),
  backend = "cmdstanr",
  iter = 2000,
  warmup = 1000,
  seed = 123,
  refresh = 100,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  
  file = here("models", "pa_m1")
)

saveRDS(pa_m1, here("models", "pa_m1.rds"))

summary(pa_m1)



# ============================================================
# Model 2: Study design/fire speed interaction
# ============================================================

pa_m2 <- brm(
  model_response ~
    study_design * fire_speed_z +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_mod_dat,
  family = bernoulli(link = "logit"),
  prior = pa_priors,
  
  chains = 4,
  cores = 4,
  threads = threading(8),
  backend = "cmdstanr",
  iter = 2000,
  warmup = 1000,
  seed = 123,
  refresh = 100,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  
  file = here("models", "pa_m2")
)

saveRDS(pa_m2, here("models", "pa_m2.rds"))

summary(pa_m2)
