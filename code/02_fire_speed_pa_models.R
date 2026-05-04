### Fast fires impacts on biodiversity

### Modeling presence/absence

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
# Final PA modeling prep
# ============================================================

pa_mod_dat <- pa_dat %>%
  filter(
    !is.na(model_response),
    !is.na(study_design),
    !is.na(unique_project_ID),
    !is.na(site_id),
    !is.na(species),
    !is.na(broad_taxon),
    !is.na(ecosystem_group),
    !is.na(log_effort_pa_z),
    !is.na(log_effort2_pa_z)
  ) %>%
  mutate(
    model_response = as.integer(model_response),
    
    study_design = factor(
      study_design,
      levels = c("Before", "After_Unburnt", "After_Burnt")
    ),
    
    broad_taxon = fct_lump_min(as.factor(broad_taxon), min = 100),
    broad_taxon = fct_relevel(broad_taxon, "Bird", after = 0),
    
    unique_project_ID = as.factor(unique_project_ID),
    site_id = as.factor(site_id),
    species = as.factor(species)
  )


# ============================================================
# Burned-only fire-speed dataset
# ============================================================

pa_fire_dat <- pa_mod_dat %>%
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
  )

pa_fire_dat <- pa_fire_dat %>%
  filter(!is.na(fire_severity_model))

# ============================================================
# Quick checks
# ============================================================

pa_mod_dat %>%
  summarise(
    n = n(),
    n_projects = n_distinct(unique_project_ID),
    n_sites = n_distinct(site_id),
    n_species = n_distinct(species),
    n_taxa = n_distinct(broad_taxon),
    prevalence = mean(model_response)
  )

pa_fire_dat %>%
  summarise(
    n = n(),
    n_projects = n_distinct(unique_project_ID),
    n_sites = n_distinct(site_id),
    n_species = n_distinct(species),
    prevalence = mean(model_response),
    fire_speed_min = min(fire_speed_model, na.rm = TRUE),
    fire_speed_max = max(fire_speed_model, na.rm = TRUE),
    days_since_fire_min = min(days_since_fire_model, na.rm = TRUE),
    days_since_fire_max = max(days_since_fire_model, na.rm = TRUE)
  )

table(pa_mod_dat$study_design, useNA = "ifany")
table(pa_fire_dat$fire_severity_model, useNA = "ifany")


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
    ecosystem_group +
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
# Model 2: Using burned data only - fire speed
# ============================================================

pa_m2 <- brm(
  model_response ~
    fire_speed_z +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
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


# ============================================================
# Model 3: Using burned data only - fire speed + fire severity
# ============================================================

pa_m3 <- brm(
  model_response ~
    fire_speed_z +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
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
  
  file = here("models", "pa_m3")
)

saveRDS(pa_m3, here("models", "pa_m3.rds"))

summary(pa_m3)


# ============================================================
# Model 4: Using burned data only - fire speed + fire severity + fire speed/taxon differences
# ============================================================

pa_m4 <- brm(
  model_response ~
    fire_speed_z * broad_taxon +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
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
  
  file = here("models", "pa_m4")
)

saveRDS(pa_m4, here("models", "pa_m4.rds"))

summary(pa_m4)


# ============================================================
# Model 4_1: Using burned data only - fire speed + fire severity + fire speed/taxon differences + number of previous fires
# ============================================================

pa_m4_1 <- brm(
  model_response ~
    fire_speed_z * broad_taxon +
    fire_severity_model +
    no_fires_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
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
  
  file = here("models", "pa_m4_1")
)

saveRDS(pa_m4_1, here("models", "pa_m4_1.rds"))

summary(pa_m4_1)

# ============================================================
# Model 5: Using burned data only - fire speed + fire severity + fire speed/ecosystem differences
# ============================================================

pa_m5 <- brm(
  model_response ~
    fire_speed_z * ecosystem_group +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
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
  
  file = here("models", "pa_m5")
)

saveRDS(pa_m5, here("models", "pa_m5.rds"))

summary(pa_m5)


# ============================================================
# Model 6: Using burned data only - fire speed + fire severity + fire speed/ecosystem differences + fire speed/taxon differences
# ============================================================

pa_m6 <- brm(
  model_response ~
    fire_speed_z * broad_taxon +
    fire_speed_z * ecosystem_group +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
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
  
  file = here("models", "pa_m6")
)

saveRDS(pa_m6, here("models", "pa_m6.rds"))

summary(pa_m6)


# ============================================================
# Loo comparison
# ============================================================


# Compute LOO
loo_pa_m1 <- loo(pa_m1)
loo_pa_m2 <- loo(pa_m2)
loo_pa_m3 <- loo(pa_m3)
loo_pa_m4 <- loo(pa_m4)
loo_pa_m4_1 <- loo(pa_m4_1)
loo_pa_m5 <- loo(pa_m5)
loo_pa_m6 <- loo(pa_m6)

# Save LOO objects
saveRDS(loo_pa_m1, here("models", "loo", "loo_pa_m1.rds"))
saveRDS(loo_pa_m2, here("models", "loo", "loo_pa_m2.rds"))
saveRDS(loo_pa_m3, here("models", "loo", "loo_pa_m3.rds"))
saveRDS(loo_pa_m4, here("models", "loo", "loo_pa_m4.rds"))
saveRDS(loo_pa_m4_1, here("models", "loo", "loo_pa_m4_1.rds"))
saveRDS(loo_pa_m5, here("models", "loo", "loo_pa_m5.rds"))
saveRDS(loo_pa_m6, here("models", "loo", "loo_pa_m6.rds"))

# Compare models (burned only dataset models)
loo_compare(
  loo_pa_m2,
  loo_pa_m3,
  loo_pa_m4,
  loo_pa_m4_1,
  loo_pa_m5,
  loo_pa_m6
)



# ============================================================
# Model 7: Using burned data only - m4_1 plus days since fire/taxon interaction
# ============================================================

pa_m7 <- brm(
  model_response ~
    fire_speed_z * broad_taxon +
    days_since_fire_z * broad_taxon +
    fire_severity_model +
    no_fires_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
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
  
  file = here("models", "pa_m7")
)

saveRDS(pa_m7, here("models", "pa_m7.rds"))

summary(pa_m7)


# Compute LOO
loo_pa_m7 <- loo(pa_m7)

# Save LOO objects
saveRDS(loo_pa_m7, here("models", "loo", "loo_pa_m7.rds"))

# Compare models (burned only dataset models)
loo_compare(
  loo_pa_m4_1,
  loo_pa_m7
)







################################################################################
################################################################################



### Re-run final models with more iterations




# Baseline BACI model
pa_m1_final <- brm(
  model_response ~
    study_design +
    ecosystem_group +
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
  iter = 6000,
  warmup = 3000,
  seed = 123,
  refresh = 100,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  
  file = here("models", "pa_m1_final")
)


# Fire severity/fire speed
pa_m3_final <- brm(
  model_response ~
    fire_speed_z +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
  family = bernoulli(link = "logit"),
  prior = pa_priors,
  
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
  
  file = here("models", "pa_m3_final")
)


# Taxon interaction w/ fire speed
pa_m4_final <- brm(
  model_response ~
    fire_speed_z * broad_taxon +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
  family = bernoulli(link = "logit"),
  prior = pa_priors,
  
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
  
  file = here("models", "pa_m4_final")
)


# Taxon recovery
pa_m7_final <- brm(
  model_response ~
    fire_speed_z * broad_taxon +
    days_since_fire_z * broad_taxon +
    fire_severity_model +
    days_since_fire_z +
    ecosystem_group +
    broad_taxon +
    log_effort_pa_z +
    log_effort2_pa_z +
    (1 | unique_project_ID) +
    (1 | site_id) +
    (1 | species),
  
  data = pa_fire_dat,
  family = bernoulli(link = "logit"),
  prior = pa_priors,
  
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
  
  file = here("models", "pa_m7_final")
)




# Summaries

summary(pa_m1_final)
summary(pa_m3_final)
summary(pa_m4_final)
summary(pa_m7_final)

pp_check(pa_m4_final)
