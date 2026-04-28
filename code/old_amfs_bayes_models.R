### Fast fires impacts on biodiversity

### This script filters the merged AMFS/FIRED dataset and runs some models.
### Matt Bitters
### matthew.bitters@colorado.edu


# Install and load required packages
packages <- c("here", "dplyr", "tidyr", "ggplot2", "brms", "rstan")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(brms)
library(rstan)

# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


################################################################################
################################################################################

### Read, factor, and filter data

# Read in merged AMFS/FIRED dataset
amfs <- read.csv("data/derived/amfs_analysis_ready.csv")
str(amfs)
head(amfs)
nrow(amfs)





amfs %>%
  group_by(unique_project_ID) %>%
  summarise(
    method = first(survey_method),
    response_variable = first(response_variable),
    min_val = min(response_value, na.rm=TRUE),
    max_val = max(response_value, na.rm=TRUE),
    n = n()
  ) %>%
  arrange(response_variable)



# Standardize response within project and taxon
amfs_std <- amfs %>%
  group_by(unique_project_ID, broad_taxon) %>%
  mutate(
    n_nonmiss = sum(!is.na(response_value)),
    sd_resp   = sd(response_value, na.rm = TRUE),
    response_std = if_else(
      n_nonmiss >= 2 & sd_resp > 0,
      (response_value - mean(response_value, na.rm = TRUE)) / sd_resp,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  select(-n_nonmiss, -sd_resp)

amfs_std_model <- amfs_std %>%
  filter(!is.na(response_std))

# Check to see how many are removed
amfs_std %>%
  summarise(
    n = n(),
    n_std_ok = sum(!is.na(response_std)),
    n_std_na = sum(is.na(response_std))
  )














# Filter to abundance (test)
dat_abund <- amfs %>%
  filter(response_variable == "presence/abundance")
nrow(dat_abund)


# Filter to abundance (test)
dat_abund <- amfs %>%
  filter(unique_project_ID == "UP_082") %>%
  filter(response_value <= 1)
nrow(dat_abund)

################################################################################
################################################################################

### Priors

# Negative binomial
priors_nb <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),              # fixed effects slopes
  prior(exponential(1), class = "sd"),             # group-level SDs
  prior(exponential(1), class = "shape")           # only for negbinomial()
)

# Binomial
priors_binom <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd")
)

# standardized gaussian 
priors_std <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

################################################################################
################################################################################

### Models

dat_abund$response_variable_text <- as.integer(dat_abund$response_variable_text)
unique(dat_abund$response_variable_text)
unique(dat_abund$survey_effort_traps)

# m1
m1 <- brm(
  formula = response_variable_text | trials(survey_effort_traps) ~ log1p(fsr_km2_dy) + fire_severity + no.burns +
    (1 | broad_taxon) +
    (1 | site_name),
  data = dat_abund,
  family = binomial,
  prior = priors_binom,
  chains = 4, 
  cores = 4,
  iter = 6000, 
  warmup = 3000,
  control = list(adapt_delta = 0.99),
  save_pars = save_pars(all = TRUE)
)

summary(m1)
ranef(m1)$broad_taxon
ranef(m1)$site_name


mean(dat_abund$response_value)
var(dat_abund$response_value)





# standardized reponse model

amfs_std_model <- amfs_std_model %>%
  mutate(
    fsr_z = scale(log1p(fsr_km2_dy)),
    burns_z = scale(no.burns),
    days_z = scale(log_days_since_fire)
  )


m_std <- brm(
  response_std ~ 
    fsr_z +
    burnt_unburnt +
    fire_severity +
    burns_z +
    days_z +
    (1 + fsr_z | unique_project_ID) +
    (1 | broad_taxon),
  data = amfs_std_model,
  family = gaussian(),
  prior = priors_std,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.99),
  save_pars = save_pars(all = TRUE)
)

summary(m_std)
ranef(m_std)$broad_taxon
ranef(m_std)$unique_project_ID


