library(usethis)
# code to prepare datasets taken from `SpatialEpi`

library(SpatialEpi)
usethis::use_data(pennLC, overwrite = T)
usethis::use_data(scotland, overwrite = T)
usethis::use_data(NYleukemia, overwrite = T)

# code to prepare case-control data

## load individual-level training data
library(data.table)
library(dplyr)
library(sf)

individual.data_path <- 'data-raw/matchdata_isisegysc.csv'
survey = data.table::fread(file = individual.data_path)

survey_egypt <- survey %>%
  dplyr::select(case, coledu, age, married, student, lowstat,
         population_density, total_population_2006, christian_2006_pct, university_2006_pct,
         agriculture_2006_pct, mursi_vote_2012_pct, sqrt_killed_at_rabaa, unemployment_2013q4_pct,
         sqrt_protest_post_Mubarak, adm2_pcode)

usethis::use_data(survey_egypt, overwrite = T)

# load Egypt shape file
shapefile_path <- "data-raw/EGY_adm2.shp"
shape_egypt = sf::st_read(shapefile_path)
usethis::use_data(shape_egypt, overwrite = T)

# get denominators data for prevalence adjustment
denom_mena = data.table::fread(file = 'data-raw/mena_pops.csv')
usethis::use_data(denom_mena, overwrite = T)

##get data.list file for estimation example

data(survey_egypt)
data(shape_egypt)
data(denom_mena)
denom_egypt <- denom_mena$male_18_sunni[which(denom_mena$country=="Egypt")]

data.list =
  data.prep(
    shape = shape_egypt,
    survey = survey_egypt,
    shape_large.area_id_name = "ADM1_EN",
    shape_large.area_id_num = NA,
    shape_small.area_id_name = "ADM2_EN",
    shape_small.area_id_num = "ADM2_PCODE",
    survey_small.area_id_name = NA,
    survey_small.area_id_num = "adm2_pcode",
    drop.incomplete.records = T,
    colnames_X = c(
      "coledu",
      "age",
      "married",
      "student",
      "lowstat",
      "population_density",
      "total_population_2006",
      "christian_2006_pct",
      "university_2006_pct",
      "agriculture_2006_pct",
      "mursi_vote_2012_pct",
      "sqrt_killed_at_rabaa",
      "unemployment_2013q4_pct",
      "sqrt_protest_post_Mubarak"
    ),
    interactions_list = list(age2 = "age*age", coledu_lowstat = "coledu*lowstat"),
    scale_X = "1sd",
    colname_y = "case",
    contamination = T,
    pi = 1000 / denom_egypt,
    large_area_shape = T
  )

usethis::use_data(data.list, overwrite = T)

##get monitor object for fit() example

fit_object = fit(
  data = data.list,
  show_code = T,
  contamination = T,
  offset = T,
  beta_prior = "cauchy",
  small_area_prior = "BYM2",
  intercept_scale = 10,
  large_area_prior = "random",
  iter = 25000,
  warmup = 22500,
  thin = 4,
  cores = 4,
  chains = 4,
  control = list(max_treedepth = 25, adapt_delta = 0.99),
  verbose = T
)

mon_egypt <- rstan::monitor(rstan::extract(fit_object, permuted = FALSE, inc_warmup = TRUE))

usethis::use_data(mon_egypt, overwrite = T)
