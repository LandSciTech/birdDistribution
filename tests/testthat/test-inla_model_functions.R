library(tidyverse)

# get test data set
test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))

# Using km for projection units because it helps with INLA model calculations to have smaller numbers
AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "

test_area <- sf::st_convex_hull(sf::st_union(test_dat$all_surveys)) %>%
  sf::st_transform(AEA_proj)

atl_sq <- sf::st_make_grid(
  test_area,
  cellsize = units::set_units(5 * 5, km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE
) %>%
  sf::st_as_sf() %>%
  sf::st_intersection(test_area) %>%
  na.omit() %>%
  mutate(sq_id = 1:nrow(.)) %>%
  rename(geometry = x)

bcr_poly <- test_area %>% sf::st_buffer(-10000)

# Data frame of covariates with priors
sd_linear <- 1
cov_df <- data.frame(
  covariate = paste0("PC", 1:8), model = "linear",
  mean = 0, prec = 1 / sd_linear^2
)
cov_df <- cov_df %>%
  mutate(beta = 1) %>%
  bind_rows(cov_df %>% mutate(beta = 2, prec = 1 / (sd_linear / 2)^2))


test_that("overall flow works", {
  mod <- fit_inla(
    sp_code = test_dat$species_to_model$Species_Code_BSC[1],
    analysis_data = test_dat,
    proj_use = AEA_proj,
    study_poly = test_area,
    covariates = cov_df,
    mod_dir = tempdir(),
    save_mod = FALSE
  )

  pred <- predict_inla(
    dat = test_dat$ONGrid,
    analysis_data = test_dat,
    mod = mod,
    sp_code = test_dat$species_to_model$Species_Code_BSC[1],
    covariates = cov_df,
    do_crps = FALSE
  )

  dir_use <- tempdir()
  map_inla_preds(
    sp_code = test_dat$species_to_model$Species_Code_BSC[1],
    analysis_data = test_dat,
    pred,
    proj_use = AEA_proj,
    study_poly = test_area,
    atlas_squares = atl_sq %>% sf::st_transform(AEA_proj),
    bcr_poly = bcr_poly,
    map_dir = dir_use
  )
  if (interactive()) {
    list.files(dir_use, pattern = "png$", full.names = TRUE) %>% walk(shell.exec)
  }
})
