skip_if_not(interactive(), message = "gee tests need to be run interactively")
# library(rgee)

# initialize rgee
# Use this to save a token, doesn't need to be repeated every time
# rgee::ee$Authenticate(auth_mode='notebook')
rgee::ee$Initialize(project = "ee-sarahendicott-eccc") # <-- EDIT THIS FOR YOUR PROJECT

# Optionally
# make a request to verify you are connected.
# rgee::ee$String('Hello from the Earth Engine servers!')$getInfo()


# get test data set
test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))

# Format data to match BAM covariate download code
visit <- test_dat$all_surveys %>%
  sf::st_drop_geometry() %>%
  # rename to fit WT format
  dplyr::mutate(
    project_id = Project_Name, location_id = survey_ID,
    latitude = Latitude, longitude = Longitude,
    year = lubridate::year(Date_Time), Obs_Index, .keep = "none"
  )

visit_yr <- visit %>%
  dplyr::select(project_id, location_id, latitude, longitude, year) %>%
  dplyr::mutate(year = as.integer(year)) %>%
  unique() %>%
  dplyr::mutate(id = paste(project_id, location_id, latitude, longitude, year, sep = "_")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(crs = 5072)

# Get variable list from BAM google drive

googledrive::drive_auth(email = "sarah.endicott@canada.ca")

var_list_pth <- file.path(tempdir(), "NationalModels_V5_VariableList.xlsx")
googledrive::drive_download("https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn",
  path = var_list_pth, overwrite = TRUE
)

BAM_var_list <- readxl::read_xlsx(var_list_pth,
  sheet = "ExtractionLookup"
)

var_list_csv_pth <- file.path(tempdir(), "NationalModels_V5_VariableList.csv")

BAM_var_list %>%
  # Not currently working for some US and Alaska data
  filter(Extent %in% c("Canada", "global")) %>%
  group_by(Source, TemporalResolution) %>%
  slice(1) %>%
  mutate(Complete = 0) %>%
  write.csv(file = var_list_csv_pth, row.names = FALSE)


test_that("Data can download", {
  # takes ~6 mins to run
  out_dir <- file.path(tempdir(), "test_covars")
  tictoc::tic()
  # Function uses GEE and the BAM variables list to extract BAM variables that are
  # available from GEE including for year matched time series variables
  gee_points_extract(
    meth_path = var_list_csv_pth,
    loc.yr = visit_yr %>% group_by(year) %>% slice(1),
    out_save = out_dir,
    gd_dl_dir = file.path(tempdir(), "test_gd_dl"),
    do.test.run = FALSE
  )

  tictoc::toc()

  expect_length(list.files(out_dir), 3)

  expect_equal(read.csv(var_list_csv_pth) %>% pull(Complete) %>% min(), 1)
})
