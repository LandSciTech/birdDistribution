#' Download location and year information from WildTrax
#'
#' @param proj_id project_id in wildTrax
#' @param sens_id type of sensor, either "ARU", "CAM, or "PC"
#'
#' @return data frame of locations and years
#' @export
#'
wt_dl_loc_year <- function(proj_id, sens_id) {
  rec_id <- ifelse(sens_id == "ARU", "recording", "point_count")
  date_col <- ifelse(sens_id == "ARU", "recording_date_time", "survey_date")

  loc <- wildrtrax::wt_download_report(
    project_id = proj_id, sensor_id = sens_id,
    reports = c("location"), weather_cols = FALSE
  )

  if (is.null(nrow(loc))) {
    return(NULL)
  }

  loc <- loc %>%
    select(location_id, latitude, longitude) %>%
    unique()

  rec <- wt_download_report(
    project_id = proj_id, sensor_id = sens_id,
    reports = c(rec_id), weather_cols = FALSE
  ) %>%
    mutate(year = as.numeric(lubridate::year(.data[[date_col]]))) %>%
    select(project_id, location_id, year) %>%
    unique()

  dat <- full_join(loc, rec, by = join_by(location_id))
  dat
}
