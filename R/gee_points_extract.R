#' Extract data from GEE
#'
#' Extract data from GEE for buffered point locations for a set of covariates
#' described in a table
#'
#' @param meth_path path to table of extraction methods and sources. Will be
#'   updated with download progress. Must match BAM National Models format.
#'   https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn/edit?usp=sharing&ouid=104837701987164094932&rtpof=true&sd=true
#' @param loc.yr sf object with a unique row id and a year column identifying
#'   each unique site and year
#' @param out_save directory where outputs should be saved
#' @param projection.trans EPSG code for projection of output
#' @param local.buffer Buffer for smaller radius TODO: improve so this is only
#'   set in spreadsheet
#' @param neighbour.buffer Buffer for larger radius
#' @param do.test.run Should a test run of downloading 20 points from GEE be
#'   attempted first?
#' @param year_get Default is NA which means match the year column. If "max"
#'   then only the most recent year is extracted
#' @param gd_dl_dir Directory to download Google Drive folders to.
#' @param use_GD Use Google Drive directly as the source of covariates instead
#'   of the BAM method that assumes you have Google Drive linked to the local G
#'   folder. This requires use of the googledrive package which will require
#'   authentication the first time it is used.
#'
#' @return Files are saved to the out_save directory. And the spreadsheet at
#'   meth_path is updated with Complete = 1 and Running = 0 for variables that
#'   were downloaded successfully
#' @export
#'
#' @examples
gee_points_extract <- function(
    meth_path, loc.yr, out_save,
    projection.trans = 5072, local.buffer = 200,
    neighbour.buffer = 2000, do.test.run = TRUE, year_get = NA,
    gd_dl_dir = "data/raw-data/BAM_GoogleDrive_Covariates", use_GD = TRUE) {
  meth <- read.csv(meth_path,
    na.strings = c("NA", "")
  )

  # Stack Category missing for SCANFI so replace the NAs with unique numbers
  meth <- meth %>%
    mutate(StackCategory = ifelse(is.na(StackCategory), 1:n() + 900, StackCategory))

  if (!dir.exists(out_save)) dir.create(out_save)

  loc.n <- loc.yr %>%
    dplyr::mutate(year = as.integer(year)) %>%
    sf::st_transform(crs = projection.trans)

  loc.buff <- sf::st_buffer(loc.n, local.buffer)
  loc.buff2 <- sf::st_buffer(loc.n, neighbour.buffer)

  if (do.test.run) {
    loc.buff.small <- loc.buff %>% slice(1:20)
    gee_data <- ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")
    # Not working because SpatialResolution is a character
    # gee_data1<-ee_extract(gee_data, loc.buff.small, fun = ee$Reducer$mean(),
    #   scale = as.integer(meth.gee$SpatialResolution[3]))

    # Since National file says to use native I am extracting that from the object and using it.
    nominal_scale <- gee_data$projection()$nominalScale()$getInfo()

    # set scale to default for now
    gee_data1 <- ee_extract(gee_data, loc.buff.small,
      fun = ee$Reducer$mean(),
      scale = nominal_scale
    )
    message("Test run successful")
  }



  # 2.Extract Temporally Static #===============================================

  ## 2.1 Get list of static layers to run
  meth.gee <- dplyr::filter(
    meth, Source == "Google Earth Engine", Use == 1,
    TemporalResolution == "static", Complete == 0
  )

  if (!nrow(meth.gee) == 0) {
    ## 2.3. Create a plain dataframe
    loc.gee.static <- loc.n %>% sf::st_drop_geometry()

    ## 2.4. Make method loop
    for (i in 1:nrow(meth.gee)) {
      print(i)
      if (meth.gee$GEEtype[i] == "image") {
        if (!is.na(meth.gee$GEEBand[i])) {
          img.i <- ee$Image(meth.gee$Link[i])$
            select(meth.gee$GEEBand[i])
        } else {
          img.i <- ee$Image(meth.gee$Link[i])
        }
      }
      if (meth.gee$GEEtype[i] == "imagecollection") {
        img.i <- ee$ImageCollection(meth.gee$Link[i])$
          select(meth.gee$GEEBand[i])$
          toBands()
      }

      loc.gee.i <- switch(as.character(meth.gee$RadiusExtent[i]),
        "200" = loc.buff, # buffer width!
        "2000" = loc.buff2, # buffer width!
        "NA" = loc.n
      )

      n <- 1000 # number of calls per query/
      # GEE only accepts 5000 at time, but payload are restricted then 1000 was chosen
      loc.gee.i <- loc.gee.i %>%
        select(id) %>%
        group_by(row_number() %/% n) %>%
        dplyr::group_map(~.x)


      gee.data.static <- data.frame()

      for (j in 1:length(loc.gee.i)) {
        # Issue with ee_extract so I created a modified version. See Issue:
        # https://github.com/r-spatial/rgee/issues/367
        print("check mean or cv")
        if (meth.gee$RadiusFunction[i] == "mean") {
          gee.data.static.extract <- ee_extract(img.i, loc.gee.i[[j]],
            fun = ee$Reducer$mean(),
            scale = as.integer(meth.gee$GEEScale[i])
          )
        }

        if (meth.gee$RadiusFunction[i] == "cv") {
          gee.data.static.extract <- ee_extract(img.i, loc.gee.i[[j]],
            fun = ee$Reducer$stdDev(),
            scale = as.integer(meth.gee$GEEScale[i])
          )
        }

        if (ncol(gee.data.static.extract) < 2) { # Corrects for extraction when all data NA
          gee.data.static.extract <- gee.data.static.extract %>%
            mutate(namecol = NA)
          names(gee.data.static.extract)[2] <- meth.gee$GEEBand[i]
          message("all data is NA for: ", meth.gee$Link[i])
        }

        gee.data.static <- rbind(gee.data.static, gee.data.static.extract)
        print(paste0(j, " of ", length(loc.gee.i)))
      }
      print("left join")
      names(gee.data.static)[2] <- meth.gee$Label[i]
      loc.gee.static <- left_join(loc.gee.static, gee.data.static)

      print(paste0("Finished ", i, " of ", nrow(meth.gee)))

      # Update meth to track progress
      meth$Complete[which(meth$Label == meth.gee$Label[i])] <- 1
      meth$Running[which(meth$Label == meth.gee$Label[i])] <- 0
    }

    # Zerofill
    zerocols <- meth.gee %>%
      dplyr::filter(Zerofill == 1)
    loc.gee.static <- loc.gee.static %>%
      mutate_at(vars(zerocols$Label), ~ replace(., is.na(.), 0))


    write.csv(loc.gee.static,
      file = file.path(out_save, "data_covariates_GEE_static.csv"),
      row.names = FALSE
    )

    write.csv(meth, file = meth_path, row.names = FALSE)
  } else {
    message("All static variables have already been completed")
  }


  # 3. Extract Temporally matched #==============================================

  ## 3.1. Get list of temporally matched layers to run
  meth.gee <- dplyr::filter(
    meth, Source == "Google Earth Engine", Use == 1,
    TemporalResolution == "match", Complete == 0
  )


  if (!nrow(meth.gee) == 0) {
    ## 3.2. Plain dataframe for joining to output
    loc.gee.match <- loc.n %>% sf::st_drop_geometry()


    ## 3.3. Set up to loop

    for (i in 1:nrow(meth.gee)) {
      if (!is.na(year_get) && year_get == "max") {
        loc.gee.match$year <- meth.gee$GEEYearMax[i]
        meth.gee$GEEYearMin[i] <- meth.gee$GEEYearMax[i]
      }

      # Identify years of imagery
      years.gee <- seq(meth.gee$GEEYearMin[i], meth.gee$GEEYearMax[i])
      # considers having all years

      # Match year of data to year of data
      dt <- data.table::data.table(year = years.gee, val = years.gee)
      data.table::setattr(dt, "sorted", "year")
      data.table::setkey(dt, year)

      loc.n.i <- switch(as.character(meth.gee$RadiusExtent[i]),
        "200" = loc.buff, # buffer width!
        "2000" = loc.buff2, # buffer width!
        "NA" = loc.n
      )

      loc.n.i$yearrd <- dt[J(loc.n.i$year), roll = "nearest"]$val

      # Set up to loop through years
      loc.j <- data.frame()
      for (j in 1:length(years.gee)) {
        loc.n.yr <- dplyr::filter(loc.n.i, yearrd == years.gee[j]) %>%
          select(id)

        if (nrow(loc.n.yr) > 0) {
          # Set start & end date for image filtering---
          start.k <- paste0(
            years.gee[j] + meth.gee$YearMatch[i], "-",
            meth.gee$GEEMonthMin[i], "-01"
          )

          if (meth.gee$GEEMonthMax[i] > meth.gee$GEEMonthMin[i]) {
            end.k <- paste0(
              years.gee[j] + meth.gee$YearMatch[i], "-",
              meth.gee$GEEMonthMax[i], "-28"
            )
          }
          if (meth.gee$GEEMonthMax[i] < meth.gee$GEEMonthMin[i]) {
            end.k <- paste0(years.gee[j], "-", meth.gee$GEEMonthMax[i], "-28")
          }

          # Get the image
          img.i <- ee$ImageCollection(meth.gee$Link[i])$
            filter(ee$Filter$date(start.k, end.k))$select(meth.gee$GEEBand[i])$mean()

          # Extract
          if (nrow(loc.n.yr) > 1000) { # 5000 entries is the limit for google earth engine
            # Some maps exceed the payload size, therefore we selected 1000
            n <- 1000
            split <- loc.n.yr %>%
              group_by(row_number() %/% n) %>%
              group_map(~.x)
          }

          if (nrow(loc.n.yr) <= 1000) {
            split <- list(loc.n.yr)
          }

          gee.data <- data.frame()
          for (h in 1:length(split)) {
            if (meth.gee$Extraction[i] == "point") {
              gee.ext.data <- ee_extract(img.i, split[[h]],
                scale = as.integer(meth.gee$GEEScale[i])
              )
            }

            if (meth.gee$Extraction[i] == "radius") {
              if (meth.gee$RadiusFunction[i] == "mean") {
                gee.ext.data <- ee_extract(img.i, split[[h]],
                  fun = ee$Reducer$mean(),
                  scale = as.integer(meth.gee$GEEScale[i])
                )
              }

              if (meth.gee$RadiusFunction[i] == "mode") {
                gee.ext.data <- ee_extract(img.i, split[[h]],
                  fun = ee$Reducer$mode(),
                  scale = as.integer(meth.gee$GEEScale[i])
                )
              }
            }

            if (ncol(gee.ext.data) < 2) { # Corrects for extraction when all data NA
              gee.ext.data <- gee.ext.data %>%
                mutate(namecol = NA)
              names(gee.ext.data)[2] <- meth.gee$GEEBand[i]

              message("all data is NA for: ", meth.gee$Link[i], " and year ", years.gee[j])
            }
            gee.data <- rbind(gee.data, gee.ext.data)
          }

          # Collapse data across years
          names(gee.data)[2] <- meth.gee$Label[i]
          loc.j <- rbind(loc.j, gee.data)
        }

        print(paste0("Finished year ", j, " of ", length(years.gee)))
      }
      loc.gee.match <- left_join(loc.gee.match, loc.j)
      print(paste0("Finished ", i, " of ", nrow(meth.gee)))

      meth$Complete[which(meth$Label == meth.gee$Label[i])] <- 1
      meth$Running[which(meth$Label == meth.gee$Label[i])] <- 0
    }

    ## 3.4. Save
    write.csv(loc.gee.match,
      file = file.path(out_save, "data_covariates_GEE_match.csv"),
      row.names = FALSE
    )

    write.csv(meth, file = meth_path, row.names = FALSE)
  } else {
    message("All temporally matched variables have already been completed")
  }

  # Extract Google Drive #======================================================
  meth.gd <- meth %>% dplyr::filter(
    Source == "Google Drive" | Source == "SCANFI",
    Use == 1, Complete == 0
  )

  if (!nrow(meth.gd) == 0) {
    # data frame for adding output to
    loc.gd <- data.frame(loc.n) %>% dplyr::select(-geometry)

    # This loop is set so covariates can be stack together in the StackCategory column
    loop <- sort(unique(meth.gd$StackCategory))
    loop.i.count <- 0
    # was set to 38, don't know why...
    start_value <- 0

    for (i in loop) {
      if (i <= start_value) next
      # Filter to stack category
      meth.gd.i <- dplyr::filter(meth.gd, StackCategory == i)

      # Determine if temporally static
      if (meth.gd.i$TemporalResolution[1] == "static") {
        if (any(!is.na(meth.gd.i$LocalLink))) {
          meth.gd.i <- meth.gd.i %>%
            mutate(Link = ifelse(!is.na(meth.gd.i$LocalLink), meth.gd.i$LocalLink,
              meth.gd.i$Link
            ))
        }
        if (use_GD & any(str_detect(meth.gd.i$Link, "G:/Shared drives/BAM_NationalModels"))) {
          if (!dir.exists(gd_dl_dir)) dir.create(gd_dl_dir)
          meth.gd.i <- meth.gd.i %>%
            mutate(
              gd_path = Link %>%
                stringr::str_remove("G:/Shared drives/BAM_NationalModels/NationalModels5.0/CovariateRasters/"),
              gd_folder = gd_path %>% str_extract("(.*)/", group = 1),
              gd_file = gd_path %>% str_extract("/(.*)", group = 1)
            ) %>%
            rowwise() %>%
            mutate(
              LocalLink = gd_dl_file(gd_folder, gd_file, gd_dl_dir),
              Link = LocalLink
            )
        }

        rast.i <- terra::rast(meth.gd.i$Link)
        names(rast.i) <- meth.gd.i$Label
        loc.cov <- data.frame()
        loc.cov <- extract.gd(meth.gd.i, rast.i, loc.n, loc.cov)
      }
      # Determine if temporally matching
      else if (meth.gd.i$TemporalResolution[1] == "match") {
        # Get list of individual files
        files.i <- get.files(meth.gd.i, use_GD, gd_dl_dir)

        # if the raster is high resolution get rid of points outside the extent
        resolution.threshold <- 100
        meth.gd.i$SpatialResolution <- as.integer(meth.gd.i$SpatialResolution)
        if (meth.gd.i$SpatialResolution[1] < resolution.threshold &
          !is.na(meth.gd.i$SpatialResolution[1])) {
          rast.sample <- terra::rast(files.i$Link[1])
          loc.sample <- loc.n %>%
            sf::st_transform(terra::crs(rast.sample)) %>%
            terra::extract(x = rast.sample, ID = FALSE)
          colnames(loc.sample) <- "sample"

          # Rebuffer

          loc.n.loop <- cbind(loc.n, loc.sample) %>%
            dplyr::filter(!is.na(sample)) %>%
            select(-sample) %>%
            sf::st_transform(terra::crs(rast.sample))

          loc.buff.loop <- loc.n.loop %>%
            sf::st_transform(projection.trans) %>%
            sf::st_buffer(local.buffer) %>% # buffer width!
            sf::st_transform(terra::crs(rast.sample))

          loc.buff2.loop <- loc.n.loop %>%
            sf::st_transform(projection.trans) %>%
            sf::st_buffer(neighbour.buffer) %>% # buffer width!
            sf::st_transform(terra::crs(rast.sample))
        } else {
          loc.n.loop <- loc.n
          loc.buff.loop <- loc.buff
          loc.buff2.loop <- loc.buff2
        }

        # Match year of file to year of data
        # http://adomingues.github.io/2015/09/24/finding-closest-element-to-a-number-in-a-list/
        dt <- data.table::data.table(
          year = unique(files.i$year),
          val = unique(files.i$year)
        )
        data.table::setattr(dt, "sorted", "year")
        data.table::setkey(dt, year)

        # there is a buffer in covariate?
        loc.n.i <- switch(as.character(meth.gd.i$RadiusExtent[1]),
          "200" = loc.buff.loop, # buffer width!
          "2000" = loc.buff2.loop, # buffer width!
          "NA" = loc.n.loop
        )

        loc.n.i$year.rd <- dt[J(as.numeric(loc.n.i$year)), roll = "nearest"]$val

        # Set up to loop through years
        loc.cov <- data.frame()
        yrs <- unique(files.i$year)
        for (j in 1:length(yrs)) {
          # Read in raster in a given year
          files.j <- dplyr::filter(files.i, year == yrs[j])
          rast.i <- terra::rast(files.j$Link)
          names(rast.i) <- meth.gd.i$Label

          # NA out -9999s for Landfire
          if (str_detect(files.j$Link[1], "LandFire")) {
            rast.i <- subst(x = rast.i, from = -9999, to = NA)
          }

          # select the data of the year and extract data.
          loc.n.j <- dplyr::filter(loc.n.i, year.rd == yrs[j])
          loc.cov <- extract.gd(meth.gd.i, rast.i, loc.n.j, loc.cov)

          # NA out -9999s for Landfire I AM NOT SURE WHY I DID THIS I AM GOING TO SILENCE
          # if(str_detect(files.j$Link[1], "LandFire")){
          # loc.cov<- loc.cov %>%
          #  mutate(across(everything(), function(x){
          #   replace(x, which(x<=0), NA)}))
          # }

          print(paste0("Finished year ", j, " of ", length(yrs)))
        }
      }

      # Fix column names for CONUS
      if (any(str_detect(meth.gd.i$Link, "CONUS"))) {
        colnames(loc.cov) <- c(paste0(meth.gd.i$Label, "_conus"), "id")
        nms <- c(colnames(loc.gd), paste0(meth.gd.i$Label, "_conus"))
      } else {
        nms <- c(colnames(loc.gd), meth.gd.i$Label)
      }
      # Add output to main file
      loc.gd <- left_join(loc.gd, loc.cov)
      colnames(loc.gd) <- nms

      # Report status
      loop.i.count <- loop.i.count + 1
      print(paste0("Finished stack category ", loop.i.count, " of ", length(loop)))
      meth$Complete[meth$StackCategory == i] <- 1
      meth$Run[meth$StackCategory == i] <- 0
    }

    ## 2.5. Wrangle  and Save

    # Merge AK & CONUS columns for landfire if exist

    if ("LFheigthcv_1km_conus" %in% colnames(loc.gd)) {
      loc.gd <- loc.gd %>%
        mutate(
          LFbiomass_1km = ifelse(!is.na(LFbiomass_1km), LFbiomass_1km, LFbiomass_1km_conus),
          LFcrownclosure_1km = ifelse(!is.na(LFcrownclosure_1km), LFcrownclosure_1km, LFcrownclosure_1km_conus),
          LFheigth_1km = ifelse(!is.na(LFheigth_1km), LFheigth_1km, LFheigth_1km_conus),
          LFheigthcv_1km = ifelse(!is.na(LFheigthcv_1km), LFheigthcv_1km, LFheigthcv_1km_conus),
          LFbiomass_5x5 = ifelse(!is.na(LFbiomass_5x5), LFbiomass_5x5, LFbiomass_5x5_conus),
          LFcrownclosure_5x5 = ifelse(!is.na(LFcrownclosure_5x5), LFcrownclosure_5x5, LFcrownclosure_5x5_conus),
          LFheigth_5x5 = ifelse(!is.na(LFheigth_5x5), LFheigth_5x5, LFheigth_5x5_conus),
          LFheigthcv_5x5 = ifelse(!is.na(LFheigthcv_5x5), LFheigthcv_5x5, LFheigthcv_5x5_conus)
        ) %>%
        dplyr::select(
          -LFbiomass_1km_conus, -LFcrownclosure_1km_conus, -LFheigth_1km_conus, -LFheigthcv_1km_conus,
          -LFbiomass_5x5_conus, -LFcrownclosure_5x5_conus, -LFheigth_5x5_conus, -LFheigthcv_5x5_conus
        )
    }

    # Save

    write.csv(loc.gd,
      file = file.path(out_save, "data_covariates_GD.csv"),
      row.names = FALSE
    )

    write.csv(meth, file = meth_path, row.names = FALSE)

  } else {
    message("All Google Drive variables have already been completed")
  }
}

#' Download file from Google Drive
#'
#' Downloads the file or all files stored in the folder to the local hard drive
#' and returns their local paths.
#'
#' @param gd_folder Folder on Google Drive within the BAM CovariateRasters folder
#' @param gd_file File or folder within name within `gd_folder`. If it is a file
#'   it is expected to be a single tif if a folder it is a folder with multiple
#'   tif files for each year.
#' @param gd_dl_dir Local directory where the files should be downloaded to.
#'
#' @returns The local file path or a list of paths if `gd_file` is a folder.
#'
gd_dl_file <- function(gd_folder, gd_file, gd_dl_dir, rast_name = "") {
  file_path <- file.path(gd_dl_dir, gd_file)
  if (dir.exists(file_path)) {
    fls <- list.files(file_path, ".tif")
    if (length(fls) > 0) {
      return(file.path(file_path, fls) %>% list())
    }
  }
  if (file.exists(file_path) && !dir.exists(file_path)) {
    return(file_path)
  }
  file_out <- googledrive::drive_get("https://drive.google.com/drive/folders/11KojcXJ4sLSBONhpWGeRVlgPXfy8VDat") %>%
    googledrive::drive_ls() %>%
    filter(name == gd_folder) %>%
    googledrive::drive_ls() %>%
    filter(name == gd_file)
  if (file_out %>% googledrive::is_folder()) {
    new_dir <- file.path(gd_dl_dir, unique(file_out$name))
    if (!dir.exists(new_dir)) dir.create(new_dir)
    if (file_out %>% slice(1) %>% googledrive::is_folder() && file_out$name[1] == "SCANFI") {
      # Get file name to download
      file_name_pattern <- case_when(
        rast_name == "scanfilcc" ~ "VegTypeClass",
        rast_name == "conifer" ~ "prcC",
        rast_name == "deciduous" ~ "prcD",
        rast_name == "lodgepolepine" ~ "Lod.*polePine",
        rast_name == "scanfilcc" ~ "nfiLandCover",
        TRUE ~ rast_name
      )
      files_out <- file_out %>%
        googledrive::drive_ls() %>%
        filter(!str_detect(name, "DIFFUSION"), !str_detect(name, "nfiLandCover"))
      beepr::beep_on_error(files_out_list <- map(1:nrow(files_out), \(x){
        dl_df <- files_out %>%
          slice(x) %>%
          googledrive::drive_ls(pattern = ".tif$") %>%
          filter(str_detect(name, fixed(file_name_pattern, ignore_case = TRUE)))
        if (nrow(dl_df) == 0) {
          stop("no file matching ", file_name_pattern, " in folder ", gd_file,
            "/", files_out$name[x],
            call. = FALSE
          )
        }
        dl_df %>% rowwise() %>%
          # return local path if already downloaded
          mutate(
            LocalLink = ifelse(
              file.exists(file.path(file_path, name)),
              file.path(file_path, name),
              googledrive::drive_download(
                id,
                path = file.path(file_path, name), overwrite = TRUE
              ) %>%
                pull(local_path)
              # tried adding sleep here but didn't work
            )
          )
      }))
      files_out2 <- list_rbind(files_out_list)
      return(files_out2$LocalLink %>% list())
    }

    files_out <- file_out %>%
      googledrive::drive_ls() %>%
      rowwise() %>%
      mutate(
        LocalLink = googledrive::drive_download(
          id,
          path = file.path(file_path, name), overwrite = TRUE
        ) %>%
          pull(local_path)
      )
    return(files_out$LocalLink %>% list())
  } else {
    file_out %>%
      googledrive::drive_download(path = file_path, overwrite = TRUE) %>%
      pull(local_path) %>%
      return()
  }
}

# get.files is a function to get file structure for Landfire, NLCD and the rest
# of the covariates that are not a singular raster but multiple files. Within
# the function is consider that files can come from two different origins google
# drive repository (i.e., meth.gd.i$Link) and one local (meth.gd.i$LocalLink).
# Paths of both variables can be change in the spreadsheet. Be aware that files
# over stored in local are going to be prioritize
#
# The function read.files is going to bring the full path for each file that
# belongs to the folder and will assign it their respective year.

read.files <- function(meth.gd.i, use_GD, gd_dl_dir) {
  # For Landfire file structure
  if (any(str_detect(meth.gd.i$Label, "LF"))) {
    files <- data.frame(
      Link = list.files(meth.gd.i$Link,
        full.names = TRUE,
        recursive = TRUE, pattern = "*.tif"
      ),
      file = list.files(meth.gd.i$Link,
        recursive = TRUE,
        pattern = "*.tif"
      )
    ) %>%
      separate(file,
        into = c("year", "region", "name", "tif"),
        remove = FALSE
      ) %>%
      mutate(
        name = case_when(
          str_detect(name, c("CBD|cbd")) ~ "biomass",
          str_detect(name, c("CC|cc")) ~ "closure",
          str_detect(name, c("CH|ch")) ~ "height"
        ),
        year = as.numeric(year)
      ) %>%
      arrange(year) %>%
      unique()

    # Just height for CV extraction
    if (any(str_detect(meth.gd.i$Label, "cv"))) {
      files <- files %>%
        dplyr::filter(name == "height") # is going to just consider height (Caution)
    }
  }

  # For NLCD file structure
  else if (any(str_detect(meth.gd.i$Label, "NLCD"))) {
    files <- data.frame(
      Link = list.files(meth.gd.i$Link,
        full.names = TRUE,
        pattern = "*.img"
      ),
      file = list.files(meth.gd.i$Link, pattern = "*.img")
    ) %>%
      mutate(year = as.numeric(str_match(file, "cd_\\s*(.*?)\\s*_lan")[, 2])) %>%
      arrange(year) %>%
      unique() %>%
      dplyr::filter(!is.na(year)) # one NA is introduce as the file is 2001-2021
  }

  # For SCANFI
  else if (any(str_detect(meth.gd.i$Label, "SCANFI"))) {
    if (use_GD & any(str_detect(meth.gd.i$Link, "^G:"))) {
      files <- meth.gd.i %>%
        mutate(
          gd_folder = "Biomass",
          gd_file = "SCANFI"
        ) %>%
        rowwise() %>%
        mutate(
          LocalLink = gd_dl_file(gd_folder, gd_file, gd_dl_dir, RasterName),
          Link = list(LocalLink),
          file = list(map_chr(Link, \(x){
            str_extract(x, pattern = "^(.+/)*(.+\\.tif)$", group = 2)
          }))
        ) %>%
        select(Link, file, Name, RasterName) %>%
        unnest(cols = c(Link, file))
    } else {
      files <- data.frame(
        Link = list.files(meth.gd.i$Link,
          full.names = TRUE,
          recursive = TRUE, ignore.case = TRUE,
          pattern = "*tif$"
        ),
        file = list.files(meth.gd.i$Link,
          full.names = FALSE,
          recursive = TRUE, ignore.case = TRUE,
          pattern = "*tif$"
        )
      )
    }
    files <- files %>%
      dplyr::filter(
        !str_detect(file, "DIFFUSION"),
        !str_detect(file, "\\((.*?)\\)")
      ) %>%
      # DIFFUSION is for a repeated version on year 2020
      # then any file that has a copy (1), (2), etc.
      mutate(
        year = as.numeric(str_extract(file, "\\d+(?=\\_v)")),
        Name = case_when(
          str_detect(file, "VegTypeClass") ~ "scanfilcc",
          str_detect(file, "prcC") ~ "conifer",
          str_detect(file, "prcD") ~ "deciduous",
          str_detect(file, "prcB") ~ "deciduous",
          str_detect(file, "LodepolePine") ~ "lodgepolepine",
          str_detect(file, "nfiLandCover") ~ "scanfilcc",
          TRUE ~ file
        )
      ) %>%
      dplyr::filter(Name == meth.gd.i$RasterName |
        str_detect(file, fixed(meth.gd.i$RasterName,
          ignore_case = TRUE
        ))) %>%
      dplyr::select(Link, Name, year) %>%
      arrange(year, Name) %>%
      unique() %>%
      dplyr::filter(!is.na(year))
  } else {
    # Everything else
    if (use_GD & any(str_detect(meth.gd.i$Link, "G:/Shared drives/BAM_NationalModels"))) {
      files <- meth.gd.i %>%
        mutate(
          gd_path = Link %>%
            stringr::str_remove("G:/Shared drives/BAM_NationalModels/NationalModels5.0/CovariateRasters/"),
          gd_folder = gd_path %>% str_extract("(.*)/", group = 1),
          gd_file = gd_path %>% str_extract("/(.*)", group = 1)
        ) %>%
        rowwise() %>%
        mutate(
          LocalLink = gd_dl_file(gd_folder, gd_file, gd_dl_dir),
          Link = list(LocalLink),
          file = list(map_chr(Link, \(x){
            str_extract(x, pattern = "^(.+/)*(.+\\.tif)$", group = 2)
          }))
        )
      files <- files %>%
        select(Link, file) %>%
        unnest(cols = c(Link, file)) %>%
        mutate(year = as.numeric(str_extract(file, "[^_]+(?=\\.tif$)"))) %>%
        arrange(year) %>%
        unique() %>%
        dplyr::filter(!is.na(year))
    } else {
      files <- data.frame(
        Link = list.files(meth.gd.i$Link,
          full.names = TRUE,
          pattern = "*tif$", recursive = TRUE
        ),
        file = list.files(meth.gd.i$Link,
          pattern = "*tif$",
          recursive = TRUE
        )
      ) %>%
        mutate(year = as.numeric(str_extract(file, "[^_]+(?=\\.tif$)"))) %>%
        arrange(year) %>%
        unique() %>%
        dplyr::filter(!is.na(year))
    }
  }

  return(files)
}

get.files <- function(meth.gd.i, use_GD, gd_dl_dir) {
  if (any(!is.na(meth.gd.i$LocalLink))) {
    meth.gd.i.local <- meth.gd.i %>%
      mutate(Link = meth.gd.i$LocalLink)

    files.local <- read.files(meth.gd.i.local, use_GD, gd_dl_dir)

    if (any(!is.na(meth.gd.i$Link))) {
      files.gd <- read.files(meth.gd.i, use_GD, gd_dl_dir)

      files.final <- files.gd %>%
        mutate(Link = ifelse(files.gd$year == files.local$year, files.local$Link,
          files.gd$Link
        ))
    } else {
      files.final <- files.local
    }
  } else {
    files.final <- read.files(meth.gd.i, use_GD, gd_dl_dir)
  }
  return(files.final)
}

# Extract.gd is similar to the previous function but it is more general and it
# is going to use the Extraction.years. It is plan to extract the data from a
# rast.i and to distinguish between point or match extracion as well as betwen
# match and static temporal resolution. Static means that is not going to change
# over time, while for match we need to do a year temporal match with our visit
# data.
#
# In the second part of the code (second else) is a function we are going to
# extract using the loc.buffer created previously (in loop, line). This  part
# works only for covariates that need match in TemporalResolution, and finally
# those how had a buffer (i.e., 200 or 2000m).


extract.gd <- function(meth.gd.i, rast.i, loc.n.j, loc.cov) {
  message("extracting data from ", names(rast.i))
  if (any(str_detect(meth.gd.i$Extraction, "point"))) {
    loc.cov <- loc.n.j %>%
      sf::st_transform(terra::crs(rast.i)) %>%
      terra::extract(x = rast.i, ID = FALSE) %>% # extracts the values of the raster
      data.table::setnames(meth.gd.i$Label) %>%
      cbind(loc.n.j %>%
        sf::st_drop_geometry() %>%
        dplyr::select(id)) %>%
      rbind(loc.cov)
  }

  if (any(str_detect(meth.gd.i$Extraction, "radius"))) { # if it is a radius: Extract point value,
    # filter out NAs, reproject to UTM, buffer,
    # reproject to raster CRS, extract value,
    # rejoin to full df
    if (meth.gd.i$TemporalResolution[1] == "static") {
      loc.cov.i <- loc.n.j %>%
        sf::st_transform(terra::crs(rast.i)) %>%
        terra::extract(x = rast.i, bind = TRUE) %>%
        sf::st_as_sf() %>%
        dplyr::rename(cov = names(rast.i)[1]) %>%
        dplyr::filter(!is.na(cov)) %>%
        dplyr::select(-cov) %>%
        # sf::st_transform(crs=projection.trans) %>% # this doesn't do anything
        sf::st_buffer(meth.gd.i$RadiusExtent[1]) %>%
        sf::st_transform(terra::crs(rast.i))

      loc.cov <- loc.cov.i %>%
        exactextractr::exact_extract(
          x = rast.i,
          ifelse(meth.gd.i$RadiusFunction == "cv", "coefficient_of_variation", meth.gd.i$RadiusFunction),
          force_df = TRUE
        ) %>%
        data.table::setnames(meth.gd.i$Label) %>%
        cbind(loc.cov.i %>%
          sf::st_drop_geometry() %>%
          dplyr::select(id))
    } else {
      loc.cov <- loc.n.j %>%
        sf::st_transform(terra::crs(rast.i)) %>%
        exactextractr::exact_extract(x = rast.i, ifelse(meth.gd.i$RadiusFunction == "cv", "coefficient_of_variation", meth.gd.i$RadiusFunction), force_df = TRUE) %>%
        data.table::setnames(meth.gd.i$Label) %>%
        cbind(loc.n.j %>%
          sf::st_drop_geometry() %>%
          dplyr::select(id)) %>%
        rbind(loc.cov)
    }
  }
  loc.cov
}
