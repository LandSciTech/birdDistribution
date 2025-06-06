% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gee_points_extract.R
\name{gd_dl_file}
\alias{gd_dl_file}
\title{Download file from Google Drive}
\usage{
gd_dl_file(gd_folder, gd_file, gd_dl_dir, rast_name = "")
}
\arguments{
\item{gd_folder}{Folder on Google Drive within the BAM CovariateRasters folder}

\item{gd_file}{File or folder name within \code{gd_folder}. If it is a file it is
expected to be a single tif if a folder it is a folder with multiple tif
files for each year. Or for SCANFI a folder with multiple folders for each
year.}

\item{gd_dl_dir}{Local directory where the files should be downloaded to.}

\item{rast_name}{Only used for SCANFI to select which file within year folder
to extract}
}
\value{
The local file path or a list of paths if \code{gd_file} is a folder.
}
\description{
Downloads the file or all files stored in the folder to the local hard drive
and returns their local paths.
}
