#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table .BY
#' @importFrom data.table .EACHI
#' @importFrom data.table .GRP
#' @importFrom data.table .I
#' @importFrom data.table .N
#' @importFrom data.table .NGRP
#' @importFrom data.table .SD
#' @importFrom data.table :=
#' @importFrom data.table data.table
## usethis namespace: end
NULL


# funs_used <- funspotr::list_files_wd() %>%
#   rowwise() %>%
#   mutate(funs = list(funspotr::spot_funs_custom(pkgs = funspotr::spot_pkgs_from_description(DESCRIPTION_path = "DESCRIPTION"),
#                                            file_path = relative_paths)))
