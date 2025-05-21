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
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_at
#' @importFrom dplyr vars
#' @importFrom dplyr full_join
#' @importFrom dplyr left_join
#' @importFrom dplyr join_by
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr bind_rows
#' @importFrom dplyr n
#' @importFrom dplyr row_number
#' @importFrom dplyr rowwise
#' @importFrom dplyr arrange
#' @importFrom dplyr case_when
#' @importFrom dplyr summarize
#' @importFrom dplyr relocate
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom purrr walk
#' @importFrom purrr list_rbind
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract
#' @importFrom stringr str_subset
#' @importFrom stringr str_replace
#' @importFrom stringr fixed
#' @importFrom tidyr unnest
## usethis namespace: end
NULL


# funs_used <- funspotr::list_files_wd() %>%
#   rowwise() %>%
#   mutate(funs = list(try(funspotr::spot_funs_custom(pkgs = funspotr::spot_pkgs_from_description(DESCRIPTION_path = "DESCRIPTION"),
#                                            file_path = relative_paths))))
# funs_used %>% filter(!is(funs, "try-error")) %>% unnest(funs) %>%
#   filter(!pkgs %in% c("base", "(unknown)")) %>%
#   select(-absolute_paths) %>%
#   filter(pkgs == "stringr") %>%
#   mutate(import_text = paste("@importFrom", pkgs, funs)) %>% pull(import_text) %>% paste0(collapse = "\n#' ") %>% cat()

