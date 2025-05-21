#' Not in
#' @export
`%!in%` <- Negate(`%in%`)


#' Fit boosted regression tree models
#' Not currently in use...
fit_brt <- function(model_data,
                    response_column = NA,
                    covariate_columns = NA) {
  mod_brt <- NULL

  ntrees <- 50
  tcomplexity <- 5
  lrate <- 0.01
  m <- 0

  while (is.null(mod_brt)) {
    m <- m + 1
    if (m < 11) {
      ntrees <- 50
      lrate <- 0.01
    } else if (m < 21) {
      lrate <- 0.001
    } else if (m < 31) {
      ntrees <- 25
      lrate <- 0.001
    } else if (m < 41) {
      ntrees <- 25
      lrate <- 0.0001
    } else {
      break
    }

    ptm <- proc.time()
    if (inherits(try(
      mod_brt <- dismo::gbm.step(
        data = model_data,
        gbm.x = covariate_columns,
        gbm.y = response_column,
        offset = model_data$log_QPAD_offset,
        family = "poisson",
        tree.complexity = tcomplexity,
        learning.rate = lrate,
        n.trees = ntrees,
        n.folds = 5,
        max.trees = 10000
      )
    ), "try-error")) {
      cat("Couldn't fit model", n, "in the iteration", m, "\n")
    }
    t <- proc.time() - ptm
  }
  if (is.null(mod_brt)) {
    next
  }
  return(mod_brt)
}
