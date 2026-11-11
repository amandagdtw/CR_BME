#' @title Leave-one-out cross validation (LOOCV) at hard data locations.
#'
#' @usage new_bme_cv(data_object, model, nugget, sill, range,
#'               nsmax = 5, nhmax = 5, n = 50,
#'               zk_range = BMEmapping:::extended_range(data_object), type)
#'
#' @param data_object A list containing the hard and soft data.
#' @param model A string specifying the variogram or covariance model to use
#'        (e.g., \code{"exp"}, \code{"sph"}, etc.).
#' @param nugget A non-negative numeric value for the nugget effect in the
#'        variogram model.
#' @param sill A numeric value representing the sill (total variance) in the
#'        variogram model.
#' @param range A positive numeric value for the range (or effective range)
#'        parameter of the variogram model.
#' @param nsmax An integer specifying the maximum number of nearby soft data
#'        points to include for estimation (default is 5).
#' @param nhmax An integer specifying the maximum number of nearby hard data
#'        points to include for estimation (default is 5).
#' @param n An integer indicating the number of points at which to evaluate the
#'        posterior density over \code{zk_range} (default is 50).
#' @param zk_range A numeric vector specifying the range over which to evaluate
#'        the unobserved value at the estimation location (\code{zk}). Although
#'        \code{zk} is unknown,  it is assumed to lie within a range similar to
#'        the observed data (\code{zh}, \code{a}, and \code{b}). It is advisable
#'        to explore the posterior distribution at a few locations using
#'        \code{prob_zk()} before finalizing this range.
#' @param type A string indicating the type of BME prediction to compute: either
#'        \code{"mean"} for the posterior mean or \code{"mode"} for the
#'        posterior mode.
#'
#' @returns A data frame containing the coordinates, observed
#'          values, BME predictions (posterior \code{mean} or \code{mode}),
#'          posterior variance (if \code{type = "mean"}), residuals, and fold
#'          indices.
#'
#' @description
#' \code{new_bme_cv} performs LOOCV to evaluate the prediction performance
#' of the Bayesian Maximum Entropy (BME) spatial interpolation method using both
#' hard and soft (interval) data. For each hard data location, the function
#' removes the observed value and predicts it using all remaining hard and soft
#' data points. This is repeated for every hard data location. The predictions
#' are either posterior means or posterior modes, depending on the \code{type}
#' argument.
#'
#' This function is useful for validating the BME interpolation method and
#' tuning variogram parameters.
#'
#' @examples
#' data("utsnowload")
#' ch <- utsnowload[2:10, c("latitude", "longitude")]
#' cs <- utsnowload[68:232, c("latitude", "longitude")]
#' zh <- utsnowload[2:10, c("hard")]
#' a <- utsnowload[68:232, c("lower")]
#' b <- utsnowload[68:232, c("upper")]
#' data_object <- bme_map(ch, cs, zh, a, b)
#' new_bme_cv(data_object, lambda = 1, type = "mean")
#'
#' @export
new_bme_cv <- function(data_object, nsmax = 5, lambda, mu, nhmax = 5, n = 50,
                       zk_range = BMEmapping:::extended_range(data_object), type) {

  type <- match.arg(type, choices = c("mean", "mode"))
  col_idx <- if (type == "mode") 1 else c(2, 3)
  col_names <- if (type == "mode") "mode" else c("mean", "variance")

  ch <- as.data.frame(data_object$ch)
  cs <- as.data.frame(data_object$cs)
  zh <- data_object$zh
  a <- data_object$a
  b <- data_object$b

  nh <- nrow(ch)
  est <- matrix(NA, nrow = nh, ncol = length(col_idx))

  for (i in seq_len(nh)) {
    data_obj <- bme_map(
      ch = ch[-i, ],
      cs = cs,
      zh = zh[-i],
      a = a,
      b = b
    )

    est[i, ] <- new_bme_estimate(
      x = ch[i, ], data_object = data_obj, lambda = lambda, mu = mu,
      nsmax = nsmax, nhmax = nhmax, n = n,
      zk_range = zk_range
    )[, col_idx]
  }

  ch_names <- if (is.null(colnames(ch))) c("coord.1", "coord.2") else colnames(ch)
  result <- cbind.data.frame(ch, zh, est, round(zh - est[, 1], 4), seq_len(nh))
  names(result) <- c(ch_names, "observed", col_names, "residual", "fold")

  return(structure(result, class = c("BMEmapping", "data.frame")))
}
