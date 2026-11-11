#' @title Posterior Density Estimation at a Single Location
#'
#' @usage new_prob_zk(x, data_object, nsmax = 5, nhmax = 5, n = 50, lambda,
#'                  zk_range = BMEmapping:::extended_range(data_object))
#'
#' @param x A two-column matrix of spatial coordinates for a single estimation
#'        location.
#' @param data_object A list containing the hard and soft data.
#' @param nq A positive numeric value for the number of quantile levels
#'        (default is 3).
#' @param nsmax An integer specifying the maximum number of nearby soft data
#'        points to include for estimation (default is 5).
#' @param nhmax An integer specifying the maximum number of nearby hard data
#'        points to include for estimation (default is 5).
#' @param lambda An integer indicating the number of points at which
#' @param mu An integer indicating the number of points at which to evaluate the
#' @param n An integer indicating the number of points at which to evaluate the
#'        posterior density over \code{zk_range} (default is 50).
#' @param zk_range A numeric vector specifying the range over which to evaluate
#'        the unobserved value at the estimation location (\code{zk}). Although
#'        \code{zk} is unknown,  it is assumed to lie within a range similar to
#'        the observed data (\code{zh}, \code{a}, and \code{b}). It is advisable
#'        to explore the posterior distribution at a few locations using
#'        \code{prob_zk()} before finalizing this range
#'
#' @return A data frame with two columns: \code{zk_i} (assumed zk values) and
#'         \code{prob_zk_i} (corresponding posterior densities).
#'
#' @description
#' Computes the posterior and plots probability density function (PDF) at a
#' single unobserved spatial location using the Bayesian Maximum Entropy (BME)
#' framework. This function integrates both hard data (precise measurements) and
#' soft data (interval or uncertain observations), together with a specified
#' variogram model, to numerically estimate the posterior density across a
#' range of possible values.
#'
#' @examples
#' data("utsnowload")
#' x <- utsnowload[1, c("latitude", "longitude")]
#' ch <- utsnowload[2:67, c("latitude", "longitude")]
#' cs <- utsnowload[68:232, c("latitude", "longitude")]
#' zh <- utsnowload[2:67, "hard"]
#' a <- utsnowload[68:232, "lower"]
#' b <- utsnowload[68:232, "upper"]
#' data_object <- bme_map(ch, cs, zh, a, b)
#' new_prob_zk(x, data_object, lambda = 1)
#'
#' @export

new_prob_zk <- function(x, data_object, nsmax = 5, nhmax = 5, lambda, mu,
                        n = 50, zk_range = BMEmapping:::extended_range(data_object)) {

  x <- BMEmapping:::clean_input(x)
  BMEmapping:::check_xx(x)

  ch <- data_object$ch
  cs <- data_object$cs
  zh <- data_object$zh
  a <- data_object$a
  b <- data_object$b

  # center and radius
  zc = c(zh, (a + b)/2)
  zr = c(rep(0, length(zh)), b - a)

  d = cbind(zc, zr)

  if (nrow(x) != 1) {
    stop("Can only compute the mapping set for a single location")
  }

  # set up zk vector
  zk_vec <- seq(from = zk_range[1], to = zk_range[2], length.out = n)

  # Fit a variogram to each column of the data and extract the variogram
  # information (model and parameters)
  vg_info <- new_vg_results(ch = ch, cs = cs, zh = zh, a = a, b = b)

  # sorting nhmax hard data locations closest to the estimation location
  new_ch <- BMEmapping:::ch_nhmax(x = x, ch = ch, nhmax = nhmax)
  ch <- new_ch[[1]]
  index_h <- new_ch[[2]]
  zh <- zh[index_h]

  # sorting nsmax hard data locations closest to the estimation location
  new_cs <- BMEmapping:::cs_nsmax(x = x, cs = cs, nsmax = nsmax)
  cs <- new_cs[[1]]
  index_s <- new_cs[[2]]
  a <- a[index_s]
  b <- b[index_s]

  # Use the variogram information to compute the covariance matrices for the
  # hard, soft and estimations and average them
  #average_matrix <- BMEmapping:::covmat_avg(ch = ch, cs = cs, x = x, vg_info = vg_info)

  # variogram information for center
  model_c = vg_info$model[1]
  nugget_c = vg_info$nugget[1]
  sill_c = vg_info$sill[1]
  range_c = vg_info$range[1] * 100

  # variogram information for radius
  model_r = vg_info$model[2]
  nugget_r = vg_info$nugget[2]
  sill_r = vg_info$sill[2]
  range_r = vg_info$range[2] *100

  zs <- (a + b)/2

  #mu =  var(zc) / var(zr) #

  # Build the covariance matrices
  cov_h_h <- (lambda  * BMEmapping:::covmat(ch, ch, model_c, nugget_c, sill_c, range_c)) +
    (mu * BMEmapping:::covmat(ch, ch, model_r, nugget_r, sill_r, range_r))

  cov_h_k <- (lambda  * BMEmapping:::covmat(ch, x,  model_c, nugget_c, sill_c, range_c)) +
    (mu * BMEmapping:::covmat(ch, x,  model_r, nugget_r, sill_r, range_r))

  cov_h_s <- (lambda  * BMEmapping:::covmat(ch, cs, model_c, nugget_c, sill_c, range_c)) +
    (mu * BMEmapping:::covmat(ch, cs, model_r, nugget_r, sill_r, range_r))

  cov_k_h <- (lambda  * BMEmapping:::covmat(x,  ch, model_c, nugget_c, sill_c, range_c)) +
    (mu * BMEmapping:::covmat(x,  ch, model_r, nugget_r, sill_r, range_r))

  cov_k_k <- (lambda  * BMEmapping:::covmat(x,  x,  model_c, nugget_c, sill_c, range_c)) +
    (mu * BMEmapping:::covmat(x,  x,  model_r, nugget_r, sill_r, range_r))

  cov_s_h <- (lambda  * BMEmapping:::covmat(cs, ch, model_c, nugget_c, sill_c, range_c)) +
    (mu * BMEmapping:::covmat(cs, ch, model_r, nugget_r, sill_r, range_r))

  cov_s_k <- (lambda  * BMEmapping:::covmat(cs, x, model_c, nugget_c, sill_c, range_c)) +
    (mu * BMEmapping:::covmat(cs, x, model_r, nugget_r, sill_r, range_r))

  cov_s_s <- (lambda  * BMEmapping:::covmat(cs, cs, model_c, nugget_c, sill_c, range_c)) +
    (mu * BMEmapping:::covmat(cs, cs, model_r, nugget_r, sill_r, range_r))

  cov_k_s <- t(cov_s_k)
  cov_s_h <- t(cov_h_s)
  cov_k_h <- t(cov_h_k)

  # Composite covariances
  cov_kh_kh <- rbind(cbind(cov_k_k, cov_k_h), cbind(cov_h_k, cov_h_h))
  cov_s_kh <- cbind(cov_s_k, cov_s_h)
  cov_kh_s <- t(cov_s_kh)

  ###########################################################################
  #   Part a: compute normalization constant
  ###########################################################################

  # lower and upper limits
  lower_a <- a
  upper_a <- b

  inv_cov_hs_hs <- solve(cov_h_h)

  # covariance matrix
  cov_a <- cov_s_s - cov_s_h %*% inv_cov_hs_hs %*% cov_h_s
  if (det(cov_a) <= 0) cov_a <- cov_s_s

  # mean vector
  mu_a <- c(cov_s_h %*% inv_cov_hs_hs %*% zh)

  aa <- mvtnorm::pmvnorm(
    lower = lower_a, upper = upper_a, mean = mu_a,
    sigma = cov_a
  )[1]

  # set up container
  pk <- numeric()

  # Part B: conditional mean and covariance of zk
  m_k <- c(cov_k_h %*% solve(cov_h_h, zh)) # mean
  cov_k <- cov_k_k - cov_k_h %*% inv_cov_hs_hs %*% cov_h_k # covariance

  # Part C: compute integral of soft data
  # conditional variance
  inv_cov_kh_kh <- solve(cov_kh_kh)
  cov_soft <- cov_s_s - cov_s_kh %*% inv_cov_kh_kh %*% cov_kh_s
  if (det(cov_soft) <= 0) cov_soft <- cov_s_s

  for (i in 1:n) {
    ###########################################################################
    #      zk, zkh values
    ###########################################################################

    zk <- zk_vec[i]
    zk_h <- c(zk, zh)

    ###########################################################################
    #   Part B: compute density of zk
    ###########################################################################

    # density
    f_zk <- mvtnorm::dmvnorm(x = zk, mean = m_k, sigma = cov_k)[1]

    ###########################################################################
    #   Part C: compute integral of soft data
    ###########################################################################

    # lower and upper limits
    lower_soft <- c(a - cov_s_kh %*% inv_cov_kh_kh %*% zk_h)
    upper_soft <- c(b - cov_s_kh %*% inv_cov_kh_kh %*% zk_h)

    # conditional mean
    m_soft <- rep(0, nsmax)

    # Compute multidimensional integral
    f_soft <- mvtnorm::pmvnorm(
      lower = lower_soft, upper = upper_soft, mean = m_soft,
      sigma = cov_soft
    )[1]

    if (f_soft == 0) f_soft <- 1e-4
    if (aa == 0) aa <- 1e-4

    pk[i] <- round(((1 / aa) * f_zk * f_soft), 5)
  }

  d <- data.frame("zk_i" = zk_vec, "prob_zk_i" = pk)
  df <- d[!rowSums(is.na(d)), ]

  return(structure(df, class = c("BMEmapping", "data.frame")))
}
