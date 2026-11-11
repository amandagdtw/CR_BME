# ============================================================================
# Wrapper to bme_mean and bme_mode functions
#
# Details:
# Compute the bme mean, variance and mode of an estimation location
#
# Inputs:
# -  x matrix of estimation locations
# -  ch matrix of hard data locations
# -  cs matrix of soft data locations
# -  zh vector of hard data
# -  a vector of lower bounds of soft data
# -  b vector of lower bounds of soft data
# -  model string name of covariance or variogram model
# -  nugget a non-negative value
# -  sill a non-negative value
# -  range a non-negative value
# -  nsmax number of soft data locations closer to the estimation location
# -  nhmax number of hard data locations closer to the estimation location
#
# Outputs:
# - A data frame of estimation locations with their corresponding bme mean,
#   variance and mode estimates
# ============================================================================
new_bme_estimate <- function(x, data_object, lambda, mu, nsmax = 5, nhmax = 5,
                             n = 50, zk_range = BMEmapping:::extended_range(data_object)) {
  # x <- matrix(c(x), ncol = 2)
  nk <- nrow(x)

  # set up container for estimates: mean, variance, mode
  df <- matrix(NA, ncol = 3, nrow = nk)

  for (i in 1:nk) {
    d <- new_prob_zk(
      x = x[i, ], data_object = data_object, lambda = lambda, mu = mu,
      nsmax = nsmax, nhmax = nhmax, n = n, zk_range = zk_range
    )

    delta <- d[2, 1] - d[1, 1]

    # compute mean
    zk_mean <- sum(d[, 1] * d[, 2] * delta)

    # compute variance
    zk_var <- sum((d[, 1] - zk_mean)^2 * d[, 2] * delta)

    # compute mode
    zk_mode <- d[which.max(d[, 2]), 1]

    # gather estimates
    df[i, ] <- round(c(zk_mode, zk_mean, zk_var), 4)
  }

  return(df)
}
