# ============================================================================
# Wrapper to new prob_zk function
#
# Details:
# Fits a variogram to each column of the data (hard and qunatile data) and
# extract the model parameters
#
# Inputs:
# -  combine data: hard and qunatile data [[zh,a) (zh,q1) ... (zh,qn) (zh,b)]
# -  ch hard data locations
# -  cs soft data locations
#
# Outputs:
# - A (qn + 2) by 3 data frame of variogram model and prameters. Each row
#   represents the variogram model and prameters of each column.
# ============================================================================
#' @importFrom stats as.formula
#'
new_vg_results <- function(ch, cs, zh, a, b) {
  # --- Combine coordinates ---
  coords <- rbind(as.matrix(ch), as.matrix(cs))
  coords <- as.data.frame(coords)
  colnames(coords) <- c("x", "y")

  # center and radius
  zc = c(zh, (a + b)/2)
  zr = c(rep(0, length(zh)), b - a)

  d = cbind(zc, zr)

  # --- Combine with data ---
  d_combine <- cbind(coords, d)
  colnames(d_combine) <- c("x", "y", paste0("var", seq_len(ncol(d))))

  # --- Convert to sf ---
  sf_data <- sf::st_as_sf(d_combine, coords = c("x", "y"))

  results <- list()

  for (i in seq_len(ncol(d))) {
    varname <- paste0("var", i)
    form <- stats::as.formula(paste(varname, "~ 1"))

    # --- Empirical variogram ---
    vgm_emp <- gstat::variogram(form, data = sf_data)

    # --- Fit model ---
    vgm_fit <- gstat::fit.variogram(
      vgm_emp,
      model = gstat::vgm(c("Exp", "Sph"))
    )

    # --- Extract parameters ---
    if (nrow(vgm_fit) == 2) {
      model  <- as.character(vgm_fit$model[2])
      nugget <- vgm_fit$psill[1]
      sill   <- vgm_fit$psill[2]
      range  <- vgm_fit$range[2]
    } else {
      model  <- as.character(vgm_fit$model[1])
      nugget <- 0
      sill   <- vgm_fit$psill[1]
      range  <- vgm_fit$range[1]
    }

    results[[i]] <- data.frame(
      model = model,
      nugget = nugget,
      sill = sill,
      range = range * 100
    )
  }

  # --- Combine results ---
  df <- do.call(rbind, results)
  rownames(df) <- NULL

  return(df)
}
