start_time <- Sys.time()

library(gstat)
library(BMEmapping)
library(sf)
library(caret)  # Necessary for createFolds()
library(dplyr)
library(evd)
library(mvtnorm)

#-------------------------------------------------------------------#
# Source CRBME Functions
#-------------------------------------------------------------------#
source("R/new_vg_results.R")
source("R/new_prob_zk.R")
source("R/new_bme_estimate.R")
source("R/new_bme_predict.R")
source("R/new_bme_cv.R")

#-------------------------------------------------------------------#
# 1. Setup Simulation Space & Geostatistical Background
#-------------------------------------------------------------------#
set.seed(1)
n <- 625

# Generate static spatial coordinates
coords <- data.frame(
  x = runif(n, 0, 100),
  y = runif(n, 0, 100)
)
grid_sf <- st_as_sf(coords, coords = c("x", "y"))

vgm_model <- vgm(psill = 0.5, model = "Exp", range = 60, nugget = 0.5)

g.dummy <- gstat(
  formula = z ~ 1,
  dummy = TRUE,
  beta = 0,
  model = vgm_model,
  nmax = 50
)

sim_count <- 2
sim <- predict(g.dummy, newdata = grid_sf, nsim = sim_count)

# Data partitioning indices
nh <- 1:125
ns <- 126:525
nk <- 526:625

#-------------------------------------------------------------------#
# 2. Helper Functions
#-------------------------------------------------------------------#
metrics <- function(truth, pred) {
  res <- truth - pred
  c(bias = mean(res), MAE = mean(abs(res)), MSE = mean(res^2))
}

summarize_results <- function(results_list) {
  avg_matrix <- Reduce("+", lapply(results_list, function(mat) {
    matrix(as.numeric(mat), ncol = 3)
  })) / length(results_list)

  df <- data.frame(
    method    = c("CBME", "", "CRBME", "", "variance"),
    pred_type = c("mode", "mean", "mode", "mean", ""),
    rbind(
      round(avg_matrix[1:4, ], 4),
      round(avg_matrix[5, , drop = FALSE], 2)
    )
  )
  colnames(df) <- c("method", "pred_type", "bias", "MAE", "MSE")
  df[is.na(df)] <- ""
  return(df)
}

#-------------------------------------------------------------------#
# 3. Main Simulation Loop Framework
#-------------------------------------------------------------------#
q_list <- list(c(0.05, 0.25), c(0.40, 0.60), c(0.75, 0.95))
l_values <- c(1, 3, 5)
results_all <- list()

# Define grid search space globally to reduce initialization overhead
param_grid <- data.frame(
  lambda = seq(0.1, 2, length.out = 10),
  mu     = seq(0.1, 2, length.out = 10)
)

mu <- 0.5
beta <- 1

xm <- 0.5
alpha <- 5


for (l in l_values) {
  results_l <- list()

  for (q in q_list) {
    results_q <- list()

    for (j in seq_len(sim_count)) {

      # Gaussian to Student-t
      #z_t <- qt(
      #  pnorm(sim[[paste0("sim", j)]]),
      #  df = 35
      #)

      # Transform to Gumbel field
      #z_t <- mu - beta * log(-log(pnorm(sim[[paste0("sim", j)]])))
      #t_data <- cbind(coords, z = z_t)

      # Transform Gaussian to Pareto Type I
      z_p <- xm * (1 - pnorm(sim[[paste0("sim", j)]]))^(-1 / alpha)

      u <- pnorm(sim[[paste0("sim", j)]])

      m <- mean(u)
      s <- sd(u)

      sigma_log <- sqrt(log(1 + (s / m)^2))
      mu_log <- log(m) - 0.5 * sigma_log^2

      z_l <- qlnorm(
        u,
        meanlog = mu_log,
        sdlog = sigma_log
      )


      t_data <- cbind(coords, z = z_l)

      # Split spatial subsets
      ch <- t_data[nh, c("x", "y")]; zh <- t_data[nh, "z"]
      cs <- t_data[ns, c("x", "y")]; zs <- t_data[ns, "z"]
      ck <- t_data[nk, c("x", "y")]; zk <- t_data[nk, "z"]

      # Fit variogram on available Hard Data
      df_vg        <- st_as_sf(t_data[nh, ], coords = c("x", "y"))
      vg           <- variogram(z ~ 1, data = df_vg)
      vg_model_fit <- fit.variogram(vg, model = vgm(c("Exp", "Sph")))

      model_type <- as.character(vg_model_fit[2, 1])
      nugget     <- vg_model_fit[1, 2]
      sill       <- vg_model_fit[2, 2]
      range      <- vg_model_fit[2, 3]

      #-------------------------------------------------------------#
      # Generate Realistic Soft-Interval Windows
      #-------------------------------------------------------------#
      n_soft <- length(zs)

      # Step A: Determine total length of individual interval window
      ln <- runif(n_soft, min = l - 0.5, max = l + 0.5)

      # Step B: Pick random relative positioning factor inside the window
      u <- runif(n_soft, min = q[1], max = q[2])

      # Step C: Hang the interval boundary symmetrically or skewed based on 'u'
      a <- zs - ln * u
      b <- zs + ln * (1 - u)

      # Capture descriptive statistics
      vs   <- var(zs)
      m_ln <- mean(ln)

      # Map full data inputs for final validation evaluation
      data_object <- bme_map(ch, cs, zh, a, b)

      # =========================================================
      # 5-FOLD CROSS-VALIDATION VIA PARAMETER GRID SEARCH
      # =========================================================
      folds <- createFolds(zh, k = 5)
      cv_summary <- vector("list", nrow(param_grid))

      for (g in seq_len(nrow(param_grid))) {
        current_lambda <- param_grid$lambda[g]
        current_mu     <- param_grid$mu[g]

        fold_results <- list()

        for (i in seq_along(folds)) {
          test_idx <- folds[[i]]

          ch_train <- ch[-test_idx, ]
          zh_train <- zh[-test_idx]
          ch_test  <- ch[test_idx, ]
          zh_test  <- zh[test_idx]

          # Construct training BME object isolated from test elements
          data_object_cv <- bme_map(
            ch = ch_train,
            cs = cs,
            zh = zh_train,
            a  = a,
            b  = b
          )

          # CRBME Prediction Execution
          preds <- new_bme_predict(
            x      = ch_test,
            data   = data_object_cv,
            lambda = current_lambda,
            mu     = current_mu,
            type   = "mean"
          )

          fold_results[[i]] <- data.frame(
            observed  = zh_test,
            predicted = preds$mean,
            residual  = zh_test - preds$mean
          )
        }

        # Combine metrics across all 5 folds
        combined_folds <- do.call(rbind, fold_results)

        ME   <- mean(combined_folds$residual)
        MAE  <- mean(abs(combined_folds$residual))
        RMSE <- sqrt(mean(combined_folds$residual^2))

        # Protected R-squared calculation to avoid division by zero variances
        r_sq_calc <- cor(combined_folds$observed, combined_folds$predicted)^2
        Rsquared  <- ifelse(is.na(r_sq_calc), 0, r_sq_calc)

        cv_summary[[g]] <- data.frame(
          lambda   = current_lambda,
          mu       = current_mu,
          ME       = round(ME, 4),
          MAE      = round(MAE, 4),
          RMSE     = round(RMSE, 4),
          Rsquared = round(Rsquared, 4)
        )
      }

      # Flatten list and find optimal parameters
      cv_summary_df <- do.call(rbind, cv_summary)
      best_row      <- cv_summary_df[which.min(cv_summary_df$RMSE), ]

      best_lambda <- best_row$lambda
      best_mu     <- best_row$mu

      #-------------------------------------------------------------#
      # Final Execution & Calculations (CBME / Tuned CRBME)
      #-------------------------------------------------------------#
      zk_cbme  <- BMEmapping:::bme_estimate(ck, data_object, model_type, nugget, sill,
                               range, zk_range = range(zh, a, b))[, 1:2]
      zk_crbme <- new_bme_estimate(ck, data_object, zk_range = range(zh, a, b),
                                   lambda = best_lambda, mu = best_mu)[, 1:2]

      # Metric Compilation
      DF_cbme  <- t(apply(zk_cbme, 2, function(p) metrics(zk, p)))
      DF_crbme <- t(apply(zk_crbme, 2, function(p) metrics(zk, p)))

      DF <- rbind(DF_cbme, DF_crbme, c(0, vs, m_ln))
      results_q[[j]] <- DF
    }

    q_label <- paste0("q", paste(q, collapse = "_"))
    results_l[[q_label]] <- summarize_results(results_q)
  }

  results_all[[paste0("l", l)]] <- results_l
}

#-------------------------------------------------------------------#
# 4. View Compiled Results Structure
#-------------------------------------------------------------------#
print(results_all)

#-------------------------------------------------------------------#
# 5. Save All Simulation Results as CSV
#-------------------------------------------------------------------#

results_df <- do.call(rbind, lapply(names(results_all), function(l_name) {

  do.call(rbind, lapply(names(results_all[[l_name]]), function(q_name) {

    df <- results_all[[l_name]][[q_name]]

    data.frame(
      l = l_name,
      q = q_name,
      df,
      row.names = NULL
    )
  }))
}))

# View combined results
print(results_df)


# Save as CSV
write.csv(results_df, "simulation_results.csv", row.names = FALSE)

end_time <- Sys.time()

runtime_seconds <- as.numeric(
  difftime(end_time, start_time, units = "secs")
)

cat("\nSimulation completed!\n")
cat("Total runtime:", round(runtime_seconds, 2), "seconds\n")
cat("Total runtime:", round(runtime_seconds / 60, 2), "minutes\n")
cat("Total runtime:", round(runtime_seconds / 3600, 3), "hours\n")
