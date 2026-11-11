# ================================================================
# CR-BME: Gumbel simulation + figures from ORIGINAL realizations
# EcoSta 2026; selected setting l=3, q=(0.75,0.95)
#
# Re-runs the ORIGINAL nine settings and five realizations in their
# original order. No extra random draws occur inside the loops.
# Uses 625 original simulated values for interpolated truth display.
# Evaluates accuracy ONLY at the 100 original hold-out sites.
# ================================================================

library(devtools)
library(gstat)
library(BMEmapping)
library(sf)
library(ggplot2)
library(dplyr)
library(caret)
load_all()
dir.create('figs', showWarnings = FALSE)

# ------------------------- Original setup ------------------------
set.seed(1234)
n <- 625
gum_mu <- 0.5
gum_beta <- 1
coords <- data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100))
grid_sf <- st_as_sf(coords, coords = c('x', 'y'))
nh <- 1:125
ns <- 126:525
nk <- 526:625
vgm_model <- vgm(psill = 0.5, model = 'Exp', range = 60, nugget = 0.5)
g.dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0,
                 model = vgm_model, nmax = 50)
sim_count <- 5
sim <- predict(g.dummy, newdata = grid_sf, nsim = sim_count)
to_gumbel <- function(z) gum_mu - gum_beta * log(-log(pnorm(z)))

metrics <- function(truth, pred) {
  res <- truth - pred
  c(bias = mean(res), MAE = mean(abs(res)), MSE = mean(res^2))
}
summarize_results <- function(results_list) {
  avg_matrix <- Reduce('+', lapply(results_list, function(mat)
    matrix(as.numeric(mat), ncol = 3))) / length(results_list)
  df <- data.frame(method = c('CBME', '', 'CRBME', '', 'variance'),
                   pred_type = c('mode', 'mean', 'mode', 'mean', ''),
                   rbind(round(avg_matrix[1:4, ], 4),
                         round(avg_matrix[5, , drop = FALSE], 2)))
  colnames(df) <- c('method', 'pred_type', 'bias', 'MAE', 'MSE')
  df[is.na(df)] <- ''
  df
}

# ----------------- Original full simulation loop -----------------
q_list <- list(c(0.05, 0.25), c(0.40, 0.60), c(0.75, 0.95))
l_values <- c(1, 3, 5)
param_grid <- data.frame(lambda = seq(0.5, 5, length.out = 10),
                         mu = seq(1, 5, length.out = 10))
results_all <- list()
saved_reps <- vector('list', sim_count)

for (l in l_values) {
  results_l <- list()
  for (q in q_list) {
    results_q <- list()
    for (j in seq_len(sim_count)) {
      # Exact same Gumbel realization as main simulation
      z_t <- to_gumbel(sim[[paste0('sim', j)]])
      t_data <- cbind(coords, z = z_t)
      ch <- t_data[nh, c('x', 'y')]; zh <- t_data[nh, 'z']
      cs <- t_data[ns, c('x', 'y')]; zs <- t_data[ns, 'z']
      ck <- t_data[nk, c('x', 'y')]; zk <- t_data[nk, 'z']

      # Keep original variogram-fitting procedure for reproducibility
      df_vg <- st_as_sf(t_data[nh, ], coords = c('x', 'y'))
      vg <- variogram(z ~ 1, data = df_vg)
      vg_model_fit <- fit.variogram(vg, model = vgm(c('Exp', 'Sph')))
      model_type <- as.character(vg_model_fit[2, 1])
      nugget <- vg_model_fit[1, 2]
      sill <- vg_model_fit[2, 2]
      range_par <- vg_model_fit[2, 3]

      # Exact soft intervals, retaining the original RNG sequence
      n_soft <- length(zs)
      ln <- runif(n_soft, min = l - 0.5, max = l + 0.5)
      u <- runif(n_soft, min = q[1], max = q[2])
      a <- zs - ln * u
      b <- zs + ln * (1 - u)
      vs <- var(zs)
      m_ln <- mean(ln)
      data_object <- bme_map(ch, cs, zh, a, b)

      # Original CV: 10 paired tuning settings and 5 folds
      folds <- createFolds(zh, k = 5)
      cv_summary <- vector('list', nrow(param_grid))
      for (g in seq_len(nrow(param_grid))) {
        current_lambda <- param_grid$lambda[g]
        current_mu <- param_grid$mu[g]
        fold_results <- list()
        for (i in seq_along(folds)) {
          test_idx <- folds[[i]]
          ch_train <- ch[-test_idx, , drop = FALSE]
          zh_train <- zh[-test_idx]
          ch_test <- ch[test_idx, , drop = FALSE]
          zh_test <- zh[test_idx]
          data_object_cv <- bme_map(ch = ch_train, cs = cs,
                                    zh = zh_train, a = a, b = b)
          preds <- new_bme_predict(x = ch_test, data = data_object_cv,
                                   lambda = current_lambda, mu = current_mu,
                                   type = 'mean')
          fold_results[[i]] <- data.frame(
            observed = zh_test, predicted = preds$mean,
            residual = zh_test - preds$mean)
        }
        combined_folds <- do.call(rbind, fold_results)
        ME <- mean(combined_folds$residual)
        MAE <- mean(abs(combined_folds$residual))
        RMSE <- sqrt(mean(combined_folds$residual^2))
        r_sq_calc <- cor(combined_folds$observed, combined_folds$predicted)^2
        Rsquared <- ifelse(is.na(r_sq_calc), 0, r_sq_calc)
        cv_summary[[g]] <- data.frame(
          lambda = current_lambda, mu = current_mu,
          ME = round(ME, 4), MAE = round(MAE, 4),
          RMSE = round(RMSE, 4), Rsquared = round(Rsquared, 4))
      }
      cv_summary_df <- do.call(rbind, cv_summary)
      best_row <- cv_summary_df[which.min(cv_summary_df$RMSE), ]
      best_lambda <- best_row$lambda
      best_mu <- best_row$mu

      # Original held-out predictions and summary
      zk_cbme <- bme_estimate(ck, data_object, model_type, nugget,
                              sill, range_par,
                              zk_range = base::range(zh, a, b))[, 1:2]
      zk_crbme <- new_bme_estimate(
        ck, data_object, zk_range = base::range(zh, a, b),
        lambda = best_lambda, mu = best_mu)[, 1:2]
      DF_cbme <- t(apply(zk_cbme, 2, function(p) metrics(zk, p)))
      DF_crbme <- t(apply(zk_crbme, 2, function(p) metrics(zk, p)))
      DF <- rbind(DF_cbme, DF_crbme, c(0, vs, m_ln))
      results_q[[j]] <- DF

      # Save selected case WITHOUT any extra random-number draws
      if (l == 3 && identical(q, c(0.75, 0.95))) {
        saved_reps[[j]] <- list(
          j = j, l = l, q = q, coords = coords, z_data = z_t,
          ch = ch, zh = zh, cs = cs, zs = zs, a = a, b = b,
          ck = ck, zk = zk, model_type = model_type,
          nugget = nugget, sill = sill, range_par = range_par,
          best_lambda = best_lambda, best_mu = best_mu,
          zk_cbme = zk_cbme, zk_crbme = zk_crbme,
          best_row = best_row, DF = DF)
        saveRDS(saved_reps[[j]],
                sprintf('figs/gumbel_l3_qhigh_rep%d.rds', j))
      }
    }
    q_label <- paste0('q', paste(q, collapse = '_'))
    results_l[[q_label]] <- summarize_results(results_q)
  }
  results_all[[paste0('l', l)]] <- results_l
}
print(results_all)
saveRDS(results_all, 'figs/gumbel_all_results.rds')

# ---------------------- Figure settings ---------------------------
# Choose which ORIGINAL realization to display (1 to 5)
rep_to_plot <- 1
case <- saved_reps[[rep_to_plot]]
stopifnot(!is.null(case))

theme_pres <- theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(colour = 'grey30', fill = NA),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = 'bold', size = 12),
        legend.position = 'right')

# ------------------------- Figure 1 -------------------------------
loc_df <- coords %>%
  mutate(type = factor(c(rep('Hard (125)', length(nh)),
                         rep('Soft (400)', length(ns)),
                         rep('Hold-out (100)', length(nk))),
                       levels = c('Hard (125)', 'Soft (400)',
                                  'Hold-out (100)')))
p1 <- ggplot(loc_df, aes(x, y, shape = type, colour = type)) +
  geom_point(size = 2.1, stroke = 0.8) +
  scale_shape_manual(values = c(16, 1, 4)) +
  scale_colour_manual(values = c('#004886', '#8AA7C4', '#DC7828')) +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = 'x', y = 'y', colour = NULL, shape = NULL) + theme_pres
ggsave('figs/locations.pdf', p1, width = 6.4, height = 5)

# ------------------------- Figure 2 -------------------------------
z_all <- do.call(rbind, lapply(seq_len(sim_count), function(j)
  data.frame(z = to_gumbel(sim[[paste0('sim', j)]]),
             sim = paste0('Realization ', j))))
zz <- seq(min(z_all$z), max(z_all$z), length.out = 400)
ref <- data.frame(z = zz,
                  d = exp(-(zz - gum_mu) / gum_beta -
                            exp(-(zz - gum_mu) / gum_beta)) / gum_beta)
p2 <- ggplot() +
  geom_histogram(data = subset(z_all, sim == 'Realization 1'),
                 aes(z, after_stat(density)), bins = 40,
                 fill = 'grey85', colour = 'white') +
  geom_density(data = z_all, aes(z, group = sim),
               colour = '#8AA7C4', linewidth = 0.5) +
  geom_line(data = ref, aes(z, d),
            colour = '#DC7828', linewidth = 1.2) +
  labs(x = 'z', y = 'density', title = 'Gumbel(0.5, 1) marginal') + theme_pres
ggsave('figs/marginal_density.pdf', p2, width = 6.4, height = 4.4)

# ------------------------- Figure 3 -------------------------------
# Original main simulation has no exact field values on grid nodes.
# Use IDW solely to DISPLAY interpolated truth from the 625 sites.
# Both methods predict on the grid using SAME original data/params.
n_grid <- 30
grid_xy <- expand.grid(x = seq(0, 100, length.out = n_grid),
                       y = seq(0, 100, length.out = n_grid))
data_object <- bme_map(case$ch, case$cs, case$zh, case$a, case$b)
zk_rng <- base::range(case$zh, case$a, case$b)

message('Predicting CBME at ', nrow(grid_xy), ' grid nodes')
pred_cbme <- bme_estimate(
  grid_xy, data_object, case$model_type, case$nugget,
  case$sill, case$range_par, zk_range = zk_rng)[, 2]
message('Predicting CR-BME at ', nrow(grid_xy), ' grid nodes')
pred_crbme <- new_bme_estimate(
  grid_xy, data_object, zk_range = zk_rng,
  lambda = case$best_lambda, mu = case$best_mu)[, 2]

truth_sf <- st_as_sf(data.frame(coords, z = case$z_data),
                     coords = c('x', 'y'))
grid_pred_sf <- st_as_sf(grid_xy, coords = c('x', 'y'))
truth_idw <- gstat::idw(z ~ 1, locations = truth_sf,
                        newdata = grid_pred_sf, idp = 2)$var1.pred

map_df <- rbind(
  data.frame(grid_xy, value = truth_idw,
             panel = 'Truth (interpolated)'),
  data.frame(grid_xy, value = pred_cbme, panel = 'CBME'),
  data.frame(grid_xy, value = pred_crbme, panel = 'CR-BME'))
map_df$panel <- factor(map_df$panel,
                       levels = c('Truth (interpolated)', 'CBME', 'CR-BME'))
p3 <- ggplot(map_df, aes(x, y, fill = value)) +
  geom_raster(interpolate = TRUE) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_viridis_c(option = 'D', limits = base::range(map_df$value),
                       name = 'z') +
  coord_fixed(expand = FALSE) +
  labs(title = 'Gumbel spatial prediction',
       subtitle = sprintf('l = 3, q = (0.75, 0.95), realization %d',
                          rep_to_plot), x = 'x', y = 'y') + theme_pres
ggsave('figs/prediction_maps.pdf', p3, width = 11, height = 4.5)

# ------------------ Hold-out error figures ------------------------
# Use actual predictions at 100 original hold-out locations.
# Do NOT compare grid predictions against the IDW truth surface.
cb <- case$zk_cbme[, 2]
cr <- case$zk_crbme[, 2]
zk <- case$zk
ck <- case$ck
err_df <- rbind(
  data.frame(ck, err = abs(cb - zk), panel = 'CBME'),
  data.frame(ck, err = abs(cr - zk), panel = 'CR-BME'))
err_df$panel <- factor(err_df$panel, levels = c('CBME', 'CR-BME'))
p4 <- ggplot(err_df, aes(x, y, colour = err)) +
  geom_point(size = 3) + facet_wrap(~panel, nrow = 1) +
  scale_colour_viridis_c(option = 'magma',
                         limits = base::range(err_df$err), name = '|error|') +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(title = 'Absolute errors at original hold-out sites',
       x = 'x', y = 'y') + theme_pres
ggsave('figs/holdout_error_maps.pdf', p4, width = 9, height = 4.5)

improvement_df <- data.frame(
  ck, improvement = abs(cb - zk) - abs(cr - zk))
a <- max(abs(improvement_df$improvement))
p5 <- ggplot(improvement_df, aes(x, y, colour = improvement)) +
  geom_point(size = 3) +
  scale_colour_gradient2(low = '#C65353', mid = 'white', high = '#168C72',
                         midpoint = 0, limits = c(-a, a),
                         name = 'Error reduction') +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(title = 'CR-BME improvement at hold-out sites',
       subtitle = 'Positive = CR-BME has lower absolute error',
       x = 'x', y = 'y') + theme_pres
ggsave('figs/improvement_map.pdf', p5, width = 6.5, height = 5)

# ----------------------- Accuracy summary -------------------------
comparison <- do.call(rbind, lapply(saved_reps, function(x) {
  cb <- x$zk_cbme[, 2]; cr <- x$zk_crbme[, 2]; z <- x$zk
  c(Realization = x$j,
    CBME_MAE = mean(abs(cb - z)), CRBME_MAE = mean(abs(cr - z)),
    CBME_MSE = mean((cb - z)^2), CRBME_MSE = mean((cr - z)^2),
    lambda = x$best_lambda, mu = x$best_mu)
}))
comparison <- as.data.frame(comparison)
comparison$MAE_reduction_pct <- with(comparison,
  100 * (CBME_MAE - CRBME_MAE) / CBME_MAE)
print(comparison)
cat('\nAverage across five original realizations:\n')
print(colMeans(comparison[, c('CBME_MAE', 'CRBME_MAE',
                              'CBME_MSE', 'CRBME_MSE')]))
write.csv(comparison, 'figs/gumbel_replication_metrics.csv',
          row.names = FALSE)
print(p1); print(p2); print(p3); print(p4); print(p5)
