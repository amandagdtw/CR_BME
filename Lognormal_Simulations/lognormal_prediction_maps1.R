#===================================================================#
# FIGURE 3: LOGNORMAL SPATIAL PREDICTION MAPS
# Interpolated reference vs CBME vs CRBME (mean predictions)
# Scenario with largest mean-MSE improvement: l=3, q=(0.40,0.60)
#===================================================================#
library(devtools)
library(gstat)
library(BMEmapping)
library(sf)
library(caret)
library(ggplot2)
library(viridis)

# Load your local CRBME functions
source("R/new_vg_results.R")
source("R/new_prob_zk.R")
source("R/new_bme_estimate.R")
source("R/new_bme_predict.R")

#-------------------------------------------------------------------#
# 1. SELECT CASE AND GENERATE THE ORIGINAL GAUSSIAN SIMULATIONS
#-------------------------------------------------------------------#
fig_l <- 3
fig_q <- c(0.40, 0.60)
fig_j <- 1
n_grid <- 30

set.seed(1)
n <- 625
coords <- data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100))
grid_sf <- st_as_sf(coords, coords = c("x", "y"))
vgm_model <- vgm(psill = 0.5, model = "Exp", range = 60, nugget = 0.5)
g.dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0,
                 model = vgm_model, nmax = 50)
sim_count <- 2
sim <- predict(g.dummy, newdata = grid_sf, nsim = sim_count)

nh <- 1:125
ns <- 126:525
nk <- 526:625
param_grid <- data.frame(
  lambda = seq(0.1, 2, length.out = 10),
  mu = seq(0.1, 2, length.out = 10)
)

#-------------------------------------------------------------------#
# 2. GAUSSIAN TO LOGNORMAL: EXACT TRANSFORMATION FROM ORIGINAL STUDY
#-------------------------------------------------------------------#
gaussian_f <- sim[[paste0("sim", fig_j)]]
u_prob <- pnorm(gaussian_f)
m <- mean(u_prob)
s <- sd(u_prob)
sigma_log <- sqrt(log(1 + (s / m)^2))
mu_log <- log(m) - 0.5 * sigma_log^2
z_l_f <- qlnorm(u_prob, meanlog = mu_log, sdlog = sigma_log)

coords_f <- as.data.frame(st_coordinates(sim))
names(coords_f) <- c("x", "y")
t_data <- data.frame(coords_f, z = z_l_f)

#-------------------------------------------------------------------#
# 3. HARD, SOFT, AND VALIDATION DATA
#-------------------------------------------------------------------#
ch_f <- t_data[nh, c("x", "y"), drop = FALSE]
zh_f <- t_data[nh, "z"]
cs_f <- t_data[ns, c("x", "y"), drop = FALSE]
zs_f <- t_data[ns, "z"]
ck_f <- t_data[nk, c("x", "y"), drop = FALSE]
zk_f <- t_data[nk, "z"]

#-------------------------------------------------------------------#
# 4. FIT VARIOGRAM ON HARD DATA
#-------------------------------------------------------------------#
df_vg_f <- st_as_sf(t_data[nh, ], coords = c("x", "y"))
vg_emp_f <- variogram(z ~ 1, data = df_vg_f)
vgf <- fit.variogram(vg_emp_f, model = vgm(c("Exp", "Sph")))
print(vgf)
if (nrow(vgf) < 2) stop("Variogram fit has fewer than two components")
model_f <- as.character(vgf$model[2])
nugget_f <- vgf$psill[1]
sill_f <- vgf$psill[2]
range_f <- vgf$range[2]

#-------------------------------------------------------------------#
# 5. SOFT INTERVALS: SELECTED l AND q
#-------------------------------------------------------------------#
# This seed reproduces this map's intervals, NOT necessarily the
# intervals drawn at this point in the original nested simulation.
set.seed(1234)
n_soft_f <- length(zs_f)
ln_f <- runif(n_soft_f, fig_l - 0.5, fig_l + 0.5)
u_interval <- runif(n_soft_f, fig_q[1], fig_q[2])
a_f <- zs_f - ln_f * u_interval
b_f <- zs_f + ln_f * (1 - u_interval)
data_object_f <- bme_map(ch = ch_f, cs = cs_f, zh = zh_f,
                         a = a_f, b = b_f)
zk_rng_f <- base::range(zh_f, a_f, b_f)

#-------------------------------------------------------------------#
# 6. FIVE-FOLD CV FOR CRBME, USING ORIGINAL PARAMETER GRID
#-------------------------------------------------------------------#
set.seed(2026)
folds_f <- createFolds(zh_f, k = 5)
cv_f <- vector("list", nrow(param_grid))
for (g in seq_len(nrow(param_grid))) {
  fold_results <- vector("list", length(folds_f))
  for (i in seq_along(folds_f)) {
    test_idx <- folds_f[[i]]
    data_object_cv <- bme_map(
      ch = ch_f[-test_idx, , drop = FALSE],
      cs = cs_f,
      zh = zh_f[-test_idx],
      a = a_f,
      b = b_f
    )
    preds <- new_bme_predict(
      x = ch_f[test_idx, , drop = FALSE],
      data = data_object_cv,
      lambda = param_grid$lambda[g],
      mu = param_grid$mu[g],
      type = "mean"
    )
    fold_results[[i]] <- data.frame(
      observed = zh_f[test_idx], predicted = preds$mean
    )
  }
  combined_folds <- do.call(rbind, fold_results)
  RMSE <- sqrt(mean((combined_folds$observed -
                     combined_folds$predicted)^2))
  cv_f[[g]] <- data.frame(lambda = param_grid$lambda[g],
                          mu = param_grid$mu[g], RMSE = RMSE)
  message("CV combination ", g, "/", nrow(param_grid),
          ": RMSE = ", signif(RMSE, 5))
}
cv_f_df <- do.call(rbind, cv_f)
if (!any(is.finite(cv_f_df$RMSE))) stop("No finite CV RMSE values")
best_f <- cv_f_df[which.min(replace(cv_f_df$RMSE,
                                   !is.finite(cv_f_df$RMSE), Inf)), , drop = FALSE]
lambda_f <- best_f$lambda
mu_f <- best_f$mu
cat("\nSelected Lognormal case: l =", fig_l,
    ", q =", paste(fig_q, collapse = ", "),
    ", realization =", fig_j, "\n")
cat("Best lambda =", lambda_f, "; best mu =", mu_f,
    "; CV RMSE =", best_f$RMSE, "\n")

#-------------------------------------------------------------------#
# 7. REGULAR 30 x 30 PREDICTION GRID
#-------------------------------------------------------------------#
grid_xy <- expand.grid(x = seq(0, 100, length.out = n_grid),
                       y = seq(0, 100, length.out = n_grid))

#-------------------------------------------------------------------#
# 8. KRIGED REFERENCE FIELD (NOT EXACT SIMULATED TRUTH)
#-------------------------------------------------------------------#
truth_sf <- st_as_sf(t_data, coords = c("x", "y"))
grid_sf3 <- st_as_sf(grid_xy, coords = c("x", "y"))
vg_truth <- fit.variogram(variogram(z ~ 1, data = truth_sf),
                          model = vgm(c("Exp", "Sph")))
truth_grid <- krige(z ~ 1, locations = truth_sf,
                    newdata = grid_sf3, model = vg_truth)$var1.pred

#-------------------------------------------------------------------#
# 9. CBME AND CRBME MEAN PREDICTIONS ON THE SAME GRID
#-------------------------------------------------------------------#
message("Predicting at ", nrow(grid_xy), " grid locations...")
cbme_output <- BMEmapping:::bme_estimate(grid_xy, data_object_f,
                            model_f, nugget_f, sill_f, range_f,
                            zk_range = zk_rng_f)
crbme_output <- new_bme_estimate(grid_xy, data_object_f,
                                  zk_range = zk_rng_f,
                                  lambda = lambda_f, mu = mu_f)
pred_cbme <- cbme_output[, 2]
pred_crbme <- crbme_output[, 2]
stopifnot(length(truth_grid) == nrow(grid_xy),
          length(pred_cbme) == nrow(grid_xy),
          length(pred_crbme) == nrow(grid_xy))

#-------------------------------------------------------------------#
# 10. OPTIONAL HELD-OUT VALIDATION FOR THIS SINGLE REALIZATION
#-------------------------------------------------------------------#
cbme_valid <- BMEmapping:::bme_estimate(ck_f, data_object_f,
                           model_f, nugget_f, sill_f, range_f,
                           zk_range = zk_rng_f)[, 2]
crbme_valid <- new_bme_estimate(ck_f, data_object_f,
                                zk_range = zk_rng_f,
                                lambda = lambda_f, mu = mu_f)[, 2]
validation_metrics <- data.frame(
  method = c("CBME", "CRBME"),
  MAE = c(mean(abs(zk_f - cbme_valid)),
          mean(abs(zk_f - crbme_valid))),
  MSE = c(mean((zk_f - cbme_valid)^2),
          mean((zk_f - crbme_valid)^2))
)
print(validation_metrics)

#-------------------------------------------------------------------#
# 11. MAP DATA AND COMMON COLOR SCALE
#-------------------------------------------------------------------#
map_df <- rbind(
  data.frame(grid_xy, value = truth_grid, panel = "Reference field"),
  data.frame(grid_xy, value = pred_cbme, panel = "CBME"),
  data.frame(grid_xy, value = pred_crbme, panel = "CRBME")
)
map_df$panel <- factor(map_df$panel,
                       levels = c("Reference field", "CBME", "CRBME"))
lim <- base::range(map_df$value, na.rm = TRUE)

#-------------------------------------------------------------------#
# 12. PROFESSIONAL ECOSTA PLOT
#-------------------------------------------------------------------#
theme_ecosta <- theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.6),
        axis.text = element_text(colour = "black", size = 11),
        axis.title = element_text(size = 13),
        strip.background = element_rect(fill = "grey95", colour = "black"),
        strip.text = element_text(face = "bold", size = 13),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        panel.spacing = grid::unit(0.5, "cm"))

p3 <- ggplot(map_df, aes(x = x, y = y, fill = value)) +
  geom_raster(interpolate = TRUE) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_viridis_c(option = "D", limits = lim, name = "z") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
  labs(x = "x", y = "y") + theme_ecosta
print(p3)

#-------------------------------------------------------------------#
# 13. SAVE FIGURES AND DATA
#-------------------------------------------------------------------#
dir.create("figs", showWarnings = FALSE)
ggsave("figs/lognormal_prediction_maps.pdf", plot = p3,
       width = 11, height = 4.2)
ggsave("figs/lognormal_prediction_maps.png", plot = p3,
       width = 11, height = 4.2, dpi = 400, bg = "white")
saveRDS(list(map_df = map_df, lambda = lambda_f, mu = mu_f,
             fig_l = fig_l, fig_q = fig_q, fig_j = fig_j,
             cv_results = cv_f_df, validation_metrics = validation_metrics,
             lognormal_parameters = c(meanlog = mu_log, sdlog = sigma_log),
             intervals = data.frame(a = a_f, b = b_f)),
        file = "figs/lognormal_prediction_map_results.rds")
cat("\nLognormal prediction maps saved successfully.\n")
