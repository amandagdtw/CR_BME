
#===================================================================#
# FIGURE 3: SPATIAL PREDICTION MAPS
# Truth vs CBME vs CRBME
# Selected scenario: l = 3, q = (0.75, 0.95)
#===================================================================#

library(devtools)
library(gstat)
library(BMEmapping)
library(sf)
library(caret)
library(ggplot2)
library(viridis)

# Load local package functions if needed
# load_all()

#-------------------------------------------------------------------#
# 1. SELECT SIMULATION CASE
#-------------------------------------------------------------------#

fig_l  <- 3
fig_q  <- c(0.75, 0.95)
fig_j  <- 1
n_grid <- 30

# Use the same tuning grid as the original simulation
param_grid <- expand.grid(
  lambda = seq(0.5, 5, length.out = 10),
  mu     = seq(1, 5, length.out = 10)
)

# Use original partition indices
nh <- 1:125
ns <- 126:525
nk <- 526:625

#-------------------------------------------------------------------#
# 2. TRANSFORM GAUSSIAN FIELD TO GUMBEL
#-------------------------------------------------------------------#

to_gumbel <- function(x, location = 0.5, scale = 1) {

  p <- pnorm(x)

  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)

  location - scale * log(-log(p))
}

z_t_f <- to_gumbel(
  sim[[paste0("sim", fig_j)]]
)

# Extract original coordinates from sim
coords_f <- as.data.frame(st_coordinates(sim))
names(coords_f) <- c("x", "y")

t_data <- data.frame(
  coords_f,
  z = z_t_f
)

#-------------------------------------------------------------------#
# 3. SEPARATE HARD, SOFT AND VALIDATION DATA
#-------------------------------------------------------------------#

ch_f <- t_data[nh, c("x", "y")]
zh_f <- t_data[nh, "z"]

cs_f <- t_data[ns, c("x", "y")]
zs_f <- t_data[ns, "z"]

ck_f <- t_data[nk, c("x", "y")]
zk_f <- t_data[nk, "z"]

#-------------------------------------------------------------------#
# 4. FIT VARIOGRAM USING HARD DATA
#-------------------------------------------------------------------#

df_vg_f <- st_as_sf(
  t_data[nh, ],
  coords = c("x", "y")
)

vg_emp_f <- variogram(
  z ~ 1,
  data = df_vg_f
)

vgf <- fit.variogram(
  vg_emp_f,
  model = vgm(c("Exp", "Sph"))
)

model_f  <- as.character(vgf[2, 1])
nugget_f <- vgf[1, 2]
sill_f   <- vgf[2, 2]
range_f  <- vgf[2, 3]

#-------------------------------------------------------------------#
# 5. GENERATE SOFT INTERVALS
#-------------------------------------------------------------------#

set.seed(1234)

n_soft_f <- length(zs_f)

ln_f <- runif(
  n_soft_f,
  min = fig_l - 0.5,
  max = fig_l + 0.5
)

u_f <- runif(
  n_soft_f,
  min = fig_q[1],
  max = fig_q[2]
)

a_f <- zs_f - ln_f * u_f
b_f <- zs_f + ln_f * (1 - u_f)

data_object_f <- bme_map(
  ch = ch_f,
  cs = cs_f,
  zh = zh_f,
  a = a_f,
  b = b_f
)

zk_rng_f <- base::range(zh_f, a_f, b_f)

#-------------------------------------------------------------------#
# 6. FIVE-FOLD CROSS-VALIDATION FOR CRBME
#-------------------------------------------------------------------#

set.seed(2026)

folds_f <- createFolds(
  zh_f,
  k = 5
)

cv_f <- vector(
  "list",
  nrow(param_grid)
)

for (g in seq_len(nrow(param_grid))) {

  fold_results <- vector(
    "list",
    length(folds_f)
  )

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
      observed = zh_f[test_idx],
      predicted = preds$mean
    )
  }

  combined_folds <- do.call(
    rbind,
    fold_results
  )

  RMSE <- sqrt(
    mean(
      (combined_folds$observed -
         combined_folds$predicted)^2
    )
  )

  cv_f[[g]] <- data.frame(
    lambda = param_grid$lambda[g],
    mu = param_grid$mu[g],
    RMSE = RMSE
  )
}

cv_f_df <- do.call(rbind, cv_f)

best_f <- cv_f_df[
  which.min(cv_f_df$RMSE),
  ,
  drop = FALSE
]

lambda_f <- best_f$lambda
mu_f <- best_f$mu

cat("\nSelected simulation setting:\n")
cat("l =", fig_l, "\n")
cat("q =", fig_q, "\n")
cat("Realization =", fig_j, "\n")
cat("Best lambda =", lambda_f, "\n")
cat("Best mu =", mu_f, "\n")
cat("CV RMSE =", best_f$RMSE, "\n")

#-------------------------------------------------------------------#
# 7. CREATE REGULAR PREDICTION GRID
#-------------------------------------------------------------------#

grid_xy <- expand.grid(
  x = seq(0, 100, length.out = n_grid),
  y = seq(0, 100, length.out = n_grid)
)

#-------------------------------------------------------------------#
# 8. INTERPOLATED REFERENCE FIELD
#-------------------------------------------------------------------#

# This is a kriged reference surface from all 625 simulated
# values, not the exact continuous simulated truth.

truth_sf <- st_as_sf(
  t_data,
  coords = c("x", "y")
)

grid_sf3 <- st_as_sf(
  grid_xy,
  coords = c("x", "y")
)

vg_truth <- fit.variogram(
  variogram(z ~ 1, data = truth_sf),
  model = vgm(c("Exp", "Sph"))
)

truth_grid <- krige(
  z ~ 1,
  locations = truth_sf,
  newdata = grid_sf3,
  model = vg_truth
)$var1.pred

#-------------------------------------------------------------------#
# 9. CBME AND CRBME PREDICTIONS
#-------------------------------------------------------------------#

message(
  "Predicting at ",
  nrow(grid_xy),
  " grid locations..."
)

cbme_output <- BMEmapping:::bme_estimate(
  grid_xy,
  data_object_f,
  model_f,
  nugget_f,
  sill_f,
  range_f,
  zk_range = zk_rng_f
)

crbme_output <- new_bme_estimate(
  grid_xy,
  data_object_f,
  zk_range = zk_rng_f,
  lambda = lambda_f,
  mu = mu_f
)

# Column 2 = mean prediction
pred_cbme <- cbme_output[, 2]
pred_crbme <- crbme_output[, 2]

#-------------------------------------------------------------------#
# 10. PREPARE MAP DATA
#-------------------------------------------------------------------#

map_df <- rbind(
  data.frame(
    grid_xy,
    value = truth_grid,
    panel = "Reference field"
  ),
  data.frame(
    grid_xy,
    value = pred_cbme,
    panel = "CBME"
  ),
  data.frame(
    grid_xy,
    value = pred_crbme,
    panel = "CRBME"
  )
)

map_df$panel <- factor(
  map_df$panel,
  levels = c(
    "Reference field",
    "CBME",
    "CRBME"
  )
)

lim <- base::range(
  map_df$value,
  na.rm = TRUE
)

#-------------------------------------------------------------------#
# 11. PROFESSIONAL ECOSTA PLOT
#-------------------------------------------------------------------#

theme_ecosta <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.6
    ),
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    axis.title = element_text(
      size = 13
    ),
    strip.background = element_rect(
      fill = "grey95",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 13
    ),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.spacing = grid::unit(0.5, "cm")
  )

p3 <- ggplot(
  map_df,
  aes(
    x = x,
    y = y,
    fill = value
  )
) +
  geom_raster(
    interpolate = TRUE
  ) +
  facet_wrap(
    ~panel,
    nrow = 1
  ) +
  scale_fill_viridis_c(
    option = "D",
    limits = lim,
    name = "z"
  ) +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100)
  ) +
  coord_fixed(
    xlim = c(0, 100),
    ylim = c(0, 100),
    expand = FALSE
  ) +
  labs(
    x = "x",
    y = "y"
  ) +
  theme_ecosta

print(p3)

#-------------------------------------------------------------------#
# 12. SAVE FIGURES
#-------------------------------------------------------------------#

dir.create(
  "figs",
  showWarnings = FALSE
)

ggsave(
  "figs/prediction_maps.pdf",
  plot = p3,
  width = 11,
  height = 4.2
)

ggsave(
  "figs/prediction_maps.png",
  plot = p3,
  width = 11,
  height = 4.2,
  dpi = 400,
  bg = "white"
)

# Save selected-case data to avoid recomputing the maps
saveRDS(
  list(
    map_df = map_df,
    lambda = lambda_f,
    mu = mu_f,
    fig_l = fig_l,
    fig_q = fig_q,
    fig_j = fig_j,
    cv_results = cv_f_df
  ),
  file = "figs/prediction_map_results.rds"
)

cat("\nPrediction maps saved successfully.\n")
