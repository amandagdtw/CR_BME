library(plotly)
library(ggplot2)
library(sf)
library(gstat)
library(tidyr)
library(viridis)
library(gridExtra)
library(kableExtra)
library(BMEmapping)
library(stars)
library(maps)
library(ggspatial)
library(devtools)

load_all()


###############################################################################
#                                  Data
###############################################################################

data_lower <- read.csv("final_table_lower.csv")

# data for upper bounds
data_upper <- read.csv("final_table_upper.csv")


# lower bounds data
lower_bounds <- data_lower %>%
  dplyr::filter(STATE == "KS") %>%
  dplyr::select(
    NAME,
    LATITUDE,
    LONGITUDE,
    ELEVATION,
    RT_II,
    ECO3
  )

# upper bounds data
upper_bounds <- data_upper %>%
  dplyr::filter(STATE == "KS") %>%
  dplyr::select(
    NAME,
    RT_II,
    ECO3
  )

# lower bounds and upper bounds data
data <- merge(
  x = lower_bounds,
  y = upper_bounds,
  by = "NAME"
) %>%
  transform(
    RT_II.x = ifelse(
      RT_II.x <= RT_II.y,
      RT_II.x,
      RT_II.y
    ),
    RT_II.y = ifelse(
      RT_II.x > RT_II.y,
      RT_II.x,
      RT_II.y
    )
  ) %>%
  mutate(
    range = RT_II.y - RT_II.x,
    center = 1 / 2 * (RT_II.y + RT_II.x),
    type = ifelse(range > 23, "SNWD", "WESD")
  ) %>%
  transform(
    RT_II.x = ifelse(range > 23, RT_II.x, center),
    RT_II.y = ifelse(range > 23, RT_II.y, center)
  ) %>%
  arrange(desc(type)) %>%
  dplyr::select(
    NAME,
    LATITUDE,
    LONGITUDE,
    ELEVATION,
    center,
    RT_II.x,
    RT_II.y,
    ECO3.x,
    type
  )

names(data) <- c(
  "station",
  "latitude",
  "longitude",
  "elevation",
  "center",
  "lower",
  "upper",
  "eco3",
  "type"
)

dim(data)
table(data$type)

ggplot(
  data,
  aes(x = elevation, y = (center))
) +
  geom_point(
    aes(color = type, shape = type),
    size = 1.5,
    alpha = 0.8
  )


###############################################################################
#                              Data Analysis
###############################################################################


#-----------------------------------------------------------------------------
#                              kansas map
#-----------------------------------------------------------------------------

# Get kansas map polygon and convert to sf format
kansas <- map_data("state", region = "kansas")

ggplot() +
  geom_polygon(
    data = kansas,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "black"
  ) +
  geom_point(
    data = data,
    aes(
      x = longitude,
      y = latitude,
      color = type,
      shape = type
    ),
    size = 3,
    show.legend = TRUE
  ) +
  coord_quickmap() +
  scale_x_continuous(
    name = "longitude",
    labels = function(x)
      paste0(
        abs(x),
        "° ",
        ifelse(x < 0, "W", "E")
      )
  ) +
  scale_y_continuous(
    name = "latitude",
    labels = function(y)
      paste0(
        abs(y),
        "° ",
        ifelse(y < 0, "S", "N")
      )
  ) +
  annotation_north_arrow(
    location = "tr", # bottom-left
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.2, "cm"),
    width = unit(1.2, "cm")
  ) +
  labs(
    color = "measurement",
    shape = "measurement"
  ) +
  theme(
    text = element_text(size = 18),

    # Axis labels
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),

    # Legend
    legend.title = element_text(size = 19),
    legend.text = element_text(size = 18),

    panel.background = element_rect(
      fill = "white",
      color = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#-----------------------------------------------------------------------------
#                         Summary statistics table
#-----------------------------------------------------------------------------

# Hard and soft subsets
hard_data <- subset(data, type == "WESD")
soft_data <- subset(data, type == "SNWD")

# Function to compute summary statistics
summ_stats <- function(x) {
  c(
    n = length(x),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    iqr = IQR(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

# Compute statistics
center_h <- summ_stats(hard_data$center)
center_u <- summ_stats(soft_data$upper)
center_l <- summ_stats(soft_data$lower)
elev <- summ_stats(data$elevation)

# Create table
summary_table <- data.frame(
  Variable = c(
    "center_h",
    "center_su",
    "center_sl",
    "elevation"
  ),

  n = c(
    center_h["n"],
    center_u["n"],
    center_l["n"],
    elev["n"]
  ),

  `Mean ± SD` = c(
    sprintf(
      "%.1f ± %.1f",
      center_h["mean"],
      center_h["sd"]
    ),
    sprintf(
      "%.1f ± %.1f",
      center_u["mean"],
      center_u["sd"]
    ),
    sprintf(
      "%.1f ± %.1f",
      center_l["mean"],
      center_l["sd"]
    ),
    sprintf(
      "%.1f ± %.1f",
      elev["mean"],
      elev["sd"]
    )
  ),

  `Median (IQR)` = c(
    sprintf(
      "%.1f (%.1f)",
      center_h["median"],
      center_h["iqr"]
    ),
    sprintf(
      "%.1f (%.1f)",
      center_u["median"],
      center_u["iqr"]
    ),
    sprintf(
      "%.1f (%.1f)",
      center_l["median"],
      center_l["iqr"]
    ),
    sprintf(
      "%.1f (%.1f)",
      elev["median"],
      elev["iqr"]
    )
  ),

  `Min-Max` = c(
    sprintf(
      "%.1f - %.1f",
      center_h["min"],
      center_h["max"]
    ),
    sprintf(
      "%.1f - %.1f",
      center_u["min"],
      center_u["max"]
    ),
    sprintf(
      "%.1f - %.1f",
      center_l["min"],
      center_l["max"]
    ),
    sprintf(
      "%.1f - %.1f",
      elev["min"],
      elev["max"]
    )
  ),

  check.names = FALSE
)

summary_table



#-----------------------------------------------------------------------------
#                            Variogram fitting
#-----------------------------------------------------------------------------

SDD <- sf::st_as_sf(
  data,
  coords = c("longitude", "latitude")
)

vgm_emp <- variogram(
  center ~ 1,
  SDD
)

cutoff <- 1.5 *
  vgm_emp[which.max(vgm_emp$gamma), "dist"]

vgm_emp <- variogram(
  center ~ 1,
  cutoff = cutoff,
  SDD
)

vg_model <- fit.variogram(
  vgm_emp,
  model = vgm(c("Exp", "Sph"))
)

plot(vgm_emp, vg_model)

vg_model


# ============================================================
# Data partitioning
# ============================================================
SD <- data

nh <- as.vector(table(SD$type)[2])

# hard data indices
idx_h <- 1:nh
idx_s <- (nh + 1):(dim(SD)[1] - 6)

# hard data locations
ch <- SD[idx_h, c("longitude", "latitude")]

# hard data values
zh <- SD[idx_h, "center"]

# soft data locations
cs <- SD[idx_s, c("longitude", "latitude")]

# lower and upper bounds of soft data (intervals)
a <- SD[idx_s, "lower"]
b <- SD[idx_s, "upper"]


# ============================================================
# Create BME mapping object
# ============================================================

data_object <- bme_map(ch, cs, zh, a, b)


# ============================================================
# Spatial plot of data (Figure 1)
# ============================================================

plot(data_object)


# ============================================================
# Estimation locations
# ============================================================

idx_k <- (dim(SD)[1] - 5):dim(SD)[1]

xk <- SD[idx_k, c("longitude", "latitude")]

xk


# ============================================================
# Posterior density estimation and plots selected locations
# using default zk_range
# ============================================================

p_1 <- prob_zk(
  xk[1, ],
  data_object,
  model,
  nugget,
  sill,
  range
)

p_2 <- prob_zk(
  xk[2, ],
  data_object,
  model,
  nugget,
  sill,
  range
)

p_3 <- prob_zk(
  xk[3, ],
  data_object,
  model,
  nugget,
  sill,
  range
)

p_4 <- prob_zk(
  xk[4, ],
  data_object,
  model,
  nugget,
  sill,
  range
)

p_5 <- prob_zk(
  xk[5, ],
  data_object,
  model,
  nugget,
  sill,
  range
)

p_6 <- prob_zk(
  xk[6, ],
  data_object,
  model,
  nugget,
  sill,
  range
)

p_df <- cbind.data.frame(
  p_1,
  p_2[, 2],
  p_3[, 2],
  p_4[, 2],
  p_5[, 2],
  p_6[, 2]
)

names(p_df) <- c(
  "zk_i",
  "p1",
  "p2",
  "p3",
  "p4",
  "p5",
  "p6"
)


# ============================================================
# Function to generate ggplot for a given column name
# ============================================================

# Function to generate ggplot for a given column name
plot_prob_curve <- function(df, pi_col) {
  ggplot(df, aes(x = zk_i, y = .data[[pi_col]])) +
    geom_line(color = "darkblue", linewidth = 0.7) +
    labs(x = "Spatial Field Value (z)", y = "Posterior Density f(z)") +
    theme_minimal() +
    theme(
      panel.background = ggplot2::element_rect(fill = "white", color = "black")
    )
}


# ============================================================
# Generate individual plots
# ============================================================

p1 <- plot_prob_curve(p_df, "p1")
p2 <- plot_prob_curve(p_df, "p2")
p3 <- plot_prob_curve(p_df, "p3")
p4 <- plot_prob_curve(p_df, "p4")
p5 <- plot_prob_curve(p_df, "p5")
p6 <- plot_prob_curve(p_df, "p6")


# ============================================================
# Arrange in a 2x2 grid
# ============================================================

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  ncol = 2
)


# ============================================================
# CBME vs CRBME Comparison
# ============================================================

# CBME

cbme_mode <- bme_cv(
  data_object,
  model,
  nugget,
  sill,
  range,
  type = "mode",
  zk_range = c(10, 40)
)

cbme_mean <- bme_cv(
  data_object,
  model,
  nugget,
  sill,
  range,
  type = "mean",
  zk_range = c(10, 40)
)


# ============================================================
# 5-FOLD CROSS-VALIDATION VIA PARAMETER GRID SEARCH
# ============================================================

best_lambda <- 2
best_mu <- 0.5


q_1 <- new_prob_zk(
  xk[1, ],
  data_object,
  lambda = best_lambda,
  mu = best_mu
)

q_2 <- new_prob_zk(
  xk[2, ],
  data_object,
  lambda = best_lambda,
  mu = best_mu
)

q_3 <- new_prob_zk(
  xk[3, ],
  data_object,
  lambda = best_lambda,
  mu = best_mu
)

q_4 <- new_prob_zk(
  xk[4, ],
  data_object,
  lambda = best_lambda,
  mu = best_mu
)

q_5 <- new_prob_zk(
  xk[5, ],
  data_object,
  lambda = best_lambda,
  mu = best_mu
)

q_6 <- new_prob_zk(
  xk[6, ],
  data_object,
  lambda = best_lambda,
  mu = best_mu
)

q_df <- cbind.data.frame(
  q_1,
  q_2[, 2],
  q_3[, 2],
  q_4[, 2],
  q_5[, 2],
  q_6[, 2]
)

names(q_df) <- c(
  "zk_i",
  "p1",
  "p2",
  "p3",
  "p4",
  "p5",
  "p6"
)


# ============================================================
# Generate individual plots
# ============================================================

q1 <- plot_prob_curve(q_df, "p1")
q2 <- plot_prob_curve(q_df, "p2")
q3 <- plot_prob_curve(q_df, "p3")
q4 <- plot_prob_curve(q_df, "p4")
q5 <- plot_prob_curve(q_df, "p5")
q6 <- plot_prob_curve(q_df, "p6")


# ============================================================
# Arrange in a 2x2 grid
# ============================================================

grid.arrange(
  q1, q2, q3,
  q4, q5, q6,
  ncol = 2
)


# ============================================================
# CRBME
# ============================================================

best_lambda <- 1
best_mu <- 0.5

crbme_mode <- new_bme_cv(
  data_object,
  lambda = best_lambda,
  mu = best_mu,
  type = "mode",
  zk_range = c(0, 50)
)

crbme_mean <- new_bme_cv(
  data_object,
  lambda = best_lambda,
  mu = best_mu,
  type = "mean",
  zk_range = c(0, 50)
)


# ============================================================
# Summary
# ============================================================

CBME_mode <- summary(cbme_mode)
CBME_mean <- summary(cbme_mean)

CRBME_mode <- summary(crbme_mode)
CRBME_mean <- summary(crbme_mean)


# ============================================================
# Extract metrics
# ============================================================

get_metrics <- function(df) {
  setNames(df$Value, df$Metric)
}


cbme_mode_get <- get_metrics(CBME_mode)
crbme_mode_get <- get_metrics(CRBME_mode)

cbme_mean_get <- get_metrics(CBME_mean)
crbme_mean_get <- get_metrics(CRBME_mean)


# ============================================================
# Create summary table
# ============================================================

results_table <- data.frame(

  Prediction = c(
    "BME Mode",
    "",
    "",
    "BME Mean",
    ""
  ),

  Approach = c(
    "CBME",
    "CRBME",
    "",
    "CBME",
    "CRBME"
  ),

  ME = c(
    cbme_mode_get["ME"],
    crbme_mode_get["ME"],
    NA,
    cbme_mean_get["ME"],
    crbme_mean_get["ME"]
  ),

  MAE = c(
    cbme_mode_get["MAE"],
    crbme_mode_get["MAE"],
    NA,
    cbme_mean_get["MAE"],
    crbme_mean_get["MAE"]
  ),

  RMSE = c(
    cbme_mode_get["RMSE"],
    crbme_mode_get["RMSE"],
    NA,
    cbme_mean_get["RMSE"],
    crbme_mean_get["RMSE"]
  ),

  Rsquared = c(
    cbme_mode_get["Rsquared"],
    crbme_mode_get["Rsquared"],
    NA,
    cbme_mean_get["Rsquared"],
    crbme_mean_get["Rsquared"]
  )
)


# ============================================================
# Round metrics
# ============================================================

results_table <- results_table %>%
  mutate(
    ME = round(ME, 4),
    MAE = round(MAE, 3),
    RMSE = round(RMSE, 3),
    Rsquared = round(Rsquared, 3)
  )

# ============================================================
# Display results
# ============================================================

results_table


