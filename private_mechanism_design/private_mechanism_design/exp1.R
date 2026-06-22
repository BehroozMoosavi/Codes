# ================================================================
# EXPERIMENT 1
# Posterior-score panel across four privacy channels
#
# Optimal Binary Mechanism under Locally Private Signals
#
# The script plots the deterministic single-signal posterior score
#
#     S_i(y, lambda)
#       = E[X_i | Y_i = y]
#         + lambda E[J(X_i) | Y_i = y]
#
# under
#
#     X_i ~ Uniform[-1,1].
#
# For a Uniform[a,b] prior,
#
#     J(x) = x - (1 - F(x))/f(x) = 2x - b.
#
# With a = -1 and b = 1,
#
#     J(x) = 2x - 1.
#
# CHANNELS
#
#   1. Gaussian:
#          Y = X + N(0, 0.4^2).
#
#   2. Laplace:
#          Y = X + Laplace(0, 0.5).
#
#   3. Logistic:
#          Y = X + Logistic(0, 0.5).
#
# INTERPRETATION
#
#   * The Gaussian and logistic channels have strictly
#     increasing likelihood ratios and therefore strictly increasing
#     posterior expectations of increasing functions.
#
#   * The Laplace channel has weak MLRP. Once y >= 1 or y <= -1,
#     the posterior is exactly independent of y, so the score has
#     exactly flat global tails.
#
# OUTPUT
#
#   figures/fig_exp1_posterior_score_panel.pdf
#   tables/exp1_posterior_score_curves.csv
#   tables/exp1_zero_crossings.csv
#   tables/exp1_sanity_report.csv
#
# Run from top to bottom.
# ================================================================


rm(list = ls())


# ================================================================
# 0. USER OPTIONS
# ================================================================

SHOW_TITLE <- FALSE

out_dir <- paste0(
  "/Users/behroozmoosavi/Desktop/privacy/codes/",
  "R_codes/experiment/exp1"
)

fig_dir  <- file.path(out_dir, "figures")
data_dir <- file.path(out_dir, "tables")

for (d in c(out_dir, fig_dir, data_dir)) {
  if (!dir.exists(d)) {
    dir.create(
      d,
      recursive = TRUE
    )
  }
}


# ================================================================
# 1. PACKAGES
# ================================================================

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


# ================================================================
# 2. PLOT THEME AND SAVE FUNCTION
# ================================================================

theme_paper <- function(base_size = 13) {
  theme_bw(
    base_size = base_size,
    base_family = "serif"
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = base_size,
        hjust = 0.5
      ),
      axis.title = element_text(
        size = base_size + 1
      ),
      axis.text = element_text(
        size = base_size - 1
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(
        size = base_size
      ),
      panel.grid.major = element_line(
        color = "gray88",
        linewidth = 0.30
      ),
      panel.grid.minor = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = base_size
      ),
      strip.background = element_rect(
        fill = "gray93",
        color = "black",
        linewidth = 0.35
      ),
      panel.border = element_rect(
        color = "black",
        linewidth = 0.45
      ),
      plot.margin = margin(
        6,
        9,
        6,
        6
      )
    )
}


save_figure <- function(
    plot_object,
    filename_base,
    width = 7.6,
    height = 6.4
) {
  output_file <- file.path(
    fig_dir,
    paste0(
      filename_base,
      ".pdf"
    )
  )
  
  ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  
  invisible(
    output_file
  )
}


# ================================================================
# 3. PRIOR AND VIRTUAL VALUE
# ================================================================

x_lower <- -1
x_upper <-  1

prior_density <- function(x) {
  ifelse(
    x >= x_lower &
      x <= x_upper,
    1 / (x_upper - x_lower),
    0
  )
}


virtual_value <- function(x) {
  2 * x -
    x_upper
}


# ================================================================
# 4. CHANNEL PARAMETERS
# ================================================================

gaussian_sd <- 0.40

laplace_scale <- 0.50
logistic_scale <- 0.50


channels <- c(
  "Gaussian",
  "Laplace",
  "Logistic"
)


channel_display_names <- c(
  "Gaussian",
  "Laplace",
  "Logistic"
)


# ================================================================
# 5. CONDITIONAL LOG-DENSITIES
#
# Log-densities are used rather than densities to improve numerical
# stability in the signal tails.
# ================================================================

log_kernel_gaussian <- function(
    y,
    x
) {
  dnorm(
    y,
    mean = x,
    sd = gaussian_sd,
    log = TRUE
  )
}


log_kernel_laplace <- function(
    y,
    x
) {
  -log(
    2 *
      laplace_scale
  ) -
    abs(
      y -
        x
    ) /
    laplace_scale
}


log_kernel_logistic <- function(
    y,
    x
) {
  dlogis(
    y -
      x,
    location = 0,
    scale = logistic_scale,
    log = TRUE
  )
}


get_log_kernel <- function(channel) {
  switch(
    channel,
    "Gaussian" =
      log_kernel_gaussian,
    "Laplace" =
      log_kernel_laplace,
    "Logistic" =
      log_kernel_logistic,
    stop(
      paste(
        "Unknown channel:",
        channel
      )
    )
  )
}


# ================================================================
# 6. STABLE POSTERIOR MOMENTS
# ================================================================

posterior_moments <- function(
    signal_grid,
    channel,
    type_grid
) {
  dx <- type_grid[2] -
    type_grid[1]
  
  prior_values <- prior_density(
    type_grid
  )
  
  virtual_values <- virtual_value(
    type_grid
  )
  
  log_kernel <- get_log_kernel(
    channel
  )
  
  rows <- lapply(
    signal_grid,
    function(y) {
      log_weights <- log_kernel(
        y,
        type_grid
      )
      
      positive_prior <- prior_values > 0
      
      log_weights[positive_prior] <-
        log_weights[positive_prior] +
        log(
          prior_values[positive_prior]
        )
      
      log_weights[!positive_prior] <-
        -Inf
      
      maximum_log_weight <- max(
        log_weights
      )
      
      if (!is.finite(maximum_log_weight)) {
        return(
          data.frame(
            y = y,
            xhat = NA_real_,
            Jhat = NA_real_
          )
        )
      }
      
      stabilized_weights <- exp(
        log_weights -
          maximum_log_weight
      )
      
      denominator <- sum(
        stabilized_weights
      ) *
        dx
      
      if (
        !is.finite(denominator) ||
        denominator <= .Machine$double.eps
      ) {
        return(
          data.frame(
            y = y,
            xhat = NA_real_,
            Jhat = NA_real_
          )
        )
      }
      
      posterior_x <- sum(
        type_grid *
          stabilized_weights
      ) *
        dx /
        denominator
      
      posterior_J <- sum(
        virtual_values *
          stabilized_weights
      ) *
        dx /
        denominator
      
      data.frame(
        y = y,
        xhat = posterior_x,
        Jhat = posterior_J
      )
    }
  )
  
  bind_rows(
    rows
  )
}


# ================================================================
# 7. GRIDS AND MULTIPLIERS
# ================================================================

type_grid <- seq(
  x_lower,
  x_upper,
  length.out = 5001
)

# The range [-4,4] displays the exact Laplace plateaus and gives
# a clear view of the logistic limiting behavior.
signal_grid <- seq(
  -4,
  4,
  length.out = 2001
)

lambda_grid <- c(
  0,
  0.3,
  0.7,
  1.5
)


# ================================================================
# 8. BUILD POSTERIOR-SCORE CURVES
# ================================================================

build_channel_scores <- function(channel) {
  posterior_table <- posterior_moments(
    signal_grid = signal_grid,
    channel = channel,
    type_grid = type_grid
  )
  
  expand_grid(
    y = posterior_table$y,
    lambda = lambda_grid
  ) %>%
    left_join(
      posterior_table,
      by = "y"
    ) %>%
    mutate(
      score =
        xhat +
        lambda *
        Jhat,
      lambda_factor = factor(
        lambda,
        levels = lambda_grid
      ),
      channel = channel
    )
}


score_data <- bind_rows(
  lapply(
    channels,
    build_channel_scores
  )
) %>%
  mutate(
    channel = factor(
      channel,
      levels = channels,
      labels = channel_display_names
    )
  )


# ================================================================
# 9. ZERO-CROSSING FUNCTION
# ================================================================

find_first_zero_crossing <- function(
    y,
    score,
    tolerance = 1e-12
) {
  valid <- is.finite(y) &
    is.finite(score)
  
  y <- y[valid]
  score <- score[valid]
  
  ordering <- order(
    y
  )
  
  y <- y[ordering]
  score <- score[ordering]
  
  exact_zero <- which(
    abs(score) <= tolerance
  )
  
  if (length(exact_zero) > 0) {
    return(
      y[exact_zero[1]]
    )
  }
  
  sign_change <- which(
    score[-length(score)] *
      score[-1] <
      0
  )
  
  if (length(sign_change) == 0) {
    return(
      NA_real_
    )
  }
  
  index <- sign_change[1]
  
  y_left <- y[index]
  y_right <- y[index + 1]
  
  score_left <- score[index]
  score_right <- score[index + 1]
  
  y_left -
    score_left *
    (
      y_right -
        y_left
    ) /
    (
      score_right -
        score_left
    )
}


zero_crossings <- score_data %>%
  group_by(
    channel,
    lambda_factor,
    lambda
  ) %>%
  summarize(
    zero_y = find_first_zero_crossing(
      y,
      score
    ),
    zero_score = 0,
    .groups = "drop"
  )


# ================================================================
# 10. ANALYTIC TAIL REFERENCES
# ================================================================

# ------------------------------------------------
# Laplace and logistic right-tail limit
#
# For either additive family with scale s,
#
#     p(x | y -> +infinity)
#       proportional to exp(x/s) f(x).
#
# Under Uniform[-1,1],
#
#     E[X | y -> +infinity]
#       = coth(1/s) - s.
#
# The left-tail limit is its negative.
# ------------------------------------------------

exponential_tilt_mean <- function(scale) {
  theta <- 1 /
    scale
  
  1 /
    tanh(theta) -
    1 /
    theta
}


laplace_right_x_limit <- exponential_tilt_mean(
  laplace_scale
)

laplace_left_x_limit <- -laplace_right_x_limit

logistic_right_x_limit <- exponential_tilt_mean(
  logistic_scale
)

logistic_left_x_limit <- -logistic_right_x_limit


# ------------------------------------------------
# Gaussian tail limit
#
# For the additive Gaussian channel Y = X + N(0, sigma^2) the
# conditional mean E[Y | X = x] = x is strictly increasing in x,
# so as y -> +/- infinity the posterior concentrates at the
# corresponding support endpoint.
# ------------------------------------------------

gaussian_right_x_limit <- x_upper
gaussian_left_x_limit  <- x_lower


tail_reference <- function(
    channel,
    lambda
) {
  x_right <- switch(
    channel,
    "Gaussian" =
      gaussian_right_x_limit,
    "Laplace" =
      laplace_right_x_limit,
    "Logistic" =
      logistic_right_x_limit
  )
  
  x_left <- switch(
    channel,
    "Gaussian" =
      gaussian_left_x_limit,
    "Laplace" =
      laplace_left_x_limit,
    "Logistic" =
      logistic_left_x_limit
  )
  
  J_right <- virtual_value(
    x_right
  )
  
  J_left <- virtual_value(
    x_left
  )
  
  data.frame(
    theoretical_xhat_left = x_left,
    theoretical_xhat_right = x_right,
    theoretical_score_left =
      x_left +
      lambda *
      J_left,
    theoretical_score_right =
      x_right +
      lambda *
      J_right
  )
}


# ================================================================
# 11. SANITY REPORT
# ================================================================

sanity_report <- score_data %>%
  group_by(
    channel,
    lambda
  ) %>%
  summarize(
    numerical_y_left = min(
      y
    ),
    numerical_y_right = max(
      y
    ),
    numerical_xhat_left =
      xhat[which.min(y)],
    numerical_xhat_right =
      xhat[which.max(y)],
    numerical_score_left =
      score[which.min(y)],
    numerical_score_right =
      score[which.max(y)],
    minimum_score_increment =
      min(
        diff(
          score[order(y)]
        ),
        na.rm = TRUE
      ),
    monotone_increasing =
      minimum_score_increment >=
      -1e-8,
    .groups = "drop"
  )


reference_rows <- bind_rows(
  lapply(
    seq_len(
      nrow(sanity_report)
    ),
    function(index) {
      tail_reference(
        channel = as.character(
          sanity_report$channel[index]
        ),
        lambda = sanity_report$lambda[index]
      )
    }
  )
)


sanity_report <- bind_cols(
  sanity_report,
  reference_rows
) %>%
  mutate(
    left_score_error =
      numerical_score_left -
      theoretical_score_left,
    right_score_error =
      numerical_score_right -
      theoretical_score_right
  )


cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT 1 SANITY REPORT\n"
)

cat(
  "============================================================\n\n"
)

print(
  as.data.frame(
    sanity_report
  ),
  digits = 5,
  row.names = FALSE
)


cat(
  "\nAnalytic posterior-mean tail limits:\n"
)

cat(
  sprintf(
    "  Gaussian:  left = %.6f, right = %.6f\n",
    gaussian_left_x_limit,
    gaussian_right_x_limit
  )
)

cat(
  sprintf(
    "  Laplace:          left = %.6f, right = %.6f\n",
    laplace_left_x_limit,
    laplace_right_x_limit
  )
)

cat(
  sprintf(
    "  Logistic:         left = %.6f, right = %.6f\n",
    logistic_left_x_limit,
    logistic_right_x_limit
  )
)


cat(
  paste0(
    "\nInterpretation:\n",
    "  * Laplace is exactly flat for y <= -1 and y >= 1.\n",
    "  * Logistic remains strictly increasing at every finite y,\n",
    "    but converges to the same exponentially tilted limits as\n",
    "    Laplace when the two scales are equal.\n",
    "  * Gaussian approaches the endpoint types as |y| grows.\n\n"
  )
)


# ================================================================
# 12. FIGURE
# ================================================================

lambda_colors <- c(
  "0"   = "#0072B2",
  "0.3" = "#1B9E77",
  "0.7" = "#E6AB02",
  "1.5" = "#D55E00"
)


lambda_labels <- c(
  expression(lambda == 0),
  expression(lambda == 0.3),
  expression(lambda == 0.7),
  expression(lambda == 1.5)
)


score_plot <- ggplot(
  score_data,
  aes(
    x = y,
    y = score,
    color = lambda_factor
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray25",
    linewidth = 0.45
  ) +
  geom_line(
    linewidth = 0.95
  ) +
  geom_point(
    data = zero_crossings,
    aes(
      x = zero_y,
      y = zero_score,
      color = lambda_factor
    ),
    size = 2,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ channel,
    ncol = 2
  ) +
  scale_color_manual(
    values = lambda_colors,
    labels = lambda_labels
  ) +
  coord_cartesian(
    xlim = c(
      -4,
      4
    ),
    ylim = c(
      -6,
      2
    )
  ) +
  labs(
    x = expression(
      "signal " *
        y[i]
    ),
    y = expression(
      "posterior score  " *
        S[i](
          y[i],
          lambda
        )
    )
  ) +
  theme_paper()


if (SHOW_TITLE) {
  score_plot <- score_plot +
    labs(
      title = "Posterior scores across privacy channels"
    )
}


save_figure(
  plot_object = score_plot,
  filename_base = "fig_exp1_posterior_score_panel",
  width = 7.6,
  height = 6.4
)


# ================================================================
# 13. SAVE TABLES
# ================================================================

write.csv(
  score_data %>%
    select(
      channel,
      lambda,
      y,
      xhat,
      Jhat,
      score
    ) %>%
    arrange(
      channel,
      lambda,
      y
    ),
  file.path(
    data_dir,
    "exp1_posterior_score_curves.csv"
  ),
  row.names = FALSE
)


write.csv(
  zero_crossings %>%
    select(
      channel,
      lambda,
      zero_y
    ) %>%
    mutate(
      zero_y = round(
        zero_y,
        6
      )
    ) %>%
    arrange(
      channel,
      lambda
    ),
  file.path(
    data_dir,
    "exp1_zero_crossings.csv"
  ),
  row.names = FALSE
)


write.csv(
  sanity_report,
  file.path(
    data_dir,
    "exp1_sanity_report.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 14. COMPLETION REPORT
# ================================================================

cat(
  "\nZero crossings (pivotal signals):\n"
)

print(
  zero_crossings %>%
    select(
      channel,
      lambda,
      zero_y
    ) %>%
    arrange(
      channel,
      lambda
    ) %>%
    as.data.frame(),
  digits = 5,
  row.names = FALSE
)


cat(
  "\nFigure written to:\n  ",
  file.path(
    fig_dir,
    "fig_exp1_posterior_score_panel.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "Tables written to:\n  ",
  data_dir,
  "\n",
  sep = ""
)