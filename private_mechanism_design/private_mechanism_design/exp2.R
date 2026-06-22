# ================================================================
# EXPERIMENT 2
# Knife-edge transition in the implementation decision
#
# PRIOR
#
#     X ~ Uniform[-0.5, 1.5].
#
# For a Uniform[a,b] prior,
#
#     J(x) = 2x - b,
#     E[X] = (a+b)/2,
#     E[J(X)] = a.
#
# Therefore,
#
#     E[X] = 0.5,
#     E[J(X)] = -0.5,
#
# and
#
#     E[S(Y,lambda)]
#       = E[X] + lambda E[J(X)]
#       = 0.5 - 0.5 lambda.
#
# The critical multiplier is
#
#     lambda* = 1.
#
# The implementation rule is
#
#     q(Y_1,...,Y_n)
#       = 1{
#           sum_i [
#             E[X_i | Y_i]
#             + lambda E[J(X_i) | Y_i]
#           ] >= 0
#         }.
#
# The CLT approximation is
#
#     P(build)
#       approximately
#       Phi(
#         sqrt(n) mu_S(lambda) /
#         Sigma_Z(lambda)
#       ).
#
# OUTPUT
#
#   figures/fig_exp2_knife_edge_transition.pdf
#   figures/fig_exp2_knife_edge_collapse.pdf
#
#   tables/exp2_knife_edge.csv
#   tables/exp2_dispersion.csv
#   tables/exp2_transition_sanity.csv
#
# ================================================================


rm(list = ls())


# ================================================================
# 0. USER OPTIONS
# ================================================================

FAST        <- FALSE
SHOW_TITLES <- FALSE

set.seed(20260609)


out_dir <- paste0(
  "/Users/behroozmoosavi/Desktop/privacy/codes/",
  "R_codes/experiment/exp2"
)

fig_dir  <- file.path(out_dir, "figures")
data_dir <- file.path(out_dir, "tables")


for (directory in c(out_dir, fig_dir, data_dir)) {
  if (!dir.exists(directory)) {
    dir.create(
      directory,
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


for (package_name in packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(
      package_name,
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
      legend.title = element_text(
        size = base_size
      ),
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
    height = 6.0
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
  
  invisible(output_file)
}


# ================================================================
# 3. PRIOR AND CRITICAL MULTIPLIER
# ================================================================

x_lower <- -0.5
x_upper <-  1.5


prior_mean <- (
  x_lower +
    x_upper
) / 2


mean_virtual_value <- x_lower


mean_score <- function(lambda_value) {
  prior_mean +
    lambda_value *
    mean_virtual_value
}


lambda_star <- -prior_mean /
  mean_virtual_value


stopifnot(
  abs(lambda_star - 1) <
    1e-12
)


uniform_density <- function(
    x,
    lower,
    upper
) {
  ifelse(
    x >= lower &
      x <= upper,
    1 / (upper - lower),
    0
  )
}


uniform_virtual_value <- function(
    x,
    upper
) {
  2 * x -
    upper
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


channel_display_map <- setNames(
  channel_display_names,
  channels
)


# ================================================================
# 5. CONDITIONAL LOG-DENSITIES
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


get_log_kernel <- function(channel_name) {
  switch(
    channel_name,
    
    "Gaussian" =
      log_kernel_gaussian,
    
    "Laplace" =
      log_kernel_laplace,
    
    "Logistic" =
      log_kernel_logistic,
    
    stop(
      paste(
        "Unknown channel:",
        channel_name
      )
    )
  )
}


# ================================================================
# 6. CHANNEL SIMULATION USING COMMON UNIFORMS
# ================================================================

safe_probability <- function(
    probability,
    tolerance = 1e-14
) {
  pmin(
    pmax(
      probability,
      tolerance
    ),
    1 - tolerance
  )
}


laplace_quantile <- function(
    uniform_draw,
    scale
) {
  uniform_draw <- safe_probability(
    uniform_draw
  )
  
  ifelse(
    uniform_draw < 0.5,
    
    scale *
      log(
        2 *
          uniform_draw
      ),
    
    -scale *
      log(
        2 *
          (
            1 -
              uniform_draw
          )
      )
  )
}


simulate_channel_from_uniforms <- function(
    x,
    uniform_noise,
    channel_name
) {
  uniform_noise <- safe_probability(
    uniform_noise
  )
  
  switch(
    channel_name,
    
    "Gaussian" =
      x +
      qnorm(
        uniform_noise,
        mean = 0,
        sd = gaussian_sd
      ),
    
    "Laplace" =
      x +
      laplace_quantile(
        uniform_noise,
        laplace_scale
      ),
    
    "Logistic" =
      x +
      qlogis(
        uniform_noise,
        location = 0,
        scale = logistic_scale
      ),
    
    stop(
      paste(
        "Unknown channel:",
        channel_name
      )
    )
  )
}


# ================================================================
# 7. POSTERIOR MOMENTS BY STABLE QUADRATURE
# ================================================================

posterior_moment_table <- function(
    signal_grid,
    channel_name,
    type_grid
) {
  dx <- type_grid[2] -
    type_grid[1]
  
  prior_values <- uniform_density(
    type_grid,
    x_lower,
    x_upper
  )
  
  virtual_values <- uniform_virtual_value(
    type_grid,
    x_upper
  )
  
  log_kernel <- get_log_kernel(
    channel_name
  )
  
  output_rows <- lapply(
    signal_grid,
    function(signal_value) {
      log_weights <- log_kernel(
        signal_value,
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
            y = signal_value,
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
            y = signal_value,
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
        y = signal_value,
        xhat = posterior_x,
        Jhat = posterior_J
      )
    }
  )
  
  bind_rows(
    output_rows
  )
}


make_interpolator <- function(
    posterior_table,
    column_name
) {
  force(
    posterior_table
  )
  
  force(
    column_name
  )
  
  function(signal_values) {
    approx(
      x = posterior_table$y,
      y = posterior_table[[column_name]],
      xout = signal_values,
      rule = 2,
      ties = "ordered"
    )$y
  }
}


# ================================================================
# 8. GRIDS AND POSTERIOR LOOKUP TABLES
# ================================================================

type_grid <- seq(
  x_lower,
  x_upper,
  length.out = if (FAST) {
    2001
  } else {
    4001
  }
)


signal_grid <- seq(
  -10,
  10,
  length.out = if (FAST) {
    3001
  } else {
    6001
  }
)


posterior_tables <- setNames(
  lapply(
    channels,
    function(channel_name) {
      posterior_moment_table(
        signal_grid = signal_grid,
        channel_name = channel_name,
        type_grid = type_grid
      )
    }
  ),
  channels
)


posterior_x_functions <- setNames(
  lapply(
    channels,
    function(channel_name) {
      make_interpolator(
        posterior_tables[[channel_name]],
        "xhat"
      )
    }
  ),
  channels
)


posterior_J_functions <- setNames(
  lapply(
    channels,
    function(channel_name) {
      make_interpolator(
        posterior_tables[[channel_name]],
        "Jhat"
      )
    }
  ),
  channels
)


# ================================================================
# 9. ESTIMATE SINGLE-AGENT SCORE DISPERSION
# ================================================================

dispersion_draws <- if (FAST) {
  50000
} else {
  200000
}


dispersion_type_draws <- runif(
  dispersion_draws,
  min = x_lower,
  max = x_upper
)


dispersion_uniform_noise <- runif(
  dispersion_draws
)


dispersion_rows <- lapply(
  channels,
  function(channel_name) {
    signals <- simulate_channel_from_uniforms(
      x = dispersion_type_draws,
      uniform_noise = dispersion_uniform_noise,
      channel_name = channel_name
    )
    
    posterior_x <- posterior_x_functions[[channel_name]](
      signals
    )
    
    posterior_J <- posterior_J_functions[[channel_name]](
      signals
    )
    
    data.frame(
      channel = channel_name,
      
      var_xhat = var(
        posterior_x,
        na.rm = TRUE
      ),
      
      cov_xhat_Jhat = cov(
        posterior_x,
        posterior_J,
        use = "complete.obs"
      ),
      
      var_Jhat = var(
        posterior_J,
        na.rm = TRUE
      ),
      
      mean_xhat = mean(
        posterior_x,
        na.rm = TRUE
      ),
      
      mean_Jhat = mean(
        posterior_J,
        na.rm = TRUE
      ),
      
      stringsAsFactors = FALSE
    )
  }
)


dispersion_table <- bind_rows(
  dispersion_rows
)


score_sd <- function(
    channel_name,
    lambda_value
) {
  selected_row <- dispersion_table[
    dispersion_table$channel ==
      channel_name,
    ,
    drop = FALSE
  ]
  
  if (nrow(selected_row) != 1L) {
    stop(
      sprintf(
        paste0(
          "score_sd() expected exactly one row for channel '%s', ",
          "but found %d rows."
        ),
        channel_name,
        nrow(selected_row)
      )
    )
  }
  
  variance_value <-
    selected_row$var_xhat[[1]] +
    2 *
    lambda_value *
    selected_row$cov_xhat_Jhat[[1]] +
    lambda_value^2 *
    selected_row$var_Jhat[[1]]
  
  sqrt(
    max(
      variance_value,
      1e-14
    )
  )
}


score_sd_at_critical <- setNames(
  vapply(
    channels,
    function(channel_name) {
      score_sd(
        channel_name = channel_name,
        lambda_value = lambda_star
      )
    },
    FUN.VALUE = numeric(1)
  ),
  channels
)


dispersion_table <- dispersion_table %>%
  mutate(
    Sigma_at_lambda_star =
      unname(
        score_sd_at_critical[
          channel
        ]
      ),
    
    channel_display =
      unname(
        channel_display_map[
          channel
        ]
      )
  )


write.csv(
  dispersion_table,
  file.path(
    data_dir,
    "exp2_dispersion.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "SINGLE-AGENT SCORE DISPERSION AT lambda* = 1\n"
)

cat(
  "============================================================\n\n"
)


print(
  dispersion_table %>%
    select(
      channel_display,
      mean_xhat,
      mean_Jhat,
      Sigma_at_lambda_star
    ) %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(
          .x,
          6
        )
      )
    ) %>%
    as.data.frame(),
  row.names = FALSE
)


# ================================================================
# 10. MONTE CARLO SETTINGS
# ================================================================

monte_carlo_replications <- if (FAST) {
  2000
} else {
  10000
}


population_grid <- c(
  25,
  50,
  100,
  200,
  500
)


lambda_grid <- seq(
  0,
  2,
  by = 0.05
)


# ================================================================
# 11. WILSON CONFIDENCE INTERVAL
# ================================================================

wilson_interval <- function(
    successes,
    trials,
    confidence_level = 0.95
) {
  if (
    trials <= 0 ||
    successes < 0 ||
    successes > trials
  ) {
    stop(
      "Invalid successes or trials in wilson_interval()."
    )
  }
  
  z_value <- qnorm(
    1 -
      (
        1 -
          confidence_level
      ) /
      2
  )
  
  p_hat <- successes /
    trials
  
  denominator <- 1 +
    z_value^2 /
    trials
  
  center <- (
    p_hat +
      z_value^2 /
      (
        2 *
          trials
      )
  ) /
    denominator
  
  half_width <- (
    z_value /
      denominator
  ) *
    sqrt(
      p_hat *
        (
          1 -
            p_hat
        ) /
        trials +
        z_value^2 /
        (
          4 *
            trials^2
        )
    )
  
  c(
    lower = max(
      0,
      center -
        half_width
    ),
    upper = min(
      1,
      center +
        half_width
    )
  )
}


# ================================================================
# 12. SIMULATE ONE CHANNEL-POPULATION CELL
# ================================================================

simulate_transition_cell <- function(
    channel_name,
    population_size,
    replications,
    lambda_values
) {
  type_matrix <- matrix(
    runif(
      replications *
        population_size,
      min = x_lower,
      max = x_upper
    ),
    nrow = replications,
    ncol = population_size
  )
  
  uniform_noise_matrix <- matrix(
    runif(
      replications *
        population_size
    ),
    nrow = replications,
    ncol = population_size
  )
  
  signal_matrix <- simulate_channel_from_uniforms(
    x = type_matrix,
    uniform_noise = uniform_noise_matrix,
    channel_name = channel_name
  )
  
  posterior_x_vector <- posterior_x_functions[[channel_name]](
    as.vector(
      signal_matrix
    )
  )
  
  posterior_J_vector <- posterior_J_functions[[channel_name]](
    as.vector(
      signal_matrix
    )
  )
  
  posterior_x_matrix <- matrix(
    posterior_x_vector,
    nrow = replications,
    ncol = population_size
  )
  
  posterior_J_matrix <- matrix(
    posterior_J_vector,
    nrow = replications,
    ncol = population_size
  )
  
  aggregate_xhat <- rowSums(
    posterior_x_matrix,
    na.rm = FALSE
  )
  
  aggregate_Jhat <- rowSums(
    posterior_J_matrix,
    na.rm = FALSE
  )
  
  result_rows <- lapply(
    lambda_values,
    function(lambda_value) {
      build_indicator <- (
        aggregate_xhat +
          lambda_value *
          aggregate_Jhat
      ) >= 0
      
      valid <- is.finite(
        aggregate_xhat +
          lambda_value *
          aggregate_Jhat
      )
      
      successes <- sum(
        build_indicator[valid]
      )
      
      valid_replications <- sum(
        valid
      )
      
      if (valid_replications == 0L) {
        return(
          data.frame(
            channel = channel_name,
            n = population_size,
            lambda = lambda_value,
            probability = NA_real_,
            ci_lower = NA_real_,
            ci_upper = NA_real_
          )
        )
      }
      
      probability <- successes /
        valid_replications
      
      confidence_interval <- wilson_interval(
        successes = successes,
        trials = valid_replications
      )
      
      data.frame(
        channel = channel_name,
        n = population_size,
        lambda = lambda_value,
        probability = probability,
        ci_lower = unname(
          confidence_interval["lower"]
        ),
        ci_upper = unname(
          confidence_interval["upper"]
        )
      )
    }
  )
  
  bind_rows(
    result_rows
  )
}


# ================================================================
# 13. RUN MONTE CARLO
# ================================================================

simulation_rows <- list()
simulation_counter <- 0


for (channel_name in channels) {
  cat(
    "\nSimulating channel:",
    channel_name,
    "\n"
  )
  
  for (population_size in population_grid) {
    cat(
      "  n =",
      population_size,
      "\n"
    )
    
    simulation_counter <- simulation_counter +
      1
    
    simulation_rows[[simulation_counter]] <-
      simulate_transition_cell(
        channel_name = channel_name,
        population_size = population_size,
        replications = monte_carlo_replications,
        lambda_values = lambda_grid
      )
  }
}


simulation_results <- bind_rows(
  simulation_rows
) %>%
  mutate(
    collapse_coordinate =
      sqrt(n) *
      (
        lambda -
          lambda_star
      ) /
      unname(
        score_sd_at_critical[
          channel
        ]
      ),
    
    channel_display = factor(
      unname(
        channel_display_map[
          channel
        ]
      ),
      levels = channel_display_names
    ),
    
    population_label = factor(
      paste0(
        "n = ",
        n
      ),
      levels = paste0(
        "n = ",
        population_grid
      )
    )
  )


write.csv(
  simulation_results,
  file.path(
    data_dir,
    "exp2_knife_edge.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 14. RAW CLT PREDICTION CURVES
# ================================================================

lambda_prediction_grid <- seq(
  0,
  2,
  length.out = 501
)


raw_prediction_rows <- lapply(
  channels,
  function(channel_name) {
    prediction_grid <- expand_grid(
      n = population_grid,
      lambda = lambda_prediction_grid
    )
    
    prediction_grid$channel <- channel_name
    
    prediction_grid$predicted_probability <- pnorm(
      sqrt(
        prediction_grid$n
      ) *
        mean_score(
          prediction_grid$lambda
        ) /
        vapply(
          prediction_grid$lambda,
          function(lambda_value) {
            score_sd(
              channel_name = channel_name,
              lambda_value = lambda_value
            )
          },
          FUN.VALUE = numeric(1)
        )
    )
    
    prediction_grid
  }
)


raw_prediction <- bind_rows(
  raw_prediction_rows
) %>%
  mutate(
    channel_display = factor(
      unname(
        channel_display_map[
          channel
        ]
      ),
      levels = channel_display_names
    ),
    
    population_label = factor(
      paste0(
        "n = ",
        n
      ),
      levels = paste0(
        "n = ",
        population_grid
      )
    )
  )


# ================================================================
# 15. LOCAL COLLAPSE REFERENCE CURVE
# ================================================================

collapse_u_grid <- seq(
  min(
    simulation_results$collapse_coordinate,
    na.rm = TRUE
  ),
  max(
    simulation_results$collapse_coordinate,
    na.rm = TRUE
  ),
  length.out = 601
)


collapse_reference <- data.frame(
  collapse_coordinate = collapse_u_grid,
  
  predicted_probability = pnorm(
    mean_virtual_value *
      collapse_u_grid
  )
)


# ================================================================
# 16. COLOR PALETTE
# ================================================================

population_colors <- c(
  "#0072B2",
  "#1B9E77",
  "#E6AB02",
  "#D55E00",
  "#7B3294"
)


population_color_map <- setNames(
  population_colors[
    seq_along(
      population_grid
    )
  ],
  paste0(
    "n = ",
    population_grid
  )
)


# ================================================================
# 17. FIGURE A: RAW TRANSITION
# ================================================================

raw_transition_plot <- ggplot() +
  geom_vline(
    xintercept = lambda_star,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  geom_line(
    data = raw_prediction,
    aes(
      x = lambda,
      y = predicted_probability,
      color = population_label,
      group = population_label
    ),
    linewidth = 0.85
  ) +
  geom_errorbar(
    data = simulation_results,
    aes(
      x = lambda,
      ymin = ci_lower,
      ymax = ci_upper,
      color = population_label
    ),
    width = 0,
    linewidth = 0.25,
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_point(
    data = simulation_results,
    aes(
      x = lambda,
      y = probability,
      color = population_label
    ),
    size = 1.45,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ channel_display,
    ncol = 2
  ) +
  scale_color_manual(
    values = population_color_map,
    name = "Agents"
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      2,
      by = 0.5
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    )
  ) +
  labs(
    x = expression(
      "revenue multiplier  " *
        lambda
    ),
    y = "implementation probability"
  ) +
  theme_paper()


if (SHOW_TITLES) {
  raw_transition_plot <- raw_transition_plot +
    labs(
      title = "Knife-edge transition in implementation",
      subtitle = paste0(
        "Points: Monte Carlo; curves: CLT approximation; ",
        "dashed line: lambda* = 1"
      )
    )
}


save_figure(
  plot_object = raw_transition_plot,
  filename_base = "fig_exp2_knife_edge_transition",
  width = 7.6,
  height = 6.0
)


# ================================================================
# 18. FIGURE B: LOCAL DATA COLLAPSE
# ================================================================

collapse_plot <- ggplot() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  geom_line(
    data = collapse_reference,
    aes(
      x = collapse_coordinate,
      y = predicted_probability
    ),
    color = "gray25",
    linewidth = 1.10
  ) +
  geom_point(
    data = simulation_results,
    aes(
      x = collapse_coordinate,
      y = probability,
      color = population_label
    ),
    size = 1.55,
    alpha = 0.90,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ channel_display,
    ncol = 2
  ) +
  scale_color_manual(
    values = population_color_map,
    name = "Agents"
  ) +
  coord_cartesian(
    xlim = c(
      -6,
      6
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    )
  ) +
  labs(
    x = expression(
      "rescaled multiplier  " *
        u ==
        sqrt(n) *
        (
          lambda -
            lambda^"*"
        ) /
        Sigma[Z](
          lambda^"*"
        )
    ),
    y = "implementation probability"
  ) +
  theme_paper()


if (SHOW_TITLES) {
  collapse_plot <- collapse_plot +
    labs(
      title = "Local collapse of the knife-edge transition",
      subtitle = expression(
        "Gray curve:  " *
          Phi(
            -0.5 *
              u
          )
      )
    )
}


save_figure(
  plot_object = collapse_plot,
  filename_base = "fig_exp2_knife_edge_collapse",
  width = 7.6,
  height = 6.0
)


# ================================================================
# 19. SANITY CHECKS
# ================================================================

at_critical <- simulation_results %>%
  filter(
    abs(
      lambda -
        lambda_star
    ) <
      1e-12
  ) %>%
  group_by(
    channel_display
  ) %>%
  summarize(
    probability_min = min(
      probability,
      na.rm = TRUE
    ),
    
    probability_max = max(
      probability,
      na.rm = TRUE
    ),
    
    probability_mean = mean(
      probability,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


monotonicity_check <- simulation_results %>%
  arrange(
    channel_display,
    n,
    lambda
  ) %>%
  group_by(
    channel_display,
    n
  ) %>%
  summarize(
    maximum_upward_step = {
      valid_probabilities <- probability[
        is.finite(
          probability
        )
      ]
      
      if (length(valid_probabilities) < 2L) {
        NA_real_
      } else {
        max(
          diff(
            valid_probabilities
          ),
          na.rm = TRUE
        )
      }
    },
    
    approximately_decreasing = {
      valid_probabilities <- probability[
        is.finite(
          probability
        )
      ]
      
      if (length(valid_probabilities) < 2L) {
        NA
      } else {
        max(
          diff(
            valid_probabilities
          ),
          na.rm = TRUE
        ) <= 0.03
      }
    },
    
    .groups = "drop"
  )


interpolate_lambda_at_probability <- function(
    lambda_values,
    probability_values,
    target_probability
) {
  valid <- is.finite(
    lambda_values
  ) &
    is.finite(
      probability_values
    )
  
  lambda_values <- lambda_values[
    valid
  ]
  
  probability_values <- probability_values[
    valid
  ]
  
  if (length(lambda_values) < 2L) {
    return(
      NA_real_
    )
  }
  
  ordering <- order(
    probability_values
  )
  
  probability_sorted <- probability_values[
    ordering
  ]
  
  lambda_sorted <- lambda_values[
    ordering
  ]
  
  keep <- !duplicated(
    probability_sorted
  )
  
  probability_sorted <- probability_sorted[
    keep
  ]
  
  lambda_sorted <- lambda_sorted[
    keep
  ]
  
  if (length(probability_sorted) < 2L) {
    return(
      NA_real_
    )
  }
  
  if (
    target_probability <
    min(
      probability_sorted
    ) ||
    target_probability >
    max(
      probability_sorted
    )
  ) {
    return(
      NA_real_
    )
  }
  
  approx(
    x = probability_sorted,
    y = lambda_sorted,
    xout = target_probability,
    rule = 1,
    ties = "ordered"
  )$y
}


# First calculate only one scalar crossing value per group.
transition_crossings <- simulation_results %>%
  group_by(
    channel_display,
    n
  ) %>%
  summarize(
    lambda_at_075 =
      interpolate_lambda_at_probability(
        lambda_values = lambda,
        probability_values = probability,
        target_probability = 0.75
      ),
    
    lambda_at_025 =
      interpolate_lambda_at_probability(
        lambda_values = lambda,
        probability_values = probability,
        target_probability = 0.25
      ),
    
    .groups = "drop"
  )


# Then calculate the width outside summarize().
transition_widths <- transition_crossings %>%
  mutate(
    transition_width =
      lambda_at_025 -
      lambda_at_075,
    
    scaled_width =
      sqrt(n) *
      transition_width
  )


sanity_table <- transition_widths %>%
  left_join(
    monotonicity_check,
    by = c(
      "channel_display",
      "n"
    )
  )


write.csv(
  sanity_table,
  file.path(
    data_dir,
    "exp2_transition_sanity.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT 2 SANITY REPORT\n"
)

cat(
  "============================================================\n\n"
)


cat(
  sprintf(
    "Critical multiplier lambda* = %.6f\n",
    lambda_star
  )
)


cat(
  sprintf(
    "Mean-score derivative at lambda* = E[J] = %.6f\n",
    mean_virtual_value
  )
)


cat(
  paste0(
    "At lambda*, implementation probabilities should be near 0.5\n",
    "when finite-sample skewness is small.\n\n"
  )
)


print(
  as.data.frame(
    at_critical
  ),
  digits = 5,
  row.names = FALSE
)


cat(
  "\nApproximate monotonicity in lambda:\n"
)


print(
  as.data.frame(
    monotonicity_check
  ),
  digits = 5,
  row.names = FALSE
)


cat(
  paste0(
    "\nEmpirical transition widths:\n",
    "The quantity sqrt(n) times the width should be approximately\n",
    "stable within each channel under the n^{-1/2} prediction.\n"
  )
)


print(
  as.data.frame(
    transition_widths
  ),
  digits = 5,
  row.names = FALSE
)


cat(
  sprintf(
    paste0(
      "\nCollapse reference at u = 0: Phi(0) = %.3f.\n",
      "The local reference curve is Phi(%.2f u).\n"
    ),
    pnorm(0),
    mean_virtual_value
  )
)


cat(
  paste0(
    "\nInterpretation:\n",
    "  * The raw curves should sharpen around lambda* = 1 as n grows.\n",
    "  * Within each channel, the rescaled points should approach the\n",
    "    gray curve Phi(-0.5 u).\n",
    "  * This is numerical evidence for the n^{-1/2} transition width\n",
    "    predicted by the central-limit approximation.\n"
  )
)


# ================================================================
# 20. COMPLETION REPORT
# ================================================================

cat(
  "\nFigures written to:\n"
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp2_knife_edge_transition.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp2_knife_edge_collapse.pdf"
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