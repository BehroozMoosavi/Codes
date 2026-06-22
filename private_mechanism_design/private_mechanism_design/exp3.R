# ================================================================
# EXPERIMENT 3
# Knife-edge sqrt(n) scaling of maximum reduced-form revenue
#
# Optimal Binary Mechanism under Locally Private Signals
#
# PRIOR
#
#     X ~ Uniform[0,1].
#
# For a Uniform[a,b] prior,
#
#     J(x) = 2x - b.
#
# Hence, in this experiment,
#
#     J(x) = 2x - 1,
#     E[J(X)] = 0.
#
# Let
#
#     Jhat(Y) = E[J(X) | Y].
#
# The revenue-maximizing reduced-form allocation is
#
#     q_rev(Y_1,...,Y_n)
#       = 1{
#           sum_i Jhat(Y_i) >= 0
#         },
#
# and therefore
#
#     R_n^{*,red}
#       = E[
#           (sum_i Jhat(Y_i))_+
#         ].
#
# At the knife edge E[Jhat(Y)] = 0, the theorem predicts
#
#     R_n^{*,red} / sqrt(n)
#       --> sigma_J / sqrt(2*pi),
#
# where
#
#     sigma_J^2 = Var(Jhat(Y)).
#
# The script produces two panels:
#
#   1. R_n^{*,red}/sqrt(n) against n, with the predicted
#      half-normal constant shown as a dashed line.
#
#   2. (R_n^{*,red})^2 against n, with the reference line
#
#          [sigma_J^2/(2*pi)] n.
#
# The second panel is a descriptive linear-axis visualization of
# the sqrt(n) rate. Formal inference is based on R_n itself, not on
# naively squaring a confidence interval.
#
# OUTPUT
#
#   figures/fig_exp3_sqrtn_revenue.pdf
#   figures/fig_exp3_rate_squared.pdf
#
#   tables/exp3_one_agent_moments.csv
#   tables/exp3_knife_edge_revenue.csv
#   tables/exp3_sanity_report.csv
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
  "R_codes/experiment/exp3"
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
    height = 5.4
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
# 3. PRIOR AND VIRTUAL VALUE
# ================================================================

x_lower <- 0
x_upper <- 1


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


virtual_value <- function(x) {
  2 * x -
    x_upper
}


# The knife-edge identity.
theoretical_mean_virtual_value <- x_lower


stopifnot(
  abs(
    theoretical_mean_virtual_value
  ) <
    1e-12
)


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


channel_colors <- c(
  "Gaussian" = "#0072B2",
  "Laplace" = "#D55E00",
  "Logistic" = "#1B9E77"
)


# ================================================================
# 5. CONDITIONAL LOG-DENSITIES
#
# Log-densities are used for stable posterior quadrature.
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
# 6. CHANNEL SIMULATION USING UNIFORM RANDOM NUMBERS
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
  
  virtual_values <- virtual_value(
    type_grid
  )
  
  log_kernel <- get_log_kernel(
    channel_name
  )
  
  rows <- lapply(
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
    rows
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
# 8. NUMERICAL GRIDS AND LOOKUP TABLES
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
      cat(
        "Building posterior table for",
        channel_name,
        "...\n"
      )
      
      posterior_moment_table(
        signal_grid = signal_grid,
        channel_name = channel_name,
        type_grid = type_grid
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
# 9. ONE-AGENT MOMENTS AND THEORETICAL CONSTANTS
# ================================================================

one_agent_draws <- if (FAST) {
  40000
} else {
  200000
}


# Common type and uniform-noise draws improve comparisons across
# channels, although the resulting signals have different marginals.

one_agent_types <- runif(
  one_agent_draws,
  min = x_lower,
  max = x_upper
)


one_agent_uniforms <- runif(
  one_agent_draws
)


one_agent_rows <- lapply(
  channels,
  function(channel_name) {
    signals <- simulate_channel_from_uniforms(
      x = one_agent_types,
      uniform_noise = one_agent_uniforms,
      channel_name = channel_name
    )
    
    posterior_J <- posterior_J_functions[[channel_name]](
      signals
    )
    
    posterior_J <- posterior_J[
      is.finite(
        posterior_J
      )
    ]
    
    sample_mean_Jhat <- mean(
      posterior_J
    )
    
    sample_sd_Jhat <- sd(
      posterior_J
    )
    
    half_normal_constant <-
      sample_sd_Jhat /
      sqrt(
        2 *
          pi
      )
    
    data.frame(
      channel = channel_name,
      mean_Jhat = sample_mean_Jhat,
      sd_Jhat = sample_sd_Jhat,
      variance_Jhat = sample_sd_Jhat^2,
      half_normal_constant = half_normal_constant,
      squared_rate_constant = half_normal_constant^2,
      stringsAsFactors = FALSE
    )
  }
)


one_agent_table <- bind_rows(
  one_agent_rows
) %>%
  mutate(
    channel_display = factor(
      unname(
        channel_display_map[
          channel
        ]
      ),
      levels = channel_display_names
    )
  )


write.csv(
  one_agent_table,
  file.path(
    data_dir,
    "exp3_one_agent_moments.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "ONE-AGENT POSTERIOR VIRTUAL-VALUE MOMENTS\n"
)

cat(
  "============================================================\n\n"
)


print(
  one_agent_table %>%
    select(
      channel_display,
      mean_Jhat,
      sd_Jhat,
      half_normal_constant
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

population_grid <- c(
  25,
  50,
  100,
  200,
  500,
  1000,
  2000
)


monte_carlo_replications <- if (FAST) {
  3000
} else {
  10000
}


# Computations are performed in batches to avoid allocating a very
# large B-by-n matrix when n is large.

batch_size <- if (FAST) {
  300
} else {
  500
}


# ================================================================
# 11. SIMULATE MAXIMUM REDUCED-FORM REVENUE
#
# The simulated replication-level quantity is exactly
#
#     (sum_i Jhat(Y_i))_+.
#
# No welfare indicator is included.
# ================================================================

simulate_revenue_cell <- function(
    channel_name,
    population_size,
    replications,
    batch_size
) {
  revenue_draws <- numeric(
    replications
  )
  
  completed <- 0L
  
  while (completed < replications) {
    current_batch <- min(
      batch_size,
      replications -
        completed
    )
    
    type_matrix <- matrix(
      runif(
        current_batch *
          population_size,
        min = x_lower,
        max = x_upper
      ),
      nrow = current_batch,
      ncol = population_size
    )
    
    uniform_matrix <- matrix(
      runif(
        current_batch *
          population_size
      ),
      nrow = current_batch,
      ncol = population_size
    )
    
    signal_matrix <- simulate_channel_from_uniforms(
      x = type_matrix,
      uniform_noise = uniform_matrix,
      channel_name = channel_name
    )
    
    posterior_J_vector <- posterior_J_functions[[channel_name]](
      as.vector(
        signal_matrix
      )
    )
    
    posterior_J_matrix <- matrix(
      posterior_J_vector,
      nrow = current_batch,
      ncol = population_size
    )
    
    aggregate_virtual_value <- rowSums(
      posterior_J_matrix,
      na.rm = FALSE
    )
    
    if (any(!is.finite(aggregate_virtual_value))) {
      stop(
        paste0(
          "Non-finite aggregate posterior virtual values for channel ",
          channel_name,
          " at n = ",
          population_size,
          ". Increase the signal-grid range or inspect posterior tables."
        )
      )
    }
    
    batch_revenue <- pmax(
      aggregate_virtual_value,
      0
    )
    
    indices <- (
      completed +
        1L
    ):(
      completed +
        current_batch
    )
    
    revenue_draws[
      indices
    ] <- batch_revenue
    
    completed <- completed +
      current_batch
  }
  
  revenue_mean <- mean(
    revenue_draws
  )
  
  revenue_standard_error <- sd(
    revenue_draws
  ) /
    sqrt(
      replications
    )
  
  data.frame(
    channel = channel_name,
    n = population_size,
    replications = replications,
    revenue = revenue_mean,
    revenue_se = revenue_standard_error,
    stringsAsFactors = FALSE
  )
}


# ================================================================
# 12. RUN THE MONTE CARLO EXPERIMENT
# ================================================================

simulation_rows <- list()
simulation_counter <- 0L


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
      1L
    
    simulation_rows[[simulation_counter]] <-
      simulate_revenue_cell(
        channel_name = channel_name,
        population_size = population_size,
        replications = monte_carlo_replications,
        batch_size = batch_size
      )
  }
}


simulation_results <- bind_rows(
  simulation_rows
) %>%
  mutate(
    channel_display = factor(
      unname(
        channel_display_map[
          channel
        ]
      ),
      levels = channel_display_names
    )
  ) %>%
  left_join(
    one_agent_table %>%
      select(
        channel,
        half_normal_constant,
        squared_rate_constant
      ),
    by = "channel"
  ) %>%
  mutate(
    revenue_ci_lower =
      pmax(
        0,
        revenue -
          1.96 *
          revenue_se
      ),
    
    revenue_ci_upper =
      revenue +
      1.96 *
      revenue_se,
    
    revenue_over_sqrt_n =
      revenue /
      sqrt(n),
    
    revenue_over_sqrt_n_se =
      revenue_se /
      sqrt(n),
    
    normalized_ci_lower =
      pmax(
        0,
        revenue_over_sqrt_n -
          1.96 *
          revenue_over_sqrt_n_se
      ),
    
    normalized_ci_upper =
      revenue_over_sqrt_n +
      1.96 *
      revenue_over_sqrt_n_se,
    
    revenue_squared =
      revenue^2,
    
    # Delta-method standard error for g(r)=r^2:
    #
    #     se[g(rhat)] approximately |2 rhat| se(rhat).
    revenue_squared_se =
      2 *
      abs(revenue) *
      revenue_se,
    
    revenue_squared_ci_lower =
      pmax(
        0,
        revenue_squared -
          1.96 *
          revenue_squared_se
      ),
    
    revenue_squared_ci_upper =
      revenue_squared +
      1.96 *
      revenue_squared_se,
    
    squared_reference =
      squared_rate_constant *
      n
  )


write.csv(
  simulation_results,
  file.path(
    data_dir,
    "exp3_knife_edge_revenue.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 13. FIGURE 1
# R_n / sqrt(n) APPROACHES THE HALF-NORMAL CONSTANT
# ================================================================

constant_plot <- ggplot(
  simulation_results,
  aes(
    x = n,
    y = revenue_over_sqrt_n,
    color = channel_display,
    group = channel_display
  )
) +
  geom_ribbon(
    aes(
      ymin = normalized_ci_lower,
      ymax = normalized_ci_upper,
      fill = channel_display
    ),
    alpha = 0.10,
    color = NA,
    show.legend = FALSE
  ) +
  geom_hline(
    data = one_agent_table,
    aes(
      yintercept = half_normal_constant,
      color = channel_display
    ),
    linetype = "dashed",
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    size = 2.30
  ) +
  scale_color_manual(
    values = channel_colors
  ) +
  scale_fill_manual(
    values = channel_colors
  ) +
  scale_x_continuous(
    breaks = population_grid
  ) +
  labs(
    x = "number of agents  n",
    y = expression(
      R[n]^{"*,red"} /
        sqrt(n)
    ),
    subtitle = expression(
      "Dashed lines:  " *
        sigma[J] /
        sqrt(
          2 *
            pi
        )
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  constant_plot <- constant_plot +
    labs(
      title = expression(
        "Knife-edge " *
          sqrt(n) *
          " revenue scaling"
      )
    )
}


save_figure(
  plot_object = constant_plot,
  filename_base = "fig_exp3_sqrtn_revenue",
  width = 7.6,
  height = 5.4
)


# ================================================================
# 14. FIGURE 2
# SQUARED REVENUE AGAINST n
# ================================================================

squared_plot <- ggplot(
  simulation_results,
  aes(
    x = n,
    y = revenue_squared,
    color = channel_display,
    group = channel_display
  )
) +
  geom_line(
    aes(
      y = squared_reference
    ),
    linetype = "dashed",
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(
      ymin = revenue_squared_ci_lower,
      ymax = revenue_squared_ci_upper
    ),
    width = 20,
    linewidth = 0.35
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    size = 2.30
  ) +
  scale_color_manual(
    values = channel_colors
  ) +
  scale_x_continuous(
    breaks = population_grid
  ) +
  labs(
    x = "number of agents  n",
    y = expression(
      (
        R[n]^{"*,red"}
      )^2
    ),
    subtitle = expression(
      "Dashed lines:  " *
        (
          sigma[J]^2 /
            (
              2 *
                pi
            )
        ) *
        n
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  squared_plot <- squared_plot +
    labs(
      title = "Linear-axis visualization of the knife-edge rate"
    )
}


save_figure(
  plot_object = squared_plot,
  filename_base = "fig_exp3_rate_squared",
  width = 7.6,
  height = 5.4
)


# ================================================================
# 15. SANITY REPORT
# ================================================================

largest_n <- max(
  population_grid
)


constant_check <- simulation_results %>%
  filter(
    n ==
      largest_n
  ) %>%
  transmute(
    channel_display,
    n,
    revenue_over_sqrt_n,
    predicted_constant =
      half_normal_constant,
    empirical_to_predicted =
      revenue_over_sqrt_n /
      predicted_constant
  )


slope_rows <- lapply(
  levels(
    simulation_results$channel_display
  ),
  function(channel_label) {
    selected_data <- simulation_results[
      simulation_results$channel_display ==
        channel_label,
      ,
      drop = FALSE
    ]
    
    fitted_model <- lm(
      revenue_squared ~
        0 +
        n,
      data = selected_data
    )
    
    fitted_slope <- unname(
      coef(
        fitted_model
      )[["n"]]
    )
    
    predicted_slope <- unique(
      selected_data$squared_rate_constant
    )
    
    if (length(predicted_slope) != 1L) {
      stop(
        paste(
          "Expected one predicted slope for",
          channel_label
        )
      )
    }
    
    data.frame(
      channel_display = channel_label,
      fitted_slope = fitted_slope,
      predicted_slope = predicted_slope,
      fitted_to_predicted =
        fitted_slope /
        predicted_slope,
      stringsAsFactors = FALSE
    )
  }
)


slope_check <- bind_rows(
  slope_rows
)


sanity_report <- constant_check %>%
  left_join(
    slope_check,
    by = "channel_display"
  )


write.csv(
  sanity_report,
  file.path(
    data_dir,
    "exp3_sanity_report.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT 3 SANITY REPORT\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Theoretical knife-edge identity:\n"
)

cat(
  sprintf(
    "  E[J(X)] = %.6f\n\n",
    theoretical_mean_virtual_value
  )
)


cat(
  "Estimated one-agent constants:\n"
)


print(
  one_agent_table %>%
    select(
      channel_display,
      mean_Jhat,
      sd_Jhat,
      half_normal_constant
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


cat(
  "\nAt the largest population size, R_n/sqrt(n) versus prediction:\n"
)


print(
  as.data.frame(
    constant_check
  ),
  digits = 5,
  row.names = FALSE
)


cat(
  "\nFitted slope of R_n^2 on n through the origin:\n"
)


print(
  as.data.frame(
    slope_check
  ),
  digits = 5,
  row.names = FALSE
)


cat(
  paste0(
    "\nInterpretation:\n",
    "  * mean_Jhat should be close to zero for every channel.\n",
    "  * empirical_to_predicted should approach one as n increases.\n",
    "  * fitted_to_predicted should be reasonably close to one.\n",
    "  * The second panel is a descriptive visualization of the rate;\n",
    "    the first panel corresponds directly to the asymptotic theorem.\n"
  )
)


# ================================================================
# 16. COMPLETION REPORT
# ================================================================

cat(
  "\nFigures written to:\n"
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp3_sqrtn_revenue.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp3_rate_squared.pdf"
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