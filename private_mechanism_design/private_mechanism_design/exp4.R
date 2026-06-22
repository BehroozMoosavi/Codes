# ================================================================
# EXPERIMENT 4
# Large-deviation rate and its numerical dependence on Gaussian noise
#
# Optimal Binary Mechanism under Locally Private Signals
#
# PRIOR AND RULE
#
#     X ~ Uniform[-0.5,1.5],
#     lambda = 1.5,
#     Y = X + N(0,sigma^2).
#
# For a Uniform[a,b] prior,
#
#     J(x) = 2x - b.
#
# Here,
#
#     E[X] = 0.5,
#     E[J(X)] = a = -0.5,
#
# so the one-agent posterior score
#
#     S_sigma(Y)
#       = E[X | Y]
#         + lambda E[J(X) | Y]
#
# has mean
#
#     E[S_sigma(Y)]
#       = 0.5 + 1.5(-0.5)
#       = -0.25.
#
# Thus implementation,
#
#     sum_i S_sigma(Y_i) >= 0,
#
# is a right-tail large-deviation event.
#
# For each sigma, the Cramer rate at the build boundary is
#
#     I_sigma(0)
#       = sup_{t >= 0}
#           {
#             -log E[exp(t S_sigma(Y))]
#           }.
#
# This script estimates I_sigma(0) directly from a large common-random-
# numbers simulation of the one-agent posterior score. It then separately
# simulates finite-n implementation probabilities and compares
#
#     -log P(build) / n
#
# with the estimated rate.
#
# IMPORTANT INTERPRETATION
#
# The final regression
#
#     I_hat(sigma)
#       approximately
#       beta_0 + beta_1 sigma^2
#
# is descriptive over the displayed sigma range. The script does not claim
#
#     I_sigma(0) = Theta(sigma^2)
#
# as sigma -> 0. Since the zero-noise rate is generally positive here, the
# more natural local question is whether
#
#     I_sigma(0) - I_0(0)
#
# is approximately proportional to sigma^2 over a selected range.
#
# OUTPUT
#
#   figures/fig_exp4_implementation_vs_sigma.pdf
#   figures/fig_exp4_logaccept_vs_n.pdf
#   figures/fig_exp4_rate_vs_sigma2.pdf
#   figures/fig_exp4_finite_n_rate.pdf
#
#   tables/exp4_one_agent_rate.csv
#   tables/exp4_implementation_probabilities.csv
#   tables/exp4_finite_n_rate.csv
#   tables/exp4_sanity_report.csv
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
  "R_codes/experiment/exp4"
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
  "tibble",
  "scales"
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
library(scales)


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
    width = 7.5,
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
# 3. PRIOR, VIRTUAL VALUE, AND SCORE
# ================================================================

x_lower <- -0.5
x_upper <-  1.5

lambda_value <- 1.5


prior_mean <- (
  x_lower +
    x_upper
) / 2


mean_virtual_value <- x_lower


mean_score_theory <-
  prior_mean +
  lambda_value *
  mean_virtual_value


stopifnot(
  abs(
    mean_score_theory +
      0.25
  ) <
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


virtual_value <- function(x) {
  2 * x -
    x_upper
}


# ================================================================
# 4. NOISE AND SAMPLE-SIZE GRIDS
# ================================================================

sigma_grid <- c(
  0.20,
  0.30,
  0.40,
  0.50,
  0.60,
  0.80,
  1.00,
  1.20
)


population_grid <- c(
  20,
  40,
  60,
  80,
  120,
  160
)


one_agent_draws <- if (FAST) {
  100000
} else {
  600000
}


implementation_replications <- if (FAST) {
  5000
} else {
  30000
}


batch_size <- if (FAST) {
  500
} else {
  1000
}


# ================================================================
# 5. STABLE GAUSSIAN POSTERIOR MOMENTS
# ================================================================

type_grid <- seq(
  x_lower,
  x_upper,
  length.out = if (FAST) {
    2001
  } else {
    5001
  }
)


type_step <- type_grid[2] -
  type_grid[1]


prior_values <- uniform_density(
  type_grid,
  x_lower,
  x_upper
)


virtual_values <- virtual_value(
  type_grid
)


# Signal grids are sigma-specific. The support of X is bounded, so
# extending the grid by 9 standard deviations gives negligible Gaussian
# tail truncation for the displayed range.

make_signal_grid <- function(
    sigma,
    number_of_points
) {
  seq(
    x_lower -
      9 *
      sigma,
    x_upper +
      9 *
      sigma,
    length.out = number_of_points
  )
}


posterior_table_gaussian <- function(
    sigma,
    signal_grid
) {
  rows <- lapply(
    signal_grid,
    function(signal_value) {
      log_weights <- dnorm(
        signal_value,
        mean = type_grid,
        sd = sigma,
        log = TRUE
      )
      
      positive_prior <- prior_values > 0
      
      log_weights[positive_prior] <-
        log_weights[positive_prior] +
        log(
          prior_values[positive_prior]
        )
      
      log_weights[!positive_prior] <- -Inf
      
      maximum_log_weight <- max(
        log_weights
      )
      
      stabilized_weights <- exp(
        log_weights -
          maximum_log_weight
      )
      
      denominator <- sum(
        stabilized_weights
      ) *
        type_step
      
      if (
        !is.finite(denominator) ||
        denominator <= .Machine$double.eps
      ) {
        return(
          data.frame(
            y = signal_value,
            xhat = NA_real_,
            Jhat = NA_real_,
            score = NA_real_
          )
        )
      }
      
      posterior_x <- sum(
        type_grid *
          stabilized_weights
      ) *
        type_step /
        denominator
      
      posterior_J <- sum(
        virtual_values *
          stabilized_weights
      ) *
        type_step /
        denominator
      
      posterior_score <-
        posterior_x +
        lambda_value *
        posterior_J
      
      data.frame(
        y = signal_value,
        xhat = posterior_x,
        Jhat = posterior_J,
        score = posterior_score
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


posterior_objects <- setNames(
  lapply(
    sigma_grid,
    function(sigma_value) {
      cat(
        "Building posterior table for sigma =",
        sigma_value,
        "...\n"
      )
      
      signal_grid <- make_signal_grid(
        sigma = sigma_value,
        number_of_points = if (FAST) {
          3001
        } else {
          7001
        }
      )
      
      posterior_table <- posterior_table_gaussian(
        sigma = sigma_value,
        signal_grid = signal_grid
      )
      
      list(
        table = posterior_table,
        score = make_interpolator(
          posterior_table,
          "score"
        ),
        Jhat = make_interpolator(
          posterior_table,
          "Jhat"
        )
      )
    }
  ),
  as.character(
    sigma_grid
  )
)


# ================================================================
# 6. LOG-MEAN-EXP AND CRAMER-RATE ESTIMATION
# ================================================================

log_mean_exp <- function(values) {
  maximum_value <- max(
    values
  )
  
  maximum_value +
    log(
      mean(
        exp(
          values -
            maximum_value
        )
      )
    )
}


estimate_cramer_rate <- function(
    score_draws,
    initial_upper = 5,
    maximum_upper = 100
) {
  score_draws <- score_draws[
    is.finite(
      score_draws
    )
  ]
  
  if (length(score_draws) < 100L) {
    stop(
      "Too few finite score draws in estimate_cramer_rate()."
    )
  }
  
  if (mean(score_draws) >= 0) {
    stop(
      "The one-agent score mean must be negative for this experiment."
    )
  }
  
  objective <- function(t_value) {
    log_mean_exp(
      t_value *
        score_draws
    )
  }
  
  upper <- initial_upper
  
  repeat {
    optimization <- optimize(
      f = objective,
      interval = c(
        0,
        upper
      ),
      maximum = FALSE,
      tol = 1e-10
    )
    
    # If the minimizer is not close to the upper boundary, accept it.
    if (
      optimization$minimum <
      0.95 *
      upper
    ) {
      break
    }
    
    upper <- 2 *
      upper
    
    if (upper > maximum_upper) {
      stop(
        paste0(
          "Unable to bracket the tilting parameter. ",
          "Increase maximum_upper."
        )
      )
    }
  }
  
  t_hat <- optimization$minimum
  
  log_mgf_hat <- optimization$objective
  
  rate_hat <- -log_mgf_hat
  
  exponential_terms <- exp(
    t_hat *
      score_draws -
      max(
        t_hat *
          score_draws
      )
  )
  
  # Delta-method standard error for -log(mean(exp(t*S))).
  #
  # The scale normalization cancels in sd(Z)/mean(Z).
  
  relative_standard_error <-
    sd(
      exponential_terms
    ) /
    sqrt(
      length(
        exponential_terms
      )
    ) /
    mean(
      exponential_terms
    )
  
  rate_standard_error <- relative_standard_error
  
  c(
    t_hat = t_hat,
    rate_hat = rate_hat,
    rate_se = rate_standard_error,
    mean_score = mean(
      score_draws
    ),
    sd_score = sd(
      score_draws
    )
  )
}


# ================================================================
# 7. ONE-AGENT SCORE DRAWS AND DIRECT RATE ESTIMATES
#
# Common X and standard-normal draws are reused across sigma values.
# This reduces noise in the cross-sigma comparison.
# ================================================================

one_agent_types <- runif(
  one_agent_draws,
  min = x_lower,
  max = x_upper
)


one_agent_standard_normals <- rnorm(
  one_agent_draws
)


rate_rows <- lapply(
  sigma_grid,
  function(sigma_value) {
    cat(
      "Estimating one-agent rate for sigma =",
      sigma_value,
      "...\n"
    )
    
    signals <-
      one_agent_types +
      sigma_value *
      one_agent_standard_normals
    
    score_draws <-
      posterior_objects[[
        as.character(
          sigma_value
        )
      ]]$score(
        signals
      )
    
    rate_summary <- estimate_cramer_rate(
      score_draws = score_draws
    )
    
    data.frame(
      sigma = sigma_value,
      sigma2 = sigma_value^2,
      t_hat = unname(
        rate_summary["t_hat"]
      ),
      I_hat = unname(
        rate_summary["rate_hat"]
      ),
      I_se = unname(
        rate_summary["rate_se"]
      ),
      mean_score = unname(
        rate_summary["mean_score"]
      ),
      sd_score = unname(
        rate_summary["sd_score"]
      ),
      number_of_draws = length(
        score_draws
      ),
      stringsAsFactors = FALSE
    )
  }
)


rate_table <- bind_rows(
  rate_rows
) %>%
  mutate(
    I_ci_lower =
      pmax(
        0,
        I_hat -
          1.96 *
          I_se
      ),
    
    I_ci_upper =
      I_hat +
      1.96 *
      I_se
  )


write.csv(
  rate_table,
  file.path(
    data_dir,
    "exp4_one_agent_rate.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 8. FINITE-n IMPLEMENTATION PROBABILITIES
# ================================================================

wilson_interval <- function(
    successes,
    trials,
    confidence_level = 0.95
) {
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
  
  denominator <-
    1 +
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


simulate_probability_cell <- function(
    sigma_value,
    population_size,
    replications,
    batch_size
) {
  score_function <-
    posterior_objects[[
      as.character(
        sigma_value
      )
    ]]$score
  
  build_indicators <- logical(
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
    
    noise_matrix <- matrix(
      rnorm(
        current_batch *
          population_size,
        mean = 0,
        sd = sigma_value
      ),
      nrow = current_batch,
      ncol = population_size
    )
    
    signal_matrix <-
      type_matrix +
      noise_matrix
    
    score_vector <- score_function(
      as.vector(
        signal_matrix
      )
    )
    
    score_matrix <- matrix(
      score_vector,
      nrow = current_batch,
      ncol = population_size
    )
    
    aggregate_score <- rowSums(
      score_matrix
    )
    
    if (any(!is.finite(aggregate_score))) {
      stop(
        paste0(
          "Non-finite aggregate score for sigma = ",
          sigma_value,
          ", n = ",
          population_size,
          "."
        )
      )
    }
    
    indices <- (
      completed +
        1L
    ):(
      completed +
        current_batch
    )
    
    build_indicators[
      indices
    ] <-
      aggregate_score >= 0
    
    completed <- completed +
      current_batch
  }
  
  success_count <- sum(
    build_indicators
  )
  
  # Smoothed estimate is used only for logs and finite-n rate plots.
  # The ordinary Monte Carlo estimate remains success_count/B.
  
  probability_hat <-
    success_count /
    replications
  
  probability_smoothed <-
    (
      success_count +
        0.5
    ) /
    (
      replications +
        1
    )
  
  confidence_interval <- wilson_interval(
    successes = success_count,
    trials = replications
  )
  
  data.frame(
    n = population_size,
    sigma = sigma_value,
    replications = replications,
    accept_count = success_count,
    acceptance = probability_hat,
    acceptance_smoothed = probability_smoothed,
    acc_low = unname(
      confidence_interval["lower"]
    ),
    acc_high = unname(
      confidence_interval["upper"]
    ),
    stringsAsFactors = FALSE
  )
}


simulation_grid <- expand_grid(
  sigma = sigma_grid,
  n = population_grid
)


simulation_rows <- vector(
  mode = "list",
  length = nrow(
    simulation_grid
  )
)


for (index in seq_len(
  nrow(
    simulation_grid
  )
)) {
  sigma_value <- simulation_grid$sigma[index]
  population_size <- simulation_grid$n[index]
  
  cat(
    "Simulating sigma =",
    sigma_value,
    ", n =",
    population_size,
    "...\n"
  )
  
  simulation_rows[[index]] <-
    simulate_probability_cell(
      sigma_value = sigma_value,
      population_size = population_size,
      replications = implementation_replications,
      batch_size = batch_size
    )
}


probability_table <- bind_rows(
  simulation_rows
) %>%
  left_join(
    rate_table %>%
      select(
        sigma,
        I_hat
      ),
    by = "sigma"
  ) %>%
  mutate(
    finite_n_rate =
      -log(
        acceptance_smoothed
      ) /
      n,
    
    predicted_probability =
      exp(
        -n *
          I_hat
      )
  )


write.csv(
  probability_table,
  file.path(
    data_dir,
    "exp4_implementation_probabilities.csv"
  ),
  row.names = FALSE
)


write.csv(
  probability_table %>%
    select(
      sigma,
      n,
      finite_n_rate,
      I_hat,
      acceptance,
      acceptance_smoothed,
      accept_count
    ),
  file.path(
    data_dir,
    "exp4_finite_n_rate.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 9. COLOR PALETTES
# ================================================================

sigma_colors <- setNames(
  colorRampPalette(
    c(
      "#0072B2",
      "#1B9E77",
      "#E6AB02",
      "#D55E00"
    )
  )(
    length(
      sigma_grid
    )
  ),
  as.character(
    sigma_grid
  )
)


population_colors <- setNames(
  colorRampPalette(
    c(
      "#0072B2",
      "#D55E00"
    )
  )(
    length(
      population_grid
    )
  ),
  as.character(
    population_grid
  )
)


# ================================================================
# 10. FIGURE 1
# IMPLEMENTATION PROBABILITY AGAINST sigma
# ================================================================

implementation_plot <- ggplot(
  probability_table,
  aes(
    x = sigma,
    y = acceptance,
    color = factor(n),
    group = factor(n)
  )
) +
  geom_ribbon(
    aes(
      ymin = acc_low,
      ymax = acc_high,
      fill = factor(n)
    ),
    alpha = 0.08,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 0.95
  ) +
  geom_point(
    size = 2.1
  ) +
  scale_color_manual(
    values = population_colors,
    name = "Agents"
  ) +
  scale_fill_manual(
    values = population_colors
  ) +
  scale_x_continuous(
    breaks = sigma_grid
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  labs(
    x = expression(
      "Gaussian noise level  " *
        sigma
    ),
    y = "implementation probability"
  ) +
  theme_paper()


if (SHOW_TITLES) {
  implementation_plot <- implementation_plot +
    labs(
      title = "Implementation probability and Gaussian noise"
    )
}


save_figure(
  plot_object = implementation_plot,
  filename_base = "fig_exp4_implementation_vs_sigma",
  width = 7.5,
  height = 5.4
)


# ================================================================
# 11. FIGURE 2
# EXPONENTIAL DECAY AGAINST n
#
# The dashed curves show exp(-n I_hat). They capture the exponential
# rate, not the generally nontrivial subexponential prefactor.
# ================================================================

decay_plot <- ggplot(
  probability_table,
  aes(
    x = n,
    y = acceptance_smoothed,
    color = factor(sigma),
    group = factor(sigma)
  )
) +
  geom_line(
    aes(
      y = predicted_probability
    ),
    linetype = "dashed",
    linewidth = 0.75
  ) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.65
  ) +
  geom_point(
    size = 2.0
  ) +
  scale_y_log10(
    labels = label_log()
  ) +
  scale_color_manual(
    values = sigma_colors,
    name = expression(
      sigma
    )
  ) +
  scale_x_continuous(
    breaks = population_grid
  ) +
  labs(
    x = "number of agents  n",
    y = "implementation probability (log scale)",
    subtitle = expression(
      "Dashed curves:  " *
        exp(
          -n *
            hat(I)[sigma](
              0
            )
        )
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  decay_plot <- decay_plot +
    labs(
      title = "Finite-sample implementation probabilities"
    )
}


save_figure(
  plot_object = decay_plot,
  filename_base = "fig_exp4_logaccept_vs_n",
  width = 7.5,
  height = 5.4
)


# ================================================================
# 12. FIGURE 3
# DIRECT RATE ESTIMATE AGAINST sigma^2
# ================================================================

descriptive_fit <- lm(
  I_hat ~ sigma2,
  data = rate_table
)


fit_intercept <- unname(
  coef(
    descriptive_fit
  )[["(Intercept)"]]
)


fit_slope <- unname(
  coef(
    descriptive_fit
  )[["sigma2"]]
)


fit_r_squared <- summary(
  descriptive_fit
)$r.squared


rate_plot <- ggplot(
  rate_table,
  aes(
    x = sigma2,
    y = I_hat
  )
) +
  geom_abline(
    intercept = fit_intercept,
    slope = fit_slope,
    linetype = "dashed",
    color = "black",
    linewidth = 0.75
  ) +
  geom_errorbar(
    aes(
      ymin = I_ci_lower,
      ymax = I_ci_upper
    ),
    width = 0.015,
    linewidth = 0.40,
    color = "gray35"
  ) +
  geom_point(
    size = 3.0,
    color = "#0072B2"
  ) +
  annotate(
    geom = "text",
    x = min(
      rate_table$sigma2
    ),
    y = max(
      rate_table$I_ci_upper
    ),
    label = paste0(
      "descriptive fit:  I = ",
      round(
        fit_intercept,
        4
      ),
      " + ",
      round(
        fit_slope,
        4
      ),
      " sigma^2",
      "\nR^2 = ",
      round(
        fit_r_squared,
        3
      )
    ),
    hjust = 0,
    vjust = 1,
    family = "serif",
    size = 4.1
  ) +
  labs(
    x = expression(
      sigma^2
    ),
    y = expression(
      "estimated rate  " *
        hat(I)[sigma](
          0
        )
    ),
    subtitle = paste0(
      "Affine relation shown only as a descriptive fit ",
      "over the displayed range"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  rate_plot <- rate_plot +
    labs(
      title = "Large-deviation rate and Gaussian noise"
    )
}


save_figure(
  plot_object = rate_plot,
  filename_base = "fig_exp4_rate_vs_sigma2",
  width = 7.5,
  height = 5.4
)


# ================================================================
# 13. FIGURE 4
# FINITE-n RATE -log P_n / n AGAINST n
# ================================================================

finite_rate_plot <- ggplot(
  probability_table,
  aes(
    x = n,
    y = finite_n_rate,
    color = factor(sigma),
    group = factor(sigma)
  )
) +
  geom_hline(
    data = rate_table,
    aes(
      yintercept = I_hat,
      color = factor(sigma)
    ),
    linetype = "dashed",
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 0.90
  ) +
  geom_point(
    size = 2.1
  ) +
  scale_color_manual(
    values = sigma_colors,
    name = expression(
      sigma
    )
  ) +
  scale_x_continuous(
    breaks = population_grid
  ) +
  labs(
    x = "number of agents  n",
    y = expression(
      -log(
        P(
          q^"*" ==
            1
        )
      ) /
        n
    ),
    subtitle = "Dashed lines are direct one-agent Cramer-rate estimates"
  ) +
  theme_paper()


if (SHOW_TITLES) {
  finite_rate_plot <- finite_rate_plot +
    labs(
      title = "Finite-sample convergence toward the Cramer rate"
    )
}


save_figure(
  plot_object = finite_rate_plot,
  filename_base = "fig_exp4_finite_n_rate",
  width = 7.5,
  height = 5.4
)


# ================================================================
# 14. SANITY CHECKS
# ================================================================

mean_check <- rate_table %>%
  transmute(
    sigma,
    estimated_mean_score = mean_score,
    theoretical_mean_score = mean_score_theory,
    mean_error =
      estimated_mean_score -
      theoretical_mean_score
  )


largest_n <- max(
  population_grid
)


finite_rate_check <- probability_table %>%
  filter(
    n ==
      largest_n
  ) %>%
  transmute(
    sigma,
    n,
    acceptance,
    accept_count,
    finite_n_rate,
    direct_rate = I_hat,
    finite_to_direct =
      finite_n_rate /
      direct_rate
  )


monotonicity_check <- rate_table %>%
  arrange(
    sigma
  ) %>%
  summarize(
    minimum_rate_increment = min(
      diff(
        I_hat
      )
    ),
    rate_increasing_over_grid =
      minimum_rate_increment >
      -0.005
  )


sanity_report <- rate_table %>%
  select(
    sigma,
    sigma2,
    mean_score,
    sd_score,
    t_hat,
    I_hat,
    I_se
  )


write.csv(
  sanity_report,
  file.path(
    data_dir,
    "exp4_sanity_report.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT 4 SANITY REPORT\n"
)

cat(
  "============================================================\n\n"
)


cat(
  sprintf(
    "Theoretical per-agent score mean: %.6f\n\n",
    mean_score_theory
  )
)


cat(
  "Estimated one-agent score means:\n"
)


print(
  as.data.frame(
    mean_check
  ),
  digits = 6,
  row.names = FALSE
)


cat(
  "\nDirect Cramer-rate estimates:\n"
)


print(
  rate_table %>%
    select(
      sigma,
      sigma2,
      t_hat,
      I_hat,
      I_se
    ) %>%
    as.data.frame(),
  digits = 6,
  row.names = FALSE
)


cat(
  "\nFinite-n rate at the largest n:\n"
)


print(
  as.data.frame(
    finite_rate_check
  ),
  digits = 5,
  row.names = FALSE
)


cat(
  "\nMonotonicity check over the displayed sigma grid:\n"
)


print(
  as.data.frame(
    monotonicity_check
  ),
  row.names = FALSE
)


cat(
  paste0(
    "\nDescriptive affine fit:\n",
    "  I_hat = ",
    round(
      fit_intercept,
      6
    ),
    " + ",
    round(
      fit_slope,
      6
    ),
    " sigma^2\n",
    "  R^2 = ",
    round(
      fit_r_squared,
      6
    ),
    "\n"
  )
)


cat(
  paste0(
    "\nInterpretation:\n",
    "  * Estimated score means should be close to -0.25.\n",
    "  * Each t_hat should be strictly positive.\n",
    "  * The finite-n rate -log(P_n)/n should move toward the direct\n",
    "    Cramer-rate estimate as n grows, although rare-event Monte Carlo\n",
    "    becomes noisy when the acceptance count is small.\n",
    "  * A visually linear relationship between I_hat and sigma^2 is a\n",
    "    numerical feature of this grid, not a proved asymptotic law.\n"
  )
)


# ================================================================
# 15. COMPLETION REPORT
# ================================================================

cat(
  "\nFigures written to:\n"
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp4_implementation_vs_sigma.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp4_logaccept_vs_n.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp4_rate_vs_sigma2.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp4_finite_n_rate.pdf"
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