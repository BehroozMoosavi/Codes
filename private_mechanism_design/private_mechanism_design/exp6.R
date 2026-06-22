# ================================================================
# EXPERIMENT 6
# Hierarchical unknown-endpoint model and oracle recovery
#
# Optimal Binary Mechanism under Locally Private Signals
#
# PURPOSE
#
# This experiment studies whether a computational hierarchical rule,
# which learns an unknown upper support endpoint theta from privatized
# signals, approaches a known-theta oracle.
#
# The numerical model is
#
#     theta ~ Uniform[1,2],
#     theta_0 = 1.5,
#
# and, conditional on theta,
#
#     X = x_lower + (theta-x_lower) Z,
#
# where either
#
#     Z ~ Beta(1,1)
#
# or
#
#     Z ~ Beta(2,2).
#
# The lower support endpoint is
#
#     x_lower = 0.3 > 0.
#
# Consequently,
#
#     E[J_theta(X) | theta] = x_lower = 0.3,
#
# so positive reduced-form revenue is of order n rather than collapsing
# toward zero. This makes the relative oracle-hierarchical difference
#
#     (R_oracle - R_hierarchical) / R_oracle
#
# numerically well posed.
#
# SIGNAL CHANNEL
#
#     Y = X + N(0,sigma^2),
#     sigma = 0.45.
#
# SCORE
#
#     S_theta(y,lambda)
#       = E[X | theta,y]
#         + lambda E[J_theta(X) | theta,y],
#
# with lambda = 0.5.
#
# HIERARCHICAL RULE
#
# The planner computes the posterior p(theta | y_1,...,y_n) and averages
# the conditional posterior scores over theta. The project is built when
# the aggregate hierarchical score is nonnegative.
#
# ORACLE RULE
#
# The oracle knows theta_0 = 1.5 and uses the corresponding conditional
# posterior score directly.
#
# IMPORTANT SCOPE
#
# Because theta indexes the upper endpoint of the conditional support,
# this is a nonregular model and lies outside the paper's formal
# common-support hierarchical theorem.
#
# The log-log slope of the posterior standard deviation is therefore
# reported only as a descriptive numerical concentration diagnostic.
# The dashed n^{-1/2} line is a reference slope, not a claimed
# Bernstein-von Mises theorem.
#
# OUTPUT FOR EACH PRIOR
#
#   figures/fig_exp6_revenue_levels_<prior>.pdf
#   figures/fig_exp6_revenue_relative_gap_<prior>.pdf
#   figures/fig_exp6_theta_posterior_concentration_<prior>.pdf
#
#   tables/exp6_hierarchical_summary_<prior>.csv
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
  "R_codes/experiment/exp6"
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
    width = 7.4,
    height = 5.2
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


cat(
  "Output directory:",
  normalizePath(
    out_dir
  ),
  "\nFAST =",
  FAST,
  "\n\n"
)


# ================================================================
# 3. GENERAL NUMERICAL HELPERS
# ================================================================

log_sum_exp <- function(log_values) {
  maximum_value <- max(
    log_values
  )
  
  if (!is.finite(maximum_value)) {
    return(
      -Inf
    )
  }
  
  maximum_value +
    log(
      sum(
        exp(
          log_values -
            maximum_value
        )
      )
    )
}


safe_normalize_log_weights <- function(log_weights) {
  maximum_value <- max(
    log_weights
  )
  
  if (!is.finite(maximum_value)) {
    stop(
      "All hierarchical log weights are non-finite."
    )
  }
  
  weights <- exp(
    log_weights -
      maximum_value
  )
  
  weight_sum <- sum(
    weights
  )
  
  if (
    !is.finite(weight_sum) ||
    weight_sum <= .Machine$double.eps
  ) {
    stop(
      "Unable to normalize hierarchical posterior weights."
    )
  }
  
  weights /
    weight_sum
}


ratio_delta_standard_error <- function(
    numerator_draws,
    denominator_draws
) {
  if (
    length(numerator_draws) !=
    length(denominator_draws)
  ) {
    stop(
      "Numerator and denominator draws must have equal length."
    )
  }
  
  valid <- is.finite(
    numerator_draws
  ) &
    is.finite(
      denominator_draws
    )
  
  numerator_draws <- numerator_draws[
    valid
  ]
  
  denominator_draws <- denominator_draws[
    valid
  ]
  
  number_of_draws <- length(
    numerator_draws
  )
  
  if (number_of_draws < 2L) {
    return(
      NA_real_
    )
  }
  
  denominator_mean <- mean(
    denominator_draws
  )
  
  if (
    !is.finite(denominator_mean) ||
    abs(denominator_mean) <=
    .Machine$double.eps
  ) {
    return(
      NA_real_
    )
  }
  
  ratio_estimate <- mean(
    numerator_draws
  ) /
    denominator_mean
  
  influence_values <- (
    numerator_draws -
      ratio_estimate *
      denominator_draws
  ) /
    denominator_mean
  
  sd(
    influence_values
  ) /
    sqrt(
      number_of_draws
    )
}


# ================================================================
# 4. RUN EXPERIMENT FOR ONE CONDITIONAL PRIOR
# ================================================================

run_exp6 <- function(
    prior_label,
    beta_shape_1,
    beta_shape_2,
    file_suffix
) {
  set.seed(20260609)
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "EXPERIMENT 6:",
    prior_label,
    "\n"
  )
  
  cat(
    "============================================================\n\n"
  )
  
  
  # --------------------------------------------------------------
  # 4.1 MODEL PARAMETERS
  # --------------------------------------------------------------
  
  x_lower <- 0.3
  
  theta_min  <- 1.0
  theta_max  <- 2.0
  theta_true <- 1.5
  
  gaussian_sd <- 0.45
  
  lambda_value <- 0.5
  
  
  population_grid <- c(
    25,
    50,
    100,
    200,
    400,
    800
  )
  
  
  monte_carlo_replications <- if (FAST) {
    1000
  } else {
    4000
  }
  
  
  theta_grid <- seq(
    theta_min,
    theta_max,
    length.out = if (FAST) {
      41
    } else {
      81
    }
  )
  
  
  theta_prior_weights <- rep(
    1 /
      length(
        theta_grid
      ),
    length(
      theta_grid
    )
  )
  
  
  # Ten Gaussian standard deviations beyond the largest possible
  # conditional support reduce interpolation-tail error.
  
  signal_grid <- seq(
    x_lower -
      10 *
      gaussian_sd,
    theta_max +
      10 *
      gaussian_sd,
    length.out = if (FAST) {
      1401
    } else {
      2401
    }
  )
  
  
  # --------------------------------------------------------------
  # 4.2 CONDITIONAL PRIOR, VIRTUAL VALUE, AND SAMPLER
  # --------------------------------------------------------------
  
  conditional_density <- function(
    x,
    theta
  ) {
    support_length <- theta -
      x_lower
    
    standardized_x <- (
      x -
        x_lower
    ) /
      support_length
    
    density_values <- dbeta(
      standardized_x,
      beta_shape_1,
      beta_shape_2
    ) /
      support_length
    
    density_values[
      standardized_x < 0 |
        standardized_x > 1
    ] <- 0
    
    density_values
  }
  
  
  conditional_virtual_value <- function(
    x,
    theta
  ) {
    support_length <- theta -
      x_lower
    
    standardized_x <- (
      x -
        x_lower
    ) /
      support_length
    
    density_values <- dbeta(
      standardized_x,
      beta_shape_1,
      beta_shape_2
    ) /
      support_length
    
    distribution_values <- pbeta(
      standardized_x,
      beta_shape_1,
      beta_shape_2
    )
    
    virtual_values <- x -
      (
        1 -
          distribution_values
      ) /
      pmax(
        density_values,
        1e-12
      )
    
    virtual_values[
      standardized_x < 0 |
        standardized_x > 1
    ] <- NA_real_
    
    virtual_values
  }
  
  
  sample_conditional_types <- function(
    number_of_draws,
    theta
  ) {
    x_lower +
      (
        theta -
          x_lower
      ) *
      rbeta(
        number_of_draws,
        beta_shape_1,
        beta_shape_2
      )
  }
  
  
  # --------------------------------------------------------------
  # 4.3 CONDITIONAL POSTERIOR TABLE FOR A FIXED theta
  #
  # Columns:
  #
  #     marginal = m_theta(y)
  #     xhat     = E[X | theta,y]
  #     Jhat     = E[J_theta(X) | theta,y]
  #     score    = xhat + lambda Jhat
  #
  # --------------------------------------------------------------
  
  build_conditional_posterior_table <- function(
    theta
  ) {
    type_grid <- seq(
      x_lower,
      theta,
      length.out = if (FAST) {
        701
      } else {
        1201
      }
    )
    
    type_step <- type_grid[2] -
      type_grid[1]
    
    prior_values <- conditional_density(
      type_grid,
      theta
    )
    
    virtual_values <- conditional_virtual_value(
      type_grid,
      theta
    )
    
    
    # Avoid endpoint indeterminacies in inverse-hazard calculations.
    finite_virtual <- is.finite(
      virtual_values
    )
    
    if (!all(finite_virtual)) {
      virtual_values[
        !finite_virtual
      ] <- type_grid[
        !finite_virtual
      ]
    }
    
    
    rows <- lapply(
      signal_grid,
      function(signal_value) {
        log_kernel_values <- dnorm(
          signal_value,
          mean = type_grid,
          sd = gaussian_sd,
          log = TRUE
        )
        
        positive_prior <- prior_values > 0
        
        log_weights <- rep(
          -Inf,
          length(
            type_grid
          )
        )
        
        log_weights[
          positive_prior
        ] <-
          log_kernel_values[
            positive_prior
          ] +
          log(
            prior_values[
              positive_prior
            ]
          )
        
        
        maximum_log_weight <- max(
          log_weights
        )
        
        if (!is.finite(maximum_log_weight)) {
          return(
            c(
              marginal = NA_real_,
              xhat = NA_real_,
              Jhat = NA_real_,
              score = NA_real_
            )
          )
        }
        
        
        stabilized_weights <- exp(
          log_weights -
            maximum_log_weight
        )
        
        stabilized_integral <- sum(
          stabilized_weights
        ) *
          type_step
        
        
        if (
          !is.finite(stabilized_integral) ||
          stabilized_integral <=
          .Machine$double.eps
        ) {
          return(
            c(
              marginal = NA_real_,
              xhat = NA_real_,
              Jhat = NA_real_,
              score = NA_real_
            )
          )
        }
        
        
        log_marginal <-
          maximum_log_weight +
          log(
            stabilized_integral
          )
        
        
        posterior_x <- sum(
          type_grid *
            stabilized_weights
        ) *
          type_step /
          stabilized_integral
        
        
        posterior_J <- sum(
          virtual_values *
            stabilized_weights
        ) *
          type_step /
          stabilized_integral
        
        
        posterior_score <-
          posterior_x +
          lambda_value *
          posterior_J
        
        
        c(
          marginal = exp(
            log_marginal
          ),
          xhat = posterior_x,
          Jhat = posterior_J,
          score = posterior_score
        )
      }
    )
    
    
    output <- do.call(
      rbind,
      rows
    )
    
    output <- as.data.frame(
      output
    )
    
    output$y <- signal_grid
    
    output
  }
  
  
  # --------------------------------------------------------------
  # 4.4 PRECOMPUTE CONDITIONAL TABLES
  # --------------------------------------------------------------
  
  conditional_tables <- vector(
    mode = "list",
    length = length(
      theta_grid
    )
  )
  
  
  for (theta_index in seq_along(
    theta_grid
  )) {
    cat(
      "Building table",
      theta_index,
      "of",
      length(
        theta_grid
      ),
      "for theta =",
      round(
        theta_grid[theta_index],
        4
      ),
      "\n"
    )
    
    conditional_tables[[
      theta_index
    ]] <-
      build_conditional_posterior_table(
        theta_grid[
          theta_index
        ]
      )
  }
  
  
  interpolate_table_column <- function(
    theta_index,
    column_name,
    signal_values
  ) {
    selected_table <- conditional_tables[[
      theta_index
    ]]
    
    approx(
      x = selected_table$y,
      y = selected_table[[
        column_name
      ]],
      xout = signal_values,
      rule = 2,
      ties = "ordered"
    )$y
  }
  
  
  # --------------------------------------------------------------
  # 4.5 HIERARCHICAL POSTERIOR OVER theta
  # --------------------------------------------------------------
  
  hierarchical_theta_weights <- function(
    signal_values
  ) {
    log_weights <- vapply(
      seq_along(
        theta_grid
      ),
      function(theta_index) {
        marginal_values <- interpolate_table_column(
          theta_index = theta_index,
          column_name = "marginal",
          signal_values = signal_values
        )
        
        log(
          theta_prior_weights[
            theta_index
          ]
        ) +
          sum(
            log(
              pmax(
                marginal_values,
                1e-300
              )
            )
          )
      },
      FUN.VALUE = numeric(1)
    )
    
    safe_normalize_log_weights(
      log_weights
    )
  }
  
  
  evaluate_hierarchical_rule <- function(
    signal_values
  ) {
    number_of_signals <- length(
      signal_values
    )
    
    theta_weights <- hierarchical_theta_weights(
      signal_values
    )
    
    
    # Each matrix has:
    #
    #     rows    = signals i=1,...,n
    #     columns = theta-grid points.
    
    conditional_score_matrix <- vapply(
      seq_along(
        theta_grid
      ),
      function(theta_index) {
        interpolate_table_column(
          theta_index = theta_index,
          column_name = "score",
          signal_values = signal_values
        )
      },
      FUN.VALUE = numeric(
        number_of_signals
      )
    )
    
    
    conditional_J_matrix <- vapply(
      seq_along(
        theta_grid
      ),
      function(theta_index) {
        interpolate_table_column(
          theta_index = theta_index,
          column_name = "Jhat",
          signal_values = signal_values
        )
      },
      FUN.VALUE = numeric(
        number_of_signals
      )
    )
    
    
    # Handle n=1 explicitly because vapply may simplify dimensions.
    conditional_score_matrix <- matrix(
      conditional_score_matrix,
      nrow = number_of_signals,
      ncol = length(
        theta_grid
      )
    )
    
    conditional_J_matrix <- matrix(
      conditional_J_matrix,
      nrow = number_of_signals,
      ncol = length(
        theta_grid
      )
    )
    
    
    hierarchical_scores <- as.numeric(
      conditional_score_matrix %*%
        theta_weights
    )
    
    
    hierarchical_J_values <- as.numeric(
      conditional_J_matrix %*%
        theta_weights
    )
    
    
    theta_posterior_mean <- sum(
      theta_weights *
        theta_grid
    )
    
    
    theta_posterior_variance <- sum(
      theta_weights *
        (
          theta_grid -
            theta_posterior_mean
        )^2
    )
    
    
    list(
      build = as.numeric(
        sum(
          hierarchical_scores
        ) >= 0
      ),
      
      virtual_sum = sum(
        hierarchical_J_values
      ),
      
      theta_posterior_mean =
        theta_posterior_mean,
      
      theta_posterior_sd = sqrt(
        max(
          theta_posterior_variance,
          0
        )
      )
    )
  }
  
  
  # --------------------------------------------------------------
  # 4.6 KNOWN-theta ORACLE
  # --------------------------------------------------------------
  
  true_theta_index <- which.min(
    abs(
      theta_grid -
        theta_true
    )
  )
  
  
  theta_grid_true <- theta_grid[
    true_theta_index
  ]
  
  
  evaluate_oracle_rule <- function(
    signal_values
  ) {
    oracle_scores <- interpolate_table_column(
      theta_index = true_theta_index,
      column_name = "score",
      signal_values = signal_values
    )
    
    oracle_J_values <- interpolate_table_column(
      theta_index = true_theta_index,
      column_name = "Jhat",
      signal_values = signal_values
    )
    
    list(
      build = as.numeric(
        sum(
          oracle_scores
        ) >= 0
      ),
      
      virtual_sum = sum(
        oracle_J_values
      )
    )
  }
  
  
  # --------------------------------------------------------------
  # 4.7 REGIME CHECK
  # --------------------------------------------------------------
  
  beta_mean <- beta_shape_1 /
    (
      beta_shape_1 +
        beta_shape_2
    )
  
  
  theoretical_mean_type <-
    x_lower +
    (
      theta_true -
        x_lower
    ) *
    beta_mean
  
  
  theoretical_mean_virtual_value <- x_lower
  
  
  theoretical_mean_score <-
    theoretical_mean_type +
    lambda_value *
    theoretical_mean_virtual_value
  
  
  cat(
    sprintf(
      paste0(
        "\n[%s] E[X]=%.4f, E[J]=%.4f, ",
        "mean score=%.4f\n"
      ),
      prior_label,
      theoretical_mean_type,
      theoretical_mean_virtual_value,
      theoretical_mean_score
    )
  )
  
  
  cat(
    sprintf(
      paste0(
        "[%s] theta_true=%.4f; nearest theta-grid point=%.4f\n\n"
      ),
      prior_label,
      theta_true,
      theta_grid_true
    )
  )
  
  
  if (theoretical_mean_virtual_value <= 0) {
    warning(
      paste0(
        "The selected lower endpoint does not produce the intended ",
        "positive-revenue regime."
      )
    )
  }
  
  
  # --------------------------------------------------------------
  # 4.8 SIMULATE ONE POPULATION SIZE
  # --------------------------------------------------------------
  
  simulate_population_size <- function(
    population_size
  ) {
    oracle_revenue_draws <- numeric(
      monte_carlo_replications
    )
    
    hierarchical_revenue_draws <- numeric(
      monte_carlo_replications
    )
    
    theta_sd_draws <- numeric(
      monte_carlo_replications
    )
    
    theta_mean_draws <- numeric(
      monte_carlo_replications
    )
    
    allocation_agreement_draws <- logical(
      monte_carlo_replications
    )
    
    
    for (replication_index in seq_len(
      monte_carlo_replications
    )) {
      type_draws <- sample_conditional_types(
        number_of_draws = population_size,
        theta = theta_true
      )
      
      signal_draws <-
        type_draws +
        rnorm(
          population_size,
          mean = 0,
          sd = gaussian_sd
        )
      
      
      oracle_result <- evaluate_oracle_rule(
        signal_draws
      )
      
      
      hierarchical_result <- evaluate_hierarchical_rule(
        signal_draws
      )
      
      
      # Positive reduced-form virtual revenue associated with each
      # allocation rule.
      
      oracle_revenue_draws[
        replication_index
      ] <-
        pmax(
          oracle_result$virtual_sum,
          0
        ) *
        oracle_result$build
      
      
      hierarchical_revenue_draws[
        replication_index
      ] <-
        pmax(
          hierarchical_result$virtual_sum,
          0
        ) *
        hierarchical_result$build
      
      
      theta_sd_draws[
        replication_index
      ] <-
        hierarchical_result$theta_posterior_sd
      
      
      theta_mean_draws[
        replication_index
      ] <-
        hierarchical_result$theta_posterior_mean
      
      
      allocation_agreement_draws[
        replication_index
      ] <-
        oracle_result$build ==
        hierarchical_result$build
    }
    
    
    oracle_mean <- mean(
      oracle_revenue_draws
    )
    
    
    hierarchical_mean <- mean(
      hierarchical_revenue_draws
    )
    
    
    paired_difference_draws <-
      oracle_revenue_draws -
      hierarchical_revenue_draws
    
    
    relative_gap <-
      mean(
        paired_difference_draws
      ) /
      oracle_mean
    
    
    relative_gap_standard_error <-
      ratio_delta_standard_error(
        numerator_draws =
          paired_difference_draws,
        
        denominator_draws =
          oracle_revenue_draws
      )
    
    
    data.frame(
      n = population_size,
      
      Rplus_oracle =
        oracle_mean,
      
      Rplus_oracle_se =
        sd(
          oracle_revenue_draws
        ) /
        sqrt(
          monte_carlo_replications
        ),
      
      Rplus_hier =
        hierarchical_mean,
      
      Rplus_hier_se =
        sd(
          hierarchical_revenue_draws
        ) /
        sqrt(
          monte_carlo_replications
        ),
      
      relative_gap =
        relative_gap,
      
      relative_gap_se =
        relative_gap_standard_error,
      
      relative_gap_lower =
        relative_gap -
        1.96 *
        relative_gap_standard_error,
      
      relative_gap_upper =
        relative_gap +
        1.96 *
        relative_gap_standard_error,
      
      theta_posterior_sd =
        mean(
          theta_sd_draws
        ),
      
      theta_posterior_sd_se =
        sd(
          theta_sd_draws
        ) /
        sqrt(
          monte_carlo_replications
        ),
      
      theta_posterior_mean =
        mean(
          theta_mean_draws
        ),
      
      theta_posterior_bias =
        mean(
          theta_mean_draws
        ) -
        theta_true,
      
      allocation_agreement =
        mean(
          allocation_agreement_draws
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  
  # --------------------------------------------------------------
  # 4.9 RUN OVER n
  # --------------------------------------------------------------
  
  result_rows <- lapply(
    population_grid,
    function(population_size) {
      cat(
        "Simulating n =",
        population_size,
        "\n"
      )
      
      simulate_population_size(
        population_size
      )
    }
  )
  
  
  results <- bind_rows(
    result_rows
  )
  
  
  write.csv(
    results,
    file.path(
      data_dir,
      paste0(
        "exp6_hierarchical_summary_",
        file_suffix,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # --------------------------------------------------------------
  # 4.10 FIGURE A: REVENUE LEVELS
  # --------------------------------------------------------------
  
  revenue_plot_data <- results %>%
    select(
      n,
      Rplus_oracle,
      Rplus_oracle_se,
      Rplus_hier,
      Rplus_hier_se
    ) %>%
    pivot_longer(
      cols = c(
        Rplus_oracle,
        Rplus_hier
      ),
      names_to = "method",
      values_to = "positive_revenue"
    ) %>%
    mutate(
      positive_revenue_se = ifelse(
        method ==
          "Rplus_oracle",
        Rplus_oracle_se,
        Rplus_hier_se
      ),
      
      method = recode(
        method,
        Rplus_oracle =
          "known-theta oracle",
        Rplus_hier =
          "hierarchical rule"
      ),
      
      lower =
        pmax(
          0,
          positive_revenue -
            1.96 *
            positive_revenue_se
        ),
      
      upper =
        positive_revenue +
        1.96 *
        positive_revenue_se
    )
  
  
  revenue_plot <- ggplot(
    revenue_plot_data,
    aes(
      x = n,
      y = positive_revenue,
      color = method,
      group = method
    )
  ) +
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper,
        fill = method
      ),
      alpha = 0.10,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(
      linewidth = 1.05
    ) +
    geom_point(
      size = 2.6
    ) +
    scale_x_continuous(
      breaks = population_grid
    ) +
    scale_color_manual(
      values = c(
        "known-theta oracle" =
          "#0072B2",
        
        "hierarchical rule" =
          "#D55E00"
      ),
      name = NULL
    ) +
    scale_fill_manual(
      values = c(
        "known-theta oracle" =
          "#0072B2",
        
        "hierarchical rule" =
          "#D55E00"
      )
    ) +
    labs(
      x = "number of agents  n",
      y = expression(
        "positive reduced-form revenue  " *
          R["+"]
      ),
      subtitle = paste0(
        "Both quantities grow linearly in the positive-lower-endpoint ",
        "design"
      )
    ) +
    theme_paper()
  
  
  if (SHOW_TITLES) {
    revenue_plot <- revenue_plot +
      labs(
        title = paste0(
          "Oracle and hierarchical revenue: ",
          prior_label
        )
      )
  }
  
  
  save_figure(
    plot_object = revenue_plot,
    filename_base = paste0(
      "fig_exp6_revenue_levels_",
      file_suffix
    ),
    width = 7.4,
    height = 5.2
  )
  
  
  # --------------------------------------------------------------
  # 4.11 FIGURE B: RELATIVE REVENUE DIFFERENCE
  # --------------------------------------------------------------
  
  relative_gap_plot <- ggplot(
    results,
    aes(
      x = n,
      y = relative_gap
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray35",
      linewidth = 0.5
    ) +
    geom_ribbon(
      aes(
        ymin = relative_gap_lower,
        ymax = relative_gap_upper
      ),
      alpha = 0.15,
      fill = "#D55E00"
    ) +
    geom_line(
      linewidth = 1.05,
      color = "#D55E00"
    ) +
    geom_point(
      size = 2.6,
      color = "#D55E00"
    ) +
    scale_x_continuous(
      breaks = population_grid
    ) +
    labs(
      x = "number of agents  n",
      y = paste0(
        "relative revenue difference  ",
        "(oracle - hierarchical)/oracle"
      ),
      subtitle = paste0(
        "Paired delta-method confidence band; zero denotes equal ",
        "mean positive revenue"
      )
    ) +
    theme_paper()
  
  
  if (SHOW_TITLES) {
    relative_gap_plot <- relative_gap_plot +
      labs(
        title = paste0(
          "Relative revenue difference: ",
          prior_label
        )
      )
  }
  
  
  save_figure(
    plot_object = relative_gap_plot,
    filename_base = paste0(
      "fig_exp6_revenue_relative_gap_",
      file_suffix
    ),
    width = 7.4,
    height = 5.2
  )
  
  
  # --------------------------------------------------------------
  # 4.12 FIGURE C: POSTERIOR CONCENTRATION DIAGNOSTIC
  # --------------------------------------------------------------
  
  concentration_fit_data <- results %>%
    filter(
      n >= 50,
      theta_posterior_sd > 0
    )
  
  
  concentration_model <- lm(
    log(
      theta_posterior_sd
    ) ~
      log(
        n
      ),
    data = concentration_fit_data
  )
  
  
  fitted_slope <- unname(
    coef(
      concentration_model
    )[[
      "log(n)"
    ]]
  )
  
  
  first_reference_row <- concentration_fit_data[
    1,
    ,
    drop = FALSE
  ]
  
  
  reference_constant <-
    first_reference_row$theta_posterior_sd *
    sqrt(
      first_reference_row$n
    )
  
  
  concentration_plot_data <- results %>%
    mutate(
      reference_rate =
        reference_constant /
        sqrt(
          n
        ),
      
      posterior_sd_lower =
        pmax(
          theta_posterior_sd -
            1.96 *
            theta_posterior_sd_se,
          1e-8
        ),
      
      posterior_sd_upper =
        theta_posterior_sd +
        1.96 *
        theta_posterior_sd_se
    )
  
  
  concentration_plot <- ggplot(
    concentration_plot_data,
    aes(
      x = n,
      y = theta_posterior_sd
    )
  ) +
    geom_ribbon(
      aes(
        ymin = posterior_sd_lower,
        ymax = posterior_sd_upper
      ),
      alpha = 0.12,
      fill = "#1B9E77"
    ) +
    geom_line(
      aes(
        y = reference_rate
      ),
      linetype = "dashed",
      color = "gray35",
      linewidth = 0.7
    ) +
    geom_line(
      linewidth = 1.05,
      color = "#1B9E77"
    ) +
    geom_point(
      size = 2.6,
      color = "#1B9E77"
    ) +
    annotate(
      geom = "text",
      x = max(
        results$n
      ),
      y = max(
        results$theta_posterior_sd
      ),
      label = paste0(
        "fitted slope (n >= 50) = ",
        round(
          fitted_slope,
          2
        ),
        "\nreference slope = -1/2"
      ),
      hjust = 1,
      vjust = 1,
      family = "serif",
      size = 4.1
    ) +
    scale_x_log10(
      breaks = population_grid
    ) +
    scale_y_log10() +
    labs(
      x = "number of agents  n  (log scale)",
      y = expression(
        "posterior " *
          sd(
            theta
          ) *
          "  (log scale)"
      ),
      subtitle = paste0(
        "Descriptive concentration diagnostic for the nonregular ",
        "unknown-endpoint model"
      )
    ) +
    theme_paper()
  
  
  if (SHOW_TITLES) {
    concentration_plot <- concentration_plot +
      labs(
        title = paste0(
          "Posterior concentration: ",
          prior_label
        )
      )
  }
  
  
  save_figure(
    plot_object = concentration_plot,
    filename_base = paste0(
      "fig_exp6_theta_posterior_concentration_",
      file_suffix
    ),
    width = 7.4,
    height = 5.2
  )
  
  
  # --------------------------------------------------------------
  # 4.13 CONSOLE REPORT
  # --------------------------------------------------------------
  
  cat(
    "\nResults for",
    prior_label,
    ":\n\n"
  )
  
  
  print(
    results %>%
      select(
        n,
        Rplus_oracle,
        Rplus_hier,
        relative_gap,
        relative_gap_se,
        theta_posterior_mean,
        theta_posterior_bias,
        theta_posterior_sd,
        allocation_agreement
      ) %>%
      as.data.frame(),
    digits = 5,
    row.names = FALSE
  )
  
  
  cat(
    sprintf(
      paste0(
        "\nDescriptive fitted log-log slope of posterior sd(theta): ",
        "%.4f\n"
      ),
      fitted_slope
    )
  )
  
  
  cat(
    paste0(
      "\nInterpretation:\n",
      "  * The oracle and hierarchical positive-revenue quantities are\n",
      "    expected to remain close and grow approximately linearly in n.\n",
      "  * The relative difference may be positive or negative in finite\n",
      "    samples; zero means equal expected positive revenue.\n",
      "  * The confidence interval uses the paired ratio delta method and\n",
      "    accounts for randomness in the oracle denominator.\n",
      "  * The posterior concentration slope is descriptive only. The\n",
      "    unknown-endpoint model is nonregular, so no Bernstein-von Mises\n",
      "    conclusion is asserted.\n"
    )
  )
  
  
  invisible(
    results
  )
}


# ================================================================
# 5. RUN BOTH CONDITIONAL PRIOR SHAPES
# ================================================================

results_uniform <- run_exp6(
  prior_label = "Uniform",
  beta_shape_1 = 1,
  beta_shape_2 = 1,
  file_suffix = "uniform"
)


results_beta <- run_exp6(
  prior_label = "Beta(2,2)",
  beta_shape_1 = 2,
  beta_shape_2 = 2,
  file_suffix = "beta"
)


# ================================================================
# 6. COMPLETION REPORT
# ================================================================

cat(
  "\nExperiment 6 complete.\n"
)


cat(
  "Figures written to:\n  ",
  fig_dir,
  "\n",
  sep = ""
)


cat(
  "Tables written to:\n  ",
  data_dir,
  "\n",
  sep = ""
)
