# ================================================================
# EXPERIMENT 5
# Transfers as a Fredholm inverse problem
#
# Optimal Binary Mechanism under Locally Private Signals
#
# THEORY
#
# The opponent-averaged transfer tau solves the Fredholm equation
#
#     (K tau)(x)
#       = integral k(y | x) tau(y) dy
#       = T(x),
#
# where T(x) is the envelope-implied interim transfer schedule.
#
# Completeness of the continuous channel operator implies uniqueness of
# the opponent-averaged transfer whenever an implementing transfer exists.
# Completeness does not imply that every target belongs to the range.
#
# NUMERICAL INTERPRETATION
#
# This script studies a finite-dimensional discretization of the operator.
# The type grid contains more points than the signal grid, so the matrix K
# is tall:
#
#     number of type points > number of signal points.
#
# Therefore, full column rank of K is a meaningful finite-dimensional
# analogue of injectivity. It does not prove completeness of the continuous
# operator.
#
# The target T(x) is generated from a posterior-score allocation rule.
# Opponents' aggregate scores are approximated by a Gaussian CLT formula,
# so the target is itself a numerical approximation.
#
# The recovered transfer is obtained by second-difference Tikhonov
# regularization:
#
#     tau_rho
#       = argmin_tau
#           ||K tau - T||_2^2
#           + rho ||D2 tau||_2^2
#           + ridge ||tau||_2^2.
#
# OUTPUT
#
#   figures/fig_exp5_conditioning.pdf
#   figures/fig_exp5_stability_amplification.pdf
#   figures/fig_exp5_fredholm_fit.pdf
#   figures/fig_exp5_singular_values.pdf
#
#   tables/exp5_envelope_target_vs_implemented.csv
#   tables/exp5_singular_values.csv
#   tables/exp5_summary.csv
#   tables/exp5_rho_sweep.csv
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
  "R_codes/experiment/exp5"
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


# ================================================================
# 3. MODEL AND NUMERICAL DESIGN
# ================================================================

x_lower <- 0
x_upper <- 1

n_agents <- 50
lambda_value <- 1


sigma_grid <- c(
  0.10,
  0.20,
  0.40,
  0.70,
  1.00
)


rho_main <- 1e-4


rho_sweep <- c(
  1e-5,
  1e-4,
  1e-3
)


relative_perturbation_size <- 0.01


number_of_perturbations <- if (FAST) {
  60
} else {
  200
}


number_of_sweep_perturbations <- if (FAST) {
  30
} else {
  80
}


# Tall discretization:
#
#     n_x > n_y.
#
# This avoids the automatic nonuniqueness generated by a wide matrix.

number_of_type_points <- if (FAST) {
  220
} else {
  360
}


number_of_signal_points <- if (FAST) {
  160
} else {
  240
}


number_of_quadrature_points <- if (FAST) {
  2001
} else {
  5001
}


x_grid <- seq(
  x_lower,
  x_upper,
  length.out = number_of_type_points
)


dx <- x_grid[2] -
  x_grid[1]


quadrature_x <- seq(
  x_lower,
  x_upper,
  length.out = number_of_quadrature_points
)


quadrature_dx <- quadrature_x[2] -
  quadrature_x[1]


quadrature_virtual_value <- 2 *
  quadrature_x -
  x_upper


# ================================================================
# 4. BASIC HELPERS
# ================================================================

uniform_density <- function(x) {
  ifelse(
    x >= x_lower &
      x <= x_upper,
    1 /
      (
        x_upper -
          x_lower
      ),
    0
  )
}


trapezoid_cumulative <- function(
    function_values,
    grid
) {
  if (length(function_values) != length(grid)) {
    stop(
      "function_values and grid must have the same length."
    )
  }
  
  if (length(grid) < 2L) {
    stop(
      "The grid must contain at least two points."
    )
  }
  
  increments <- diff(
    grid
  ) *
    (
      function_values[-1] +
        function_values[-length(function_values)]
    ) /
    2
  
  c(
    0,
    cumsum(
      increments
    )
  )
}


euclidean_norm <- function(vector) {
  sqrt(
    sum(
      vector^2
    )
  )
}


relative_error <- function(
    estimate,
    target
) {
  denominator <- euclidean_norm(
    target
  )
  
  if (denominator <= .Machine$double.eps) {
    return(
      NA_real_
    )
  }
  
  euclidean_norm(
    estimate -
      target
  ) /
    denominator
}


# ================================================================
# 5. SIGMA-SPECIFIC SIGNAL GRID
# ================================================================

make_signal_grid <- function(
    sigma_value
) {
  # The signal range extends ten standard deviations beyond the type
  # support. This makes Gaussian tail truncation negligible for the
  # numerical quadrature used here.
  
  seq(
    x_lower -
      10 *
      sigma_value,
    x_upper +
      10 *
      sigma_value,
    length.out = number_of_signal_points
  )
}


# ================================================================
# 6. POSTERIOR MOMENTS ON THE SIGNAL GRID
# ================================================================

posterior_moments_on_grid <- function(
    sigma_value,
    signal_grid
) {
  prior_values <- uniform_density(
    quadrature_x
  )
  
  rows <- lapply(
    signal_grid,
    function(signal_value) {
      log_weights <- dnorm(
        signal_value,
        mean = quadrature_x,
        sd = sigma_value,
        log = TRUE
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
      
      stabilized_weights <- exp(
        log_weights -
          maximum_log_weight
      )
      
      denominator <- sum(
        stabilized_weights
      ) *
        quadrature_dx
      
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
        quadrature_x *
          stabilized_weights
      ) *
        quadrature_dx /
        denominator
      
      posterior_J <- sum(
        quadrature_virtual_value *
          stabilized_weights
      ) *
        quadrature_dx /
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


# ================================================================
# 7. DISCRETIZED CHANNEL OPERATOR
# ================================================================

make_channel_matrix <- function(
    sigma_value,
    signal_grid
) {
  signal_step <- signal_grid[2] -
    signal_grid[1]
  
  outer(
    x_grid,
    signal_grid,
    function(type_value, signal_value) {
      dnorm(
        signal_value,
        mean = type_value,
        sd = sigma_value
      )
    }
  ) *
    signal_step
}


# ================================================================
# 8. SECOND-DIFFERENCE PENALTY
# ================================================================

make_second_difference_matrix <- function(
    dimension
) {
  if (dimension < 3L) {
    stop(
      "The signal grid must have at least three points."
    )
  }
  
  difference_matrix <- matrix(
    0,
    nrow = dimension - 2L,
    ncol = dimension
  )
  
  for (index in seq_len(
    dimension -
    2L
  )) {
    difference_matrix[
      index,
      index
    ] <- 1
    
    difference_matrix[
      index,
      index + 1L
    ] <- -2
    
    difference_matrix[
      index,
      index + 2L
    ] <- 1
  }
  
  difference_matrix
}


# ================================================================
# 9. REGULARIZED SOLVER
# ================================================================

regularized_solve_operator <- function(
    channel_matrix,
    second_difference_matrix,
    rho,
    ridge = 1e-10
) {
  normal_matrix <-
    crossprod(
      channel_matrix
    ) +
    rho *
    crossprod(
      second_difference_matrix
    ) +
    ridge *
    diag(
      ncol(
        channel_matrix
      )
    )
  
  solve(
    normal_matrix,
    t(
      channel_matrix
    )
  )
}


regularized_solution <- function(
    channel_matrix,
    second_difference_matrix,
    target,
    rho,
    ridge = 1e-10
) {
  solve_operator <- regularized_solve_operator(
    channel_matrix = channel_matrix,
    second_difference_matrix = second_difference_matrix,
    rho = rho,
    ridge = ridge
  )
  
  transfer <- as.numeric(
    solve_operator %*%
      target
  )
  
  implemented_target <- as.numeric(
    channel_matrix %*%
      transfer
  )
  
  list(
    operator = solve_operator,
    transfer = transfer,
    implemented_target = implemented_target
  )
}


# ================================================================
# 10. CLT-BASED ENVELOPE TARGET
#
# The interim allocation Q(x) is approximated as follows.
#
# Let
#
#     S(y) = xhat(y) + lambda Jhat(y).
#
# For an agent with signal y, the project is implemented when
#
#     S(y) + sum_{j != i} S(Y_j) >= 0.
#
# The opponents' score sum is approximated by a Gaussian distribution,
#
#     sum_{j != i} S(Y_j)
#       approximately
#       Normal(
#         (n-1) mu_S,
#         (n-1) sigma_S^2
#       ).
#
# This gives a signal-contingent pivotal probability p(y). Integrating
# p(y) against k(y|x) gives the approximate interim allocation Q(x).
# ================================================================

construct_target <- function(
    sigma_value,
    signal_grid,
    posterior_table,
    channel_matrix
) {
  signal_step <- signal_grid[2] -
    signal_grid[1]
  
  posterior_score <- posterior_table$score
  
  if (any(!is.finite(posterior_score))) {
    stop(
      paste0(
        "Non-finite posterior scores at sigma = ",
        sigma_value,
        "."
      )
    )
  }
  
  # Marginal density m(y).
  marginal_density <- vapply(
    signal_grid,
    function(signal_value) {
      sum(
        dnorm(
          signal_value,
          mean = quadrature_x,
          sd = sigma_value
        ) *
          uniform_density(
            quadrature_x
          )
      ) *
        quadrature_dx
    },
    FUN.VALUE = numeric(1)
  )
  
  marginal_mass <- sum(
    marginal_density
  ) *
    signal_step
  
  if (
    !is.finite(marginal_mass) ||
    marginal_mass <= .Machine$double.eps
  ) {
    stop(
      "The discretized marginal signal density has invalid total mass."
    )
  }
  
  marginal_density <- marginal_density /
    marginal_mass
  
  mean_score <- sum(
    posterior_score *
      marginal_density
  ) *
    signal_step
  
  score_variance <- sum(
    (
      posterior_score -
        mean_score
    )^2 *
      marginal_density
  ) *
    signal_step
  
  score_sd <- sqrt(
    max(
      score_variance,
      1e-12
    )
  )
  
  # Centered implementation threshold.
  #
  # When E[S] != 0 the opponent-sum mean (n-1)*mean_score dominates the
  # O(1) own-score, driving the pivotal probability (and hence Q(x)) to a
  # corner; the envelope transfer T(x) = x Q(x) - int Q then collapses to
  # the constant x_lo and the reconstruction panel carries no signal.
  # Centering the threshold at the opponent-sum mean (as in exp6) removes
  # this mechanical degeneracy, so the pivot is Phi(S(y) / (sqrt(n-1) sigma_S))
  # and Q(x) varies over (0,1), yielding a nondegenerate target to invert.
  pivotal_probability <- pnorm(
    posterior_score /
      (
        sqrt(
          n_agents -
            1
        ) *
          score_sd
      )
  )
  
  interim_allocation <- as.numeric(
    channel_matrix %*%
      pivotal_probability
  )
  
  cumulative_allocation <- trapezoid_cumulative(
    function_values = interim_allocation,
    grid = x_grid
  )
  
  target_transfer <-
    x_grid *
    interim_allocation -
    cumulative_allocation
  
  list(
    target = target_transfer,
    interim_allocation = interim_allocation,
    pivotal_probability = pivotal_probability,
    marginal_density = marginal_density,
    mean_score = mean_score,
    score_sd = score_sd
  )
}


# ================================================================
# 11. PERTURBATION-AMPLIFICATION DIAGNOSTIC
# ================================================================

perturbation_amplification <- function(
    solve_operator,
    target,
    transfer,
    number_of_draws,
    relative_size
) {
  target_norm <- euclidean_norm(
    target
  )
  
  transfer_norm <- euclidean_norm(
    transfer
  )
  
  if (
    target_norm <= .Machine$double.eps ||
    transfer_norm <= .Machine$double.eps
  ) {
    return(
      rep(
        NA_real_,
        number_of_draws
      )
    )
  }
  
  amplification_values <- numeric(
    number_of_draws
  )
  
  for (draw_index in seq_len(
    number_of_draws
  )) {
    random_direction <- rnorm(
      length(
        target
      )
    )
    
    random_direction_norm <- euclidean_norm(
      random_direction
    )
    
    random_direction <- random_direction /
      random_direction_norm
    
    target_perturbation <-
      relative_size *
      target_norm *
      random_direction
    
    transfer_perturbation <- as.numeric(
      solve_operator %*%
        target_perturbation
    )
    
    relative_target_error <-
      euclidean_norm(
        target_perturbation
      ) /
      target_norm
    
    relative_transfer_error <-
      euclidean_norm(
        transfer_perturbation
      ) /
      transfer_norm
    
    amplification_values[
      draw_index
    ] <-
      relative_transfer_error /
      relative_target_error
  }
  
  amplification_values
}


# ================================================================
# 12. MAIN LOOP
# ================================================================

implementation_rows <- list()
singular_value_rows <- list()
summary_rows <- list()

main_objects <- list()


for (sigma_value in sigma_grid) {
  cat(
    "\nProcessing sigma =",
    sigma_value,
    "\n"
  )
  
  signal_grid <- make_signal_grid(
    sigma_value
  )
  
  posterior_table <- posterior_moments_on_grid(
    sigma_value = sigma_value,
    signal_grid = signal_grid
  )
  
  channel_matrix <- make_channel_matrix(
    sigma_value = sigma_value,
    signal_grid = signal_grid
  )
  
  second_difference_matrix <- make_second_difference_matrix(
    ncol(
      channel_matrix
    )
  )
  
  target_object <- construct_target(
    sigma_value = sigma_value,
    signal_grid = signal_grid,
    posterior_table = posterior_table,
    channel_matrix = channel_matrix
  )
  
  target_transfer <- target_object$target
  
  solution_object <- regularized_solution(
    channel_matrix = channel_matrix,
    second_difference_matrix = second_difference_matrix,
    target = target_transfer,
    rho = rho_main
  )
  
  recovered_transfer <- solution_object$transfer
  
  implemented_target <- solution_object$implemented_target
  
  singular_values <- svd(
    channel_matrix,
    nu = 0,
    nv = 0
  )$d
  
  numerical_rank_tolerance <-
    max(
      dim(
        channel_matrix
      )
    ) *
    max(
      singular_values
    ) *
    .Machine$double.eps
  
  numerical_rank <- sum(
    singular_values >
      numerical_rank_tolerance
  )
  
  full_column_rank <-
    numerical_rank ==
    ncol(
      channel_matrix
    )
  
  condition_number <- max(
    singular_values
  ) /
    min(
      singular_values
    )
  
  amplification_draws <- perturbation_amplification(
    solve_operator = solution_object$operator,
    target = target_transfer,
    transfer = recovered_transfer,
    number_of_draws = number_of_perturbations,
    relative_size = relative_perturbation_size
  )
  
  implementation_rows[[
    as.character(
      sigma_value
    )
  ]] <- data.frame(
    x = x_grid,
    T_target = target_transfer,
    T_implemented = implemented_target,
    sigma = sigma_value,
    stringsAsFactors = FALSE
  )
  
  singular_value_rows[[
    as.character(
      sigma_value
    )
  ]] <- data.frame(
    sigma = sigma_value,
    index = seq_along(
      singular_values
    ),
    singular_value = singular_values,
    normalized_singular_value =
      singular_values /
      max(
        singular_values
      ),
    stringsAsFactors = FALSE
  )
  
  summary_rows[[
    as.character(
      sigma_value
    )
  ]] <- data.frame(
    sigma = sigma_value,
    
    number_of_type_points =
      nrow(
        channel_matrix
      ),
    
    number_of_signal_points =
      ncol(
        channel_matrix
      ),
    
    numerical_rank =
      numerical_rank,
    
    full_column_rank =
      full_column_rank,
    
    smallest_singular_value =
      min(
        singular_values
      ),
    
    largest_singular_value =
      max(
        singular_values
      ),
    
    condition_number =
      condition_number,
    
    target_relative_fit_error =
      relative_error(
        implemented_target,
        target_transfer
      ),
    
    transfer_l2_norm =
      euclidean_norm(
        recovered_transfer
      ),
    
    target_l2_norm =
      euclidean_norm(
        target_transfer
      ),
    
    roughness_l2_norm =
      euclidean_norm(
        second_difference_matrix %*%
          recovered_transfer
      ),
    
    amplification_mean =
      mean(
        amplification_draws,
        na.rm = TRUE
      ),
    
    amplification_sd =
      sd(
        amplification_draws,
        na.rm = TRUE
      ),
    
    mean_score =
      target_object$mean_score,
    
    score_sd =
      target_object$score_sd,
    
    rho =
      rho_main,
    
    stringsAsFactors = FALSE
  )
  
  main_objects[[
    as.character(
      sigma_value
    )
  ]] <- list(
    signal_grid = signal_grid,
    posterior_table = posterior_table,
    channel_matrix = channel_matrix,
    second_difference_matrix = second_difference_matrix,
    target_object = target_object,
    solution_object = solution_object
  )
}


implementation_table <- bind_rows(
  implementation_rows
)


singular_value_table <- bind_rows(
  singular_value_rows
)


summary_table <- bind_rows(
  summary_rows
)


write.csv(
  implementation_table,
  file.path(
    data_dir,
    "exp5_envelope_target_vs_implemented.csv"
  ),
  row.names = FALSE
)


write.csv(
  singular_value_table,
  file.path(
    data_dir,
    "exp5_singular_values.csv"
  ),
  row.names = FALSE
)


write.csv(
  summary_table,
  file.path(
    data_dir,
    "exp5_summary.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 13. REGULARIZATION SWEEP
# ================================================================

sweep_rows <- list()
sweep_counter <- 0L


for (sigma_value in sigma_grid) {
  base_object <- main_objects[[
    as.character(
      sigma_value
    )
  ]]
  
  channel_matrix <- base_object$channel_matrix
  
  second_difference_matrix <-
    base_object$second_difference_matrix
  
  target_transfer <-
    base_object$target_object$target
  
  for (rho_value in rho_sweep) {
    sweep_counter <- sweep_counter +
      1L
    
    solution_object <- regularized_solution(
      channel_matrix = channel_matrix,
      second_difference_matrix = second_difference_matrix,
      target = target_transfer,
      rho = rho_value
    )
    
    amplification_draws <- perturbation_amplification(
      solve_operator = solution_object$operator,
      target = target_transfer,
      transfer = solution_object$transfer,
      number_of_draws = number_of_sweep_perturbations,
      relative_size = relative_perturbation_size
    )
    
    sweep_rows[[
      sweep_counter
    ]] <- data.frame(
      sigma = sigma_value,
      rho = rho_value,
      
      target_relative_fit_error =
        relative_error(
          solution_object$implemented_target,
          target_transfer
        ),
      
      transfer_l2_norm =
        euclidean_norm(
          solution_object$transfer
        ),
      
      roughness_l2_norm =
        euclidean_norm(
          second_difference_matrix %*%
            solution_object$transfer
        ),
      
      amplification_mean =
        mean(
          amplification_draws,
          na.rm = TRUE
        ),
      
      amplification_sd =
        sd(
          amplification_draws,
          na.rm = TRUE
        ),
      
      stringsAsFactors = FALSE
    )
  }
}


sweep_table <- bind_rows(
  sweep_rows
)


write.csv(
  sweep_table,
  file.path(
    data_dir,
    "exp5_rho_sweep.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 14. COLORS
# ================================================================

sigma_colors <- c(
  "0.1" = "#0072B2",
  "0.2" = "#1B9E77",
  "0.4" = "#E6AB02",
  "0.7" = "#D55E00",
  "1"   = "#7B3294"
)


# ================================================================
# 15. FIGURE 1
# CONDITIONING OF THE DISCRETIZED OPERATOR
# ================================================================

conditioning_long <- summary_table %>%
  select(
    sigma,
    smallest_singular_value,
    condition_number
  ) %>%
  pivot_longer(
    cols = c(
      smallest_singular_value,
      condition_number
    ),
    names_to = "quantity",
    values_to = "value"
  ) %>%
  mutate(
    quantity = recode(
      quantity,
      
      smallest_singular_value =
        "smallest singular value",
      
      condition_number =
        "condition number"
    )
  )


conditioning_plot <- ggplot(
  conditioning_long,
  aes(
    x = sigma,
    y = value,
    color = quantity,
    group = quantity
  )
) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    size = 2.6
  ) +
  scale_y_log10(
    labels = label_log()
  ) +
  scale_x_continuous(
    breaks = sigma_grid
  ) +
  scale_color_manual(
    values = c(
      "smallest singular value" = "#0072B2",
      "condition number" = "#D55E00"
    ),
    name = NULL
  ) +
  labs(
    x = expression(
      "Gaussian noise level  " *
        sigma
    ),
    y = "value (log scale)",
    subtitle = paste0(
      "Finite-grid conditioning only; full column rank is ",
      "consistent with, but does not prove, continuous injectivity"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  conditioning_plot <- conditioning_plot +
    labs(
      title = "Conditioning of the discretized Fredholm operator"
    )
}


save_figure(
  plot_object = conditioning_plot,
  filename_base = "fig_exp5_conditioning",
  width = 7.4,
  height = 5.2
)


# ================================================================
# 16. FIGURE 2
# PERTURBATION AMPLIFICATION
# ================================================================

amplification_table <- summary_table %>%
  mutate(
    lower =
      pmax(
        amplification_mean -
          amplification_sd,
        1e-8
      ),
    
    upper =
      amplification_mean +
      amplification_sd
  )


amplification_plot <- ggplot(
  amplification_table,
  aes(
    x = sigma,
    y = amplification_mean
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15,
    fill = "#D55E00"
  ) +
  geom_line(
    linewidth = 1.05,
    color = "#D55E00"
  ) +
  geom_point(
    size = 2.8,
    color = "#D55E00"
  ) +
  scale_y_log10(
    labels = label_log()
  ) +
  scale_x_continuous(
    breaks = sigma_grid
  ) +
  labs(
    x = expression(
      "Gaussian noise level  " *
        sigma
    ),
    y = "relative-error amplification (log scale)",
    subtitle = paste0(
      "Mean and one-standard-deviation band over random 1% ",
      "target perturbations"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  amplification_plot <- amplification_plot +
    labs(
      title = "Numerical instability of transfer recovery"
    )
}


save_figure(
  plot_object = amplification_plot,
  filename_base = "fig_exp5_stability_amplification",
  width = 7.4,
  height = 5.2
)


# ================================================================
# 17. FIGURE 3
# TARGET AND REGULARIZED IMPLEMENTATION
# ================================================================

fit_plot <- ggplot(
  implementation_table,
  aes(
    x = x
  )
) +
  geom_line(
    aes(
      y = T_target
    ),
    linewidth = 1.25,
    color = "black"
  ) +
  geom_line(
    aes(
      y = T_implemented,
      color = factor(
        sigma
      ),
      group = factor(
        sigma
      )
    ),
    linewidth = 0.90,
    linetype = "22"
  ) +
  scale_color_manual(
    values = sigma_colors,
    name = expression(
      sigma
    )
  ) +
  labs(
    x = "type  x",
    y = expression(
      "interim transfer  " *
        T(x)
    ),
    subtitle = paste0(
      "Black: CLT-based envelope target; dashed: regularized ",
      "discretized implementation"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  fit_plot <- fit_plot +
    labs(
      title = "Regularized Fredholm transfer recovery"
    )
}


save_figure(
  plot_object = fit_plot,
  filename_base = "fig_exp5_fredholm_fit",
  width = 7.6,
  height = 5.4
)


# ================================================================
# 18. FIGURE 4
# NORMALIZED SINGULAR-VALUE SPECTRA
# ================================================================

maximum_displayed_index <- min(
  160,
  max(
    singular_value_table$index
  )
)


singular_value_plot <- ggplot(
  singular_value_table %>%
    filter(
      index <=
        maximum_displayed_index
    ),
  aes(
    x = index,
    y = normalized_singular_value,
    color = factor(
      sigma
    ),
    group = factor(
      sigma
    )
  )
) +
  geom_line(
    linewidth = 1
  ) +
  scale_y_log10(
    labels = label_scientific()
  ) +
  scale_color_manual(
    values = sigma_colors,
    name = expression(
      sigma
    )
  ) +
  labs(
    x = "singular-value index",
    y = "normalized singular value (log scale)",
    subtitle = paste0(
      "Faster decay indicates increasing numerical ill-conditioning"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  singular_value_plot <- singular_value_plot +
    labs(
      title = "Singular spectrum of the discretized operator"
    )
}


save_figure(
  plot_object = singular_value_plot,
  filename_base = "fig_exp5_singular_values",
  width = 7.6,
  height = 5.4
)


# ================================================================
# 19. SANITY REPORT
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT 5 SANITY REPORT\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Main-grid summary at rho =",
  rho_main,
  ":\n\n"
)


print(
  summary_table %>%
    select(
      sigma,
      number_of_type_points,
      number_of_signal_points,
      numerical_rank,
      full_column_rank,
      smallest_singular_value,
      condition_number,
      target_relative_fit_error,
      transfer_l2_norm,
      roughness_l2_norm,
      amplification_mean,
      amplification_sd,
      mean_score,
      score_sd
    ) %>%
    as.data.frame(),
  digits = 5,
  row.names = FALSE
)


cat(
  paste0(
    "\nInterpretation:\n",
    "  * full_column_rank = TRUE means the selected finite matrix has\n",
    "    no numerically detected column nullspace.\n",
    "  * This is a property of the selected discretization and does not\n",
    "    establish completeness of the continuous operator.\n",
    "  * A declining smallest singular value and growing condition number\n",
    "    indicate worsening numerical ill-conditioning.\n",
    "  * target_relative_fit_error measures approximate range recovery for\n",
    "    the selected regularization parameter and CLT-based target.\n",
    "  * amplification_mean measures sensitivity of the regularized\n",
    "    recovered transfer to small perturbations of the target.\n"
  )
)


cat(
  "\nRegularization sweep:\n\n"
)


sweep_display <- sweep_table %>%
  select(
    sigma,
    rho,
    target_relative_fit_error,
    transfer_l2_norm,
    roughness_l2_norm,
    amplification_mean
  ) %>%
  arrange(
    sigma,
    rho
  )


print(
  as.data.frame(
    sweep_display
  ),
  digits = 5,
  row.names = FALSE
)


cat(
  paste0(
    "\nRegularization tradeoff:\n",
    "  * Smaller rho generally improves in-sample fit but may produce a\n",
    "    larger, rougher, and less stable recovered transfer.\n",
    "  * Larger rho generally stabilizes and smooths the inverse at the\n",
    "    cost of a larger target-fit error.\n"
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
    "fig_exp5_conditioning.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp5_stability_amplification.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp5_fredholm_fit.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp5_singular_values.pdf"
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