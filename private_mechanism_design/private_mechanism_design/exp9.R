# ================================================================
# EXPERIMENT 9
#
# Laplace versus logistic under three privacy calibrations,
# with Gaussian benchmarks on the approximate-DP and GDP axes
#
# Optimal Binary Mechanism under Locally Private Signals
#
# SCIENTIFIC OBJECTS
#
# (A) BINARY-ENDPOINT TRADE-OFF CURVES
#
#     At common tight mu-GDP, the calibrated endpoint curves satisfy
#
#         T_Log(alpha) <= T_Lap(alpha)
#
#     for all alpha in [0,1].
#
#     Thus, on the endpoint state space {-1,1}, the common-mu
#     logistic experiment Blackwell-dominates the common-mu
#     Laplace experiment.
#
# (B) ENDPOINT WELFARE
#
#     For every binary endpoint decision problem,
#
#         W_Log^* >= W_Lap^*.
#
# (C) CONTINUOUS-TYPE WELFARE
#
#     X_i ~ Uniform[-1,1], lambda = 0, and
#
#         q(Y_1,...,Y_n)
#           =
#           1{
#             sum_i E[X_i | Y_i] >= 0
#           }.
#
#     The endpoint theorem does not establish a Blackwell ordering
#     of the complete continuous-type experiments. Continuous-type
#     comparisons are numerical and specification-specific.
#
# PRIVACY AXES
#
# Axis 1:
#     Equal pure-epsilon LDP scale:
#
#         b = beta = Delta/epsilon.
#
#     Only Laplace and logistic appear.
#
# Axis 2:
#     Common (epsilon,delta)-LDP boundary.
#
# Axis 3:
#     Common tight mu-GDP boundary.
#
# NUMERICAL FEATURES
#
#   * common type draws and common underlying uniforms;
#   * ratio delta-method standard errors;
#   * paired Laplace-minus-logistic standard errors;
#   * log-stable trapezoidal posterior quadrature;
#   * channel-specific adaptive signal grids;
#   * overflow diagnostics;
#   * direct trade-off and privacy-profile calibration checks;
#   * ordinary linear plot axes.
#
# ================================================================


rm(list = ls())


# ================================================================
# 0. USER OPTIONS
# ================================================================

FAST        <- FALSE
SHOW_TITLES <- FALSE

set.seed(20260612)


out_dir <- paste0(
  "/Users/behroozmoosavi/Desktop/privacy/codes/",
  "R_codes/experiment/exp9"
)

fig_dir <- file.path(
  out_dir,
  "figures"
)

data_dir <- file.path(
  out_dir,
  "tables"
)


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
  if (!requireNamespace(
    package_name,
    quietly = TRUE
  )) {
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
    width = 9,
    height = 5
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
# 3. TYPE SPACE AND CHANNEL CONSTANTS
# ================================================================

x_lower <- -1
x_upper <-  1

support_width <- x_upper - x_lower


clamp_fraction <- 0.20

clamp_lower <- x_lower +
  clamp_fraction *
  support_width

clamp_upper <- x_upper -
  clamp_fraction *
  support_width

clamped_sensitivity <- clamp_upper -
  clamp_lower


clamp_function <- function(
    x,
    lower,
    upper
) {
  pmin(
    pmax(
      x,
      lower
    ),
    upper
  )
}


channels <- c(
  "Gaussian",
  "Laplace",
  "Logistic"
)


channel_colors <- c(
  "Gaussian" = "#0072B2",
  "Laplace" = "#D55E00",
  "Logistic" = "#1B9E77"
)


privacy_axes <- c(
  "axis1",
  "axis2",
  "axis3"
)


axis_labels <- c(
  "axis1" = "Equal pure-epsilon LDP",
  "axis2" = "Common (epsilon, delta)-LDP",
  "axis3" = "Common tight mu-GDP"
)


# ================================================================
# 4. GENERAL NUMERICAL UTILITIES
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
    probability,
    scale
) {
  probability <- safe_probability(
    probability
  )
  
  ifelse(
    probability < 0.5,
    scale * log(2 * probability),
    -scale * log(2 * (1 - probability))
  )
}


ratio_estimate_and_se <- function(
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
  
  valid <- is.finite(numerator_draws) &
    is.finite(denominator_draws)
  
  numerator_draws <- numerator_draws[valid]
  denominator_draws <- denominator_draws[valid]
  
  number_of_draws <- length(
    numerator_draws
  )
  
  if (number_of_draws < 2L) {
    return(
      c(
        estimate = NA_real_,
        standard_error = NA_real_
      )
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
      c(
        estimate = NA_real_,
        standard_error = NA_real_
      )
    )
  }
  
  estimate <- mean(
    numerator_draws
  ) /
    denominator_mean
  
  influence_values <- (
    numerator_draws -
      estimate *
      denominator_draws
  ) /
    denominator_mean
  
  c(
    estimate = estimate,
    standard_error =
      sd(influence_values) /
      sqrt(number_of_draws)
  )
}


# ================================================================
# 5. STANDARDIZED LAPLACE AND LOGISTIC DENSITIES
# ================================================================

standard_laplace_density <- function(u) {
  0.5 *
    exp(
      -abs(u)
    )
}


standard_logistic_density <- function(u) {
  dlogis(
    u,
    location = 0,
    scale = 1
  )
}


standard_log_density <- function(
    u,
    family
) {
  switch(
    family,
    
    "lap" =
      -log(2) -
      abs(u),
    
    "log" =
      dlogis(
        u,
        location = 0,
        scale = 1,
        log = TRUE
      ),
    
    stop(
      "Unknown standardized family."
    )
  )
}


# ================================================================
# 6. HOCKEY-STICK PRIVACY PROFILES
# ================================================================

standard_integration_grid <- seq(
  -80,
  80,
  length.out = if (FAST) {
    8001
  } else {
    16001
  }
)


standard_grid_step <-
  standard_integration_grid[2] -
  standard_integration_grid[1]


standard_trapezoid_weights <- rep(
  1,
  length(standard_integration_grid)
)


standard_trapezoid_weights[
  c(
    1,
    length(standard_trapezoid_weights)
  )
] <- 0.5


density_mass_check <- function(family) {
  density_values <- switch(
    family,
    
    "lap" =
      standard_laplace_density(
        standard_integration_grid
      ),
    
    "log" =
      standard_logistic_density(
        standard_integration_grid
      ),
    
    stop(
      "Unknown family in density_mass_check()."
    )
  )
  
  sum(
    density_values *
      standard_trapezoid_weights
  ) *
    standard_grid_step
}


delta_standardized <- function(
    epsilon,
    separation,
    family
) {
  log_p <- standard_log_density(
    standard_integration_grid -
      separation,
    family
  )
  
  log_q <- standard_log_density(
    standard_integration_grid,
    family
  )
  
  maximum_log <- pmax(
    log_p,
    epsilon + log_q
  )
  
  positive_difference <- exp(
    maximum_log
  ) *
    pmax(
      exp(log_p - maximum_log) -
        exp(
          epsilon +
            log_q -
            maximum_log
        ),
      0
    )
  
  sum(
    positive_difference *
      standard_trapezoid_weights
  ) *
    standard_grid_step
}


# ================================================================
# 7. GAUSSIAN-DP PRIVACY PROFILE
# ================================================================

delta_gdp <- function(
    epsilon,
    mu
) {
  if (
    !is.finite(mu) ||
    mu <= 0
  ) {
    stop(
      "mu must be positive and finite."
    )
  }
  
  value <- pnorm(
    -epsilon / mu +
      mu / 2
  ) -
    exp(epsilon) *
    pnorm(
      -epsilon / mu -
        mu / 2
    )
  
  pmax(
    value,
    0
  )
}


mu_from_epsilon_delta <- function(
    epsilon,
    delta
) {
  objective <- function(mu) {
    delta_gdp(
      epsilon = epsilon,
      mu = mu
    ) -
      delta
  }
  
  lower <- 1e-8
  upper <- 1
  
  while (
    objective(upper) < 0
  ) {
    upper <- 2 * upper
    
    if (upper > 500) {
      stop(
        "Could not bracket the GDP index."
      )
    }
  }
  
  uniroot(
    f = objective,
    interval = c(
      lower,
      upper
    ),
    tol = 1e-11
  )$root
}


# ================================================================
# 8. CLOSED-FORM TIGHT GDP CALIBRATION
# ================================================================

r_laplace_star <- function(mu) {
  if (
    !is.finite(mu) ||
    mu <= 0
  ) {
    stop(
      "mu must be positive and finite."
    )
  }
  
  -2 *
    log(
      2 *
        pnorm(
          -mu / 2
        )
    )
}


r_logistic_star <- function(mu) {
  if (
    !is.finite(mu) ||
    mu <= 0
  ) {
    stop(
      "mu must be positive and finite."
    )
  }
  
  2 *
    (
      log(
        pnorm(
          mu / 2
        )
      ) -
        log(
          pnorm(
            -mu / 2
          )
        )
    )
}


# ================================================================
# 9. BINARY TRADE-OFF FUNCTIONS
# ================================================================

tradeoff_laplace <- function(
    separation,
    alpha
) {
  alpha <- pmin(
    pmax(alpha, 0),
    1
  )
  
  cutoff <- 0.5 *
    exp(
      -separation
    )
  
  output <- numeric(
    length(alpha)
  )
  
  region_1 <- alpha <= cutoff
  region_2 <- alpha > cutoff &
    alpha <= 0.5
  region_3 <- alpha > 0.5
  
  output[region_1] <-
    1 -
    exp(separation) *
    alpha[region_1]
  
  output[region_2] <-
    exp(-separation) /
    (
      4 *
        alpha[region_2]
    )
  
  output[region_3] <-
    exp(-separation) *
    (
      1 -
        alpha[region_3]
    )
  
  pmin(
    pmax(output, 0),
    1
  )
}


tradeoff_logistic <- function(
    separation,
    alpha
) {
  alpha <- pmin(
    pmax(alpha, 0),
    1
  )
  
  denominator <-
    1 -
    alpha +
    exp(separation) *
    alpha
  
  output <- (
    1 -
      alpha
  ) /
    denominator
  
  pmin(
    pmax(output, 0),
    1
  )
}


tradeoff_gdp <- function(
    mu,
    alpha
) {
  alpha_safe <- safe_probability(
    alpha
  )
  
  pnorm(
    qnorm(
      1 -
        alpha_safe
    ) -
      mu
  )
}


# ================================================================
# 10. DIRECT GDP INDEX FROM A TRADE-OFF CURVE
# ================================================================

direct_mu_from_tradeoff <- function(
    separation,
    family
) {
  tradeoff_function <- switch(
    family,
    
    "lap" =
      tradeoff_laplace,
    
    "log" =
      tradeoff_logistic,
    
    stop(
      "Unknown family in direct_mu_from_tradeoff()."
    )
  )
  
  objective <- function(alpha) {
    tradeoff_value <- tradeoff_function(
      separation = separation,
      alpha = alpha
    )
    
    qnorm(
      1 -
        alpha
    ) -
      qnorm(
        safe_probability(
          tradeoff_value
        )
      )
  }
  
  coarse_alpha_grid <- seq(
    1e-8,
    1 - 1e-8,
    length.out = if (FAST) {
      1501
    } else {
      5001
    }
  )
  
  coarse_values <- vapply(
    coarse_alpha_grid,
    objective,
    numeric(1)
  )
  
  maximum_index <- which.max(
    coarse_values
  )
  
  lower_index <- max(
    1,
    maximum_index - 2
  )
  
  upper_index <- min(
    length(coarse_alpha_grid),
    maximum_index + 2
  )
  
  optimized <- optimize(
    f = objective,
    interval = c(
      coarse_alpha_grid[lower_index],
      coarse_alpha_grid[upper_index]
    ),
    maximum = TRUE,
    tol = 1e-12
  )
  
  list(
    mu = optimized$objective,
    alpha_star = optimized$maximum
  )
}


# ================================================================
# 11. PROFILE-BASED GDP CHECK
# ================================================================

profile_mu_check <- function(
    separation,
    family
) {
  epsilon_check_grid <- seq(
    0,
    max(
      10,
      separation + 6
    ),
    length.out = if (FAST) {
      120
    } else {
      300
    }
  )
  
  channel_delta <- vapply(
    epsilon_check_grid,
    function(epsilon_value) {
      delta_standardized(
        epsilon = epsilon_value,
        separation = separation,
        family = family
      )
    },
    numeric(1)
  )
  
  maximum_violation <- function(mu) {
    gaussian_delta <- vapply(
      epsilon_check_grid,
      function(epsilon_value) {
        delta_gdp(
          epsilon = epsilon_value,
          mu = mu
        )
      },
      numeric(1)
    )
    
    max(
      channel_delta -
        gaussian_delta
    )
  }
  
  lower <- 1e-5
  upper <- 1
  
  while (
    maximum_violation(upper) > 0
  ) {
    upper <- 2 * upper
    
    if (upper > 200) {
      stop(
        "Could not bracket profile-based GDP index."
      )
    }
  }
  
  if (
    maximum_violation(lower) <= 0
  ) {
    return(lower)
  }
  
  uniroot(
    f = maximum_violation,
    interval = c(
      lower,
      upper
    ),
    tol = 1e-8
  )$root
}


# ================================================================
# 12. APPROXIMATE-DP CALIBRATION
# ================================================================

standardized_separation_for_delta <- function(
    epsilon,
    delta,
    family
) {
  objective <- function(separation) {
    delta_standardized(
      epsilon = epsilon,
      separation = separation,
      family = family
    ) -
      delta
  }
  
  lower <- max(
    epsilon *
      (1 - 1e-8),
    1e-10
  )
  
  upper <- max(
    epsilon + 1,
    1
  )
  
  while (
    objective(upper) < 0
  ) {
    upper <- 2 * upper
    
    if (upper > 500) {
      stop(
        paste(
          "Could not bracket approximate-DP separation for",
          family
        )
      )
    }
  }
  
  uniroot(
    f = objective,
    interval = c(
      lower,
      upper
    ),
    tol = 1e-9
  )$root
}


# ================================================================
# 13. CHANNEL CALIBRATION
# ================================================================

calibrate_channel <- function(
    axis_name,
    channel_name,
    epsilon_value,
    delta_value
) {
  if (
    axis_name == "axis1"
  ) {
    return(
      switch(
        channel_name,
        
        "Laplace" =
          support_width /
          epsilon_value,
        
        "Logistic" =
          support_width /
          epsilon_value,
        
        NA_real_
      )
    )
  }
  
  if (
    axis_name == "axis2"
  ) {
    equivalent_mu <- mu_from_epsilon_delta(
      epsilon = epsilon_value,
      delta = delta_value
    )
    
    return(
      switch(
        channel_name,
        
        "Gaussian" =
          support_width /
          equivalent_mu,
        
        "Clamped Gaussian" =
          clamped_sensitivity /
          equivalent_mu,
        
        "Laplace" =
          support_width /
          standardized_separation_for_delta(
            epsilon = epsilon_value,
            delta = delta_value,
            family = "lap"
          ),
        
        "Logistic" =
          support_width /
          standardized_separation_for_delta(
            epsilon = epsilon_value,
            delta = delta_value,
            family = "log"
          ),
        
        stop(
          "Unknown channel."
        )
      )
    )
  }
  
  if (
    axis_name == "axis3"
  ) {
    common_mu <- mu_from_epsilon_delta(
      epsilon = epsilon_value,
      delta = delta_value
    )
    
    return(
      switch(
        channel_name,
        
        "Gaussian" =
          support_width /
          common_mu,
        
        "Clamped Gaussian" =
          clamped_sensitivity /
          common_mu,
        
        "Laplace" =
          support_width /
          r_laplace_star(
            common_mu
          ),
        
        "Logistic" =
          support_width /
          r_logistic_star(
            common_mu
          ),
        
        stop(
          "Unknown channel."
        )
      )
    )
  }
  
  stop(
    "Unknown privacy axis."
  )
}


# ================================================================
# 14. PRIVACY AND NOISE SUMMARIES
# ================================================================

noise_standard_deviation <- function(
    channel_name,
    scale_value
) {
  switch(
    channel_name,
    
    "Gaussian" =
      scale_value,
    
    "Clamped Gaussian" =
      scale_value,
    
    "Laplace" =
      sqrt(2) *
      scale_value,
    
    "Logistic" =
      pi *
      scale_value /
      sqrt(3),
    
    stop(
      "Unknown channel."
    )
  )
}


pure_epsilon_index <- function(
    channel_name,
    scale_value
) {
  if (
    channel_name %in%
    c(
      "Laplace",
      "Logistic"
    )
  ) {
    return(
      support_width /
        scale_value
    )
  }
  
  Inf
}


tight_mu_index <- function(
    channel_name,
    scale_value
) {
  if (
    !is.finite(scale_value) ||
    scale_value <= 0
  ) {
    return(NA_real_)
  }
  
  switch(
    channel_name,
    
    "Gaussian" =
      support_width /
      scale_value,
    
    "Clamped Gaussian" =
      clamped_sensitivity /
      scale_value,
    
    "Laplace" = {
      separation <- support_width /
        scale_value
      
      -2 *
        qnorm(
          safe_probability(
            0.5 *
              exp(
                -separation / 2
              )
          )
        )
    },
    
    "Logistic" = {
      separation <- support_width /
        scale_value
      
      -2 *
        qnorm(
          safe_probability(
            plogis(
              -separation / 2
            )
          )
        )
    },
    
    stop(
      "Unknown channel."
    )
  )
}


# ================================================================
# 15. CHANNEL LOG-DENSITIES
# ================================================================

channel_log_density <- function(
    signal_value,
    type_grid,
    channel_name,
    scale_value
) {
  transformed_type <- if (
    channel_name == "Clamped Gaussian"
  ) {
    clamp_function(
      type_grid,
      clamp_lower,
      clamp_upper
    )
  } else {
    type_grid
  }
  
  switch(
    channel_name,
    
    "Gaussian" =
      dnorm(
        signal_value,
        mean = transformed_type,
        sd = scale_value,
        log = TRUE
      ),
    
    "Clamped Gaussian" =
      dnorm(
        signal_value,
        mean = transformed_type,
        sd = scale_value,
        log = TRUE
      ),
    
    "Laplace" =
      -log(
        2 *
          scale_value
      ) -
      abs(
        signal_value -
          transformed_type
      ) /
      scale_value,
    
    "Logistic" =
      dlogis(
        signal_value -
          transformed_type,
        location = 0,
        scale = scale_value,
        log = TRUE
      ),
    
    stop(
      "Unknown channel."
    )
  )
}


# ================================================================
# 16. COMMON-RANDOM-NUMBER CHANNEL SIMULATION
# ================================================================

noise_from_uniform <- function(
    uniform_matrix,
    channel_name,
    scale_value
) {
  uniform_matrix <- safe_probability(
    uniform_matrix
  )
  
  switch(
    channel_name,
    
    "Gaussian" =
      qnorm(
        uniform_matrix,
        mean = 0,
        sd = scale_value
      ),
    
    "Clamped Gaussian" =
      qnorm(
        uniform_matrix,
        mean = 0,
        sd = scale_value
      ),
    
    "Laplace" =
      laplace_quantile(
        probability = uniform_matrix,
        scale = scale_value
      ),
    
    "Logistic" =
      qlogis(
        uniform_matrix,
        location = 0,
        scale = scale_value
      ),
    
    stop(
      "Unknown channel."
    )
  )
}


simulate_channel_crn <- function(
    type_matrix,
    uniform_matrix,
    channel_name,
    scale_value
) {
  channel_noise <- noise_from_uniform(
    uniform_matrix = uniform_matrix,
    channel_name = channel_name,
    scale_value = scale_value
  )
  
  transformed_types <- if (
    channel_name == "Clamped Gaussian"
  ) {
    clamp_function(
      type_matrix,
      clamp_lower,
      clamp_upper
    )
  } else {
    type_matrix
  }
  
  transformed_types +
    channel_noise
}


# ================================================================
# 17. CHANNEL-SPECIFIC SIGNAL GRID
# ================================================================

make_signal_grid <- function(
    channel_name,
    scale_value,
    grid_size,
    tail_probability = 1e-9
) {
  tail_probability <- max(
    tail_probability,
    1e-12
  )
  
  tail_width <- switch(
    channel_name,
    
    "Gaussian" =
      qnorm(
        1 -
          tail_probability,
        mean = 0,
        sd = scale_value
      ),
    
    "Clamped Gaussian" =
      qnorm(
        1 -
          tail_probability,
        mean = 0,
        sd = scale_value
      ),
    
    "Laplace" =
      -scale_value *
      log(
        2 *
          tail_probability
      ),
    
    "Logistic" =
      qlogis(
        1 -
          tail_probability,
        location = 0,
        scale = scale_value
      ),
    
    stop(
      "Unknown channel."
    )
  )
  
  signal_center_lower <- if (
    channel_name == "Clamped Gaussian"
  ) {
    clamp_lower
  } else {
    x_lower
  }
  
  signal_center_upper <- if (
    channel_name == "Clamped Gaussian"
  ) {
    clamp_upper
  } else {
    x_upper
  }
  
  seq(
    signal_center_lower -
      tail_width,
    signal_center_upper +
      tail_width,
    length.out = grid_size
  )
}


# ================================================================
# 18. LOG-STABLE TRAPEZOIDAL POSTERIOR MEAN TABLE
# ================================================================

posterior_mean_table <- function(
    channel_name,
    scale_value,
    signal_grid,
    type_grid
) {
  type_step <- type_grid[2] -
    type_grid[1]
  
  trapezoid_weights <- rep(
    1,
    length(type_grid)
  )
  
  trapezoid_weights[
    c(
      1,
      length(trapezoid_weights)
    )
  ] <- 0.5
  
  log_prior_density <- rep(
    -log(support_width),
    length(type_grid)
  )
  
  output <- vapply(
    signal_grid,
    function(signal_value) {
      log_weights <- channel_log_density(
        signal_value = signal_value,
        type_grid = type_grid,
        channel_name = channel_name,
        scale_value = scale_value
      ) +
        log_prior_density +
        log(
          trapezoid_weights
        )
      
      maximum_log_weight <- max(
        log_weights
      )
      
      if (!is.finite(maximum_log_weight)) {
        return(NA_real_)
      }
      
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
        denominator <=
        .Machine$double.eps
      ) {
        return(NA_real_)
      }
      
      sum(
        type_grid *
          stabilized_weights
      ) *
        type_step /
        denominator
    },
    numeric(1)
  )
  
  if (
    any(
      !is.finite(output)
    )
  ) {
    stop(
      paste0(
        "Posterior table contains non-finite values for ",
        channel_name,
        " at scale ",
        format(
          scale_value,
          digits = 8
        ),
        "."
      )
    )
  }
  
  output
}


# ================================================================
# 19. CONTINUOUS-TYPE WELFARE DRAWS
# ================================================================

continuous_welfare_draws <- function(
    channel_name,
    scale_value,
    type_matrix,
    uniform_matrix,
    aggregate_true_type,
    type_grid,
    signal_grid_size
) {
  signal_grid <- make_signal_grid(
    channel_name = channel_name,
    scale_value = scale_value,
    grid_size = signal_grid_size
  )
  
  posterior_mean <- posterior_mean_table(
    channel_name = channel_name,
    scale_value = scale_value,
    signal_grid = signal_grid,
    type_grid = type_grid
  )
  
  signal_matrix <- simulate_channel_crn(
    type_matrix = type_matrix,
    uniform_matrix = uniform_matrix,
    channel_name = channel_name,
    scale_value = scale_value
  )
  
  signal_vector <- as.vector(
    signal_matrix
  )
  
  overflow_rate <- mean(
    signal_vector <
      min(signal_grid) |
      signal_vector >
      max(signal_grid)
  )
  
  interpolated_posterior <- approx(
    x = signal_grid,
    y = posterior_mean,
    xout = signal_vector,
    rule = 2,
    ties = "ordered"
  )$y
  
  if (
    any(
      !is.finite(
        interpolated_posterior
      )
    )
  ) {
    stop(
      paste0(
        "Non-finite interpolated posterior means for ",
        channel_name,
        "."
      )
    )
  }
  
  posterior_matrix <- matrix(
    interpolated_posterior,
    nrow = nrow(type_matrix),
    ncol = ncol(type_matrix)
  )
  
  aggregate_posterior <- rowSums(
    posterior_matrix
  )
  
  build_indicator <- as.numeric(
    aggregate_posterior >= 0
  )
  
  list(
    welfare_draws =
      aggregate_true_type *
      build_indicator,
    
    build_probability =
      mean(build_indicator),
    
    overflow_rate =
      overflow_rate,
    
    mean_aggregate_posterior =
      mean(aggregate_posterior)
  )
}


# ================================================================
# 20. ENDPOINT WELFARE
# ================================================================

endpoint_efficiency <- function(
    family,
    separation,
    low_type,
    high_type,
    high_type_probability
) {
  if (
    !(family %in%
      c(
        "lap",
        "log"
      ))
  ) {
    stop(
      "family must be 'lap' or 'log'."
    )
  }
  
  if (
    high_type_probability <= 0 ||
    high_type_probability >= 1
  ) {
    stop(
      "The high-type probability must lie in (0,1)."
    )
  }
  
  tradeoff_function <- if (
    family == "lap"
  ) {
    tradeoff_laplace
  } else {
    tradeoff_logistic
  }
  
  loss_function <- function(alpha) {
    high_type_probability *
      high_type *
      tradeoff_function(
        separation = separation,
        alpha = alpha
      ) +
      (
        1 -
          high_type_probability
      ) *
      abs(low_type) *
      alpha
  }
  
  interior_solution <- optimize(
    f = loss_function,
    interval = c(
      1e-12,
      1 - 1e-12
    ),
    tol = 1e-12
  )
  
  candidate_losses <- c(
    alpha_zero =
      loss_function(0),
    
    interior =
      interior_solution$objective,
    
    alpha_one =
      loss_function(1)
  )
  
  optimal_candidate <- names(
    which.min(candidate_losses)
  )
  
  optimal_alpha <- switch(
    optimal_candidate,
    
    "alpha_zero" =
      0,
    
    "interior" =
      interior_solution$minimum,
    
    "alpha_one" =
      1
  )
  
  minimum_loss <- min(
    candidate_losses
  )
  
  first_best <-
    high_type_probability *
    high_type
  
  list(
    efficiency =
      (
        first_best -
          minimum_loss
      ) /
      first_best,
    
    alpha_star =
      optimal_alpha,
    
    loss =
      minimum_loss,
    
    candidate =
      optimal_candidate
  )
}


# ================================================================
# 21. MAIN SETTINGS
# ================================================================

delta_fixed <- 1e-5


epsilon_grid <- if (FAST) {
  c(
    0.5,
    1,
    2
  )
} else {
  c(
    0.5,
    0.75,
    1,
    1.5,
    2,
    3
  )
}


number_of_agents <- if (FAST) {
  40
} else {
  100
}


monte_carlo_replications <- if (FAST) {
  2500
} else {
  10000
}


type_grid <- seq(
  x_lower,
  x_upper,
  length.out = if (FAST) {
    801
  } else {
    1601
  }
)


signal_grid_size <- if (FAST) {
  1801
} else {
  3601
}


# ================================================================
# 22. COMMON MONTE CARLO DRAWS
# ================================================================

type_matrix <- matrix(
  runif(
    monte_carlo_replications *
      number_of_agents,
    min = x_lower,
    max = x_upper
  ),
  nrow = monte_carlo_replications,
  ncol = number_of_agents
)


uniform_matrix <- matrix(
  runif(
    monte_carlo_replications *
      number_of_agents
  ),
  nrow = monte_carlo_replications,
  ncol = number_of_agents
)


aggregate_true_type <- rowSums(
  type_matrix
)


first_best_draws <- pmax(
  aggregate_true_type,
  0
)


first_best_welfare <- mean(
  first_best_draws
)


cat(
  sprintf(
    paste0(
      "[exp9] Uniform[%.0f,%.0f], ",
      "Delta = %.1f, Delta_c = %.2f, ",
      "n = %d, B = %d, delta = %.0e, ",
      "sample W_FB = %.6f\n"
    ),
    x_lower,
    x_upper,
    support_width,
    clamped_sensitivity,
    number_of_agents,
    monte_carlo_replications,
    delta_fixed,
    first_best_welfare
  )
)


# ================================================================
# 23. NUMERICAL CALIBRATION CHECKS
# ================================================================

cat(
  "\n[check 0] Standardized density masses:\n"
)


cat(
  sprintf(
    "  Laplace:  %.12f\n",
    density_mass_check("lap")
  )
)


cat(
  sprintf(
    "  Logistic: %.12f\n",
    density_mass_check("log")
  )
)


gdp_check_rows <- list()
gdp_check_counter <- 0L


for (mu_value in c(
  0.25,
  0.5,
  1,
  2
)) {
  laplace_separation <- r_laplace_star(
    mu_value
  )
  
  logistic_separation <- r_logistic_star(
    mu_value
  )
  
  laplace_direct <- direct_mu_from_tradeoff(
    separation = laplace_separation,
    family = "lap"
  )
  
  logistic_direct <- direct_mu_from_tradeoff(
    separation = logistic_separation,
    family = "log"
  )
  
  laplace_profile <- profile_mu_check(
    separation = laplace_separation,
    family = "lap"
  )
  
  logistic_profile <- profile_mu_check(
    separation = logistic_separation,
    family = "log"
  )
  
  gdp_check_counter <- gdp_check_counter + 1L
  
  gdp_check_rows[[
    gdp_check_counter
  ]] <- data.frame(
    mu_target =
      mu_value,
    
    r_laplace =
      laplace_separation,
    
    r_logistic =
      logistic_separation,
    
    mu_laplace_direct =
      laplace_direct$mu,
    
    alpha_laplace_direct =
      laplace_direct$alpha_star,
    
    mu_logistic_direct =
      logistic_direct$mu,
    
    alpha_logistic_direct =
      logistic_direct$alpha_star,
    
    mu_laplace_profile =
      laplace_profile,
    
    mu_logistic_profile =
      logistic_profile,
    
    stringsAsFactors = FALSE
  )
}


gdp_check_results <- bind_rows(
  gdp_check_rows
)


write.csv(
  gdp_check_results,
  file.path(
    data_dir,
    "exp9_gdp_calibration_checks.csv"
  ),
  row.names = FALSE
)


cat(
  "\n[check 1] Tight-GDP calibration checks:\n"
)


print(
  gdp_check_results,
  digits = 8,
  row.names = FALSE
)


# ================================================================
# 24. TRADE-OFF CURVES AT COMMON TIGHT mu
# ================================================================

tradeoff_mu_grid <- c(
  0.5,
  1,
  2
)


tradeoff_alpha_grid <- seq(
  0,
  1,
  length.out = if (FAST) {
    501
  } else {
    2001
  }
)


tradeoff_rows <- list()
tradeoff_counter <- 0L


for (mu_value in tradeoff_mu_grid) {
  laplace_separation <- r_laplace_star(
    mu_value
  )
  
  logistic_separation <- r_logistic_star(
    mu_value
  )
  
  laplace_curve <- tradeoff_laplace(
    separation = laplace_separation,
    alpha = tradeoff_alpha_grid
  )
  
  logistic_curve <- tradeoff_logistic(
    separation = logistic_separation,
    alpha = tradeoff_alpha_grid
  )
  
  gaussian_curve <- tradeoff_gdp(
    mu = mu_value,
    alpha = tradeoff_alpha_grid
  )
  
  tradeoff_counter <- tradeoff_counter + 1L
  
  tradeoff_rows[[
    tradeoff_counter
  ]] <- data.frame(
    mu = mu_value,
    alpha = tradeoff_alpha_grid,
    Laplace = laplace_curve,
    Logistic = logistic_curve,
    GDP = gaussian_curve,
    tradeoff_gap =
      laplace_curve -
      logistic_curve,
    stringsAsFactors = FALSE
  )
}


tradeoff_results <- bind_rows(
  tradeoff_rows
) %>%
  mutate(
    mu_label = factor(
      paste0(
        "mu = ",
        mu
      ),
      levels = paste0(
        "mu = ",
        tradeoff_mu_grid
      )
    )
  )


write.csv(
  tradeoff_results,
  file.path(
    data_dir,
    "exp9_tradeoff_common_mu.csv"
  ),
  row.names = FALSE
)


tradeoff_dominance_summary <- tradeoff_results %>%
  group_by(
    mu
  ) %>%
  summarize(
    minimum_gap =
      min(tradeoff_gap),
    
    maximum_gap =
      max(tradeoff_gap),
    
    endpoint_dominance_holds =
      minimum_gap >= -1e-8,
    
    maximum_numerical_violation =
      max(
        -minimum_gap,
        0
      ),
    
    fraction_strictly_positive =
      mean(
        tradeoff_gap >
          1e-8
      ),
    
    fraction_numerically_zero =
      mean(
        abs(tradeoff_gap) <=
          1e-8
      ),
    
    .groups = "drop"
  )


write.csv(
  tradeoff_dominance_summary,
  file.path(
    data_dir,
    "exp9_tradeoff_dominance_summary.csv"
  ),
  row.names = FALSE
)


cat(
  "\nEndpoint trade-off dominance summary:\n"
)


print(
  tradeoff_dominance_summary,
  digits = 10,
  row.names = FALSE
)


if (
  any(
    !tradeoff_dominance_summary$
    endpoint_dominance_holds
  )
) {
  warning(
    paste0(
      "The numerical trade-off grid reports a negative value of ",
      "T_Lap(alpha)-T_Log(alpha) beyond tolerance. ",
      "Increase the alpha-grid resolution before interpreting it ",
      "as a theoretical violation."
    )
  )
}


tradeoff_long <- tradeoff_results %>%
  select(
    mu_label,
    alpha,
    Laplace,
    Logistic,
    GDP
  ) %>%
  pivot_longer(
    cols = c(
      Laplace,
      Logistic,
      GDP
    ),
    names_to = "experiment",
    values_to = "type_II_error"
  ) %>%
  mutate(
    experiment = factor(
      experiment,
      levels = c(
        "Laplace",
        "Logistic",
        "GDP"
      )
    )
  )


tradeoff_colors <- c(
  "Laplace" = "#D55E00",
  "Logistic" = "#1B9E77",
  "GDP" = "#0072B2"
)


tradeoff_plot <- ggplot(
  tradeoff_long,
  aes(
    x = alpha,
    y = type_II_error,
    color = experiment,
    group = experiment
  )
) +
  geom_line(
    linewidth = 1.05
  ) +
  facet_wrap(
    ~ mu_label,
    nrow = 1
  ) +
  scale_color_manual(
    values = tradeoff_colors,
    name = "Experiment"
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  labs(
    x = expression(
      "Type-I error  " *
        alpha
    ),
    y = expression(
      "Optimal Type-II error  " *
        T(alpha)
    ),
    subtitle = paste0(
      "Common tight mu-GDP calibration; ",
      "the logistic endpoint curve lies weakly below Laplace"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  tradeoff_plot <- tradeoff_plot +
    labs(
      title = "Binary endpoint trade-off curves"
    )
}


save_figure(
  plot_object = tradeoff_plot,
  filename_base = "fig_exp9_tradeoff_common_mu",
  width = 10,
  height = 4.7
)


tradeoff_gap_plot <- ggplot(
  tradeoff_results,
  aes(
    x = alpha,
    y = tradeoff_gap,
    color = mu_label,
    group = mu_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray35",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  scale_color_manual(
    values = c(
      "mu = 0.5" = "#1B9E77",
      "mu = 1" = "#D55E00",
      "mu = 2" = "#7570B3"
    ),
    name = NULL
  ) +
  labs(
    x = expression(
      "Type-I error  " *
        alpha
    ),
    y = expression(
      T[Lap](alpha) -
        T[Log](alpha)
    ),
    subtitle = paste0(
      "A nonnegative gap confirms logistic endpoint ROC dominance"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  tradeoff_gap_plot <- tradeoff_gap_plot +
    labs(
      title = "Laplace-minus-logistic endpoint trade-off gap"
    )
}


save_figure(
  plot_object = tradeoff_gap_plot,
  filename_base = "fig_exp9_tradeoff_gap_common_mu",
  width = 8.8,
  height = 5
)


# ================================================================
# 25. ENDPOINT WELFARE
# ================================================================

prior_high_grid <- sort(
  unique(
    c(
      seq(
        0.05,
        0.95,
        by = 0.02
      ),
      0.5
    )
  )
)


endpoint_mu_grid <- c(
  0.5,
  1,
  2
)


endpoint_rows <- list()
endpoint_counter <- 0L


for (mu_value in endpoint_mu_grid) {
  laplace_separation <- r_laplace_star(
    mu_value
  )
  
  logistic_separation <- r_logistic_star(
    mu_value
  )
  
  for (prior_high in prior_high_grid) {
    laplace_result <- endpoint_efficiency(
      family = "lap",
      separation = laplace_separation,
      low_type = -1,
      high_type = 1,
      high_type_probability = prior_high
    )
    
    logistic_result <- endpoint_efficiency(
      family = "log",
      separation = logistic_separation,
      low_type = -1,
      high_type = 1,
      high_type_probability = prior_high
    )
    
    endpoint_counter <- endpoint_counter + 1L
    
    endpoint_rows[[
      endpoint_counter
    ]] <- data.frame(
      mu =
        mu_value,
      
      prior_high =
        prior_high,
      
      efficiency_laplace =
        laplace_result$efficiency,
      
      efficiency_logistic =
        logistic_result$efficiency,
      
      welfare_gap =
        logistic_result$efficiency -
        laplace_result$efficiency,
      
      alpha_laplace =
        laplace_result$alpha_star,
      
      alpha_logistic =
        logistic_result$alpha_star,
      
      candidate_laplace =
        laplace_result$candidate,
      
      candidate_logistic =
        logistic_result$candidate,
      
      stringsAsFactors = FALSE
    )
  }
}


endpoint_results <- bind_rows(
  endpoint_rows
) %>%
  mutate(
    mu_label = factor(
      paste0(
        "mu = ",
        mu
      ),
      levels = paste0(
        "mu = ",
        endpoint_mu_grid
      )
    )
  )


write.csv(
  endpoint_results,
  file.path(
    data_dir,
    "exp9_endpoint_welfare.csv"
  ),
  row.names = FALSE
)


endpoint_gap_summary <- endpoint_results %>%
  group_by(
    mu
  ) %>%
  summarize(
    minimum_gap =
      min(welfare_gap),
    
    maximum_gap =
      max(welfare_gap),
    
    endpoint_welfare_order_holds =
      minimum_gap >= -1e-8,
    
    maximum_numerical_violation =
      max(
        -minimum_gap,
        0
      ),
    
    fraction_logistic_strictly_better =
      mean(
        welfare_gap >
          1e-8
      ),
    
    fraction_equal =
      mean(
        abs(welfare_gap) <=
          1e-8
      ),
    
    .groups = "drop"
  )


write.csv(
  endpoint_gap_summary,
  file.path(
    data_dir,
    "exp9_endpoint_welfare_dominance_summary.csv"
  ),
  row.names = FALSE
)


cat(
  "\nEndpoint welfare-dominance summary:\n"
)


print(
  endpoint_gap_summary,
  digits = 10,
  row.names = FALSE
)


if (
  any(
    !endpoint_gap_summary$
    endpoint_welfare_order_holds
  )
) {
  warning(
    paste0(
      "The computed endpoint logistic-minus-Laplace welfare gap ",
      "falls below the selected numerical tolerance."
    )
  )
}


endpoint_plot <- ggplot(
  endpoint_results,
  aes(
    x = prior_high,
    y = welfare_gap,
    color = mu_label,
    group = mu_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray35",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  scale_color_manual(
    values = c(
      "mu = 0.5" = "#1B9E77",
      "mu = 1" = "#D55E00",
      "mu = 2" = "#7570B3"
    ),
    name = NULL
  ) +
  labs(
    x = expression(
      "Prior probability of the high type  " *
        pi
    ),
    y = expression(
      "Endpoint welfare gap  " *
        (
          W[Log] -
            W[Lap]
        ) /
        W[FB]
    ),
    subtitle = paste0(
      "Common tight mu-GDP; logistic weakly dominates ",
      "in every binary endpoint decision problem"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  endpoint_plot <- endpoint_plot +
    labs(
      title = "Binary endpoint welfare comparison"
    )
}


save_figure(
  plot_object = endpoint_plot,
  filename_base = "fig_exp9_endpoint_welfare",
  width = 8.4,
  height = 5.2
)


# ================================================================
# 26. CONTINUOUS-TYPE WELFARE EXPERIMENT
# ================================================================

summary_rows <- list()
gap_rows <- list()
overflow_rows <- list()

summary_counter <- 0L
gap_counter <- 0L
overflow_counter <- 0L


for (epsilon_value in epsilon_grid) {
  implied_mu <- mu_from_epsilon_delta(
    epsilon = epsilon_value,
    delta = delta_fixed
  )
  
  cat(
    sprintf(
      paste0(
        "\nPrivacy index epsilon = %.3f, ",
        "delta = %.0e, implied mu = %.6f\n"
      ),
      epsilon_value,
      delta_fixed,
      implied_mu
    )
  )
  
  for (axis_name in privacy_axes) {
    axis_welfare_draws <- list()
    
    for (channel_name in channels) {
      scale_value <- calibrate_channel(
        axis_name = axis_name,
        channel_name = channel_name,
        epsilon_value = epsilon_value,
        delta_value = delta_fixed
      )
      
      if (!is.finite(scale_value)) {
        next
      }
      
      cat(
        sprintf(
          "  %-5s | %-18s | scale = %.6f\n",
          axis_name,
          channel_name,
          scale_value
        )
      )
      
      simulation_result <- continuous_welfare_draws(
        channel_name = channel_name,
        scale_value = scale_value,
        type_matrix = type_matrix,
        uniform_matrix = uniform_matrix,
        aggregate_true_type = aggregate_true_type,
        type_grid = type_grid,
        signal_grid_size = signal_grid_size
      )
      
      welfare_ratio <- ratio_estimate_and_se(
        numerator_draws =
          simulation_result$welfare_draws,
        
        denominator_draws =
          first_best_draws
      )
      
      axis_welfare_draws[[
        channel_name
      ]] <- simulation_result$welfare_draws
      
      summary_counter <- summary_counter + 1L
      
      summary_rows[[
        summary_counter
      ]] <- data.frame(
        axis =
          axis_name,
        
        channel =
          channel_name,
        
        epsilon =
          epsilon_value,
        
        delta =
          delta_fixed,
        
        mu =
          implied_mu,
        
        privacy_x =
          if (
            axis_name == "axis3"
          ) {
            implied_mu
          } else {
            epsilon_value
          },
        
        scale =
          scale_value,
        
        noise_sd =
          noise_standard_deviation(
            channel_name = channel_name,
            scale_value = scale_value
          ),
        
        pure_epsilon =
          pure_epsilon_index(
            channel_name = channel_name,
            scale_value = scale_value
          ),
        
        mu_index =
          tight_mu_index(
            channel_name = channel_name,
            scale_value = scale_value
          ),
        
        efficiency =
          unname(
            welfare_ratio["estimate"]
          ),
        
        efficiency_se =
          unname(
            welfare_ratio["standard_error"]
          ),
        
        build_probability =
          simulation_result$build_probability,
        
        mean_aggregate_posterior =
          simulation_result$
          mean_aggregate_posterior,
        
        overflow_rate =
          simulation_result$overflow_rate,
        
        stringsAsFactors = FALSE
      )
      
      overflow_counter <- overflow_counter + 1L
      
      overflow_rows[[
        overflow_counter
      ]] <- data.frame(
        axis =
          axis_name,
        
        channel =
          channel_name,
        
        epsilon =
          epsilon_value,
        
        mu =
          implied_mu,
        
        scale =
          scale_value,
        
        overflow_rate =
          simulation_result$overflow_rate,
        
        stringsAsFactors = FALSE
      )
    }
    
    if (
      all(
        c(
          "Laplace",
          "Logistic"
        ) %in%
        names(axis_welfare_draws)
      )
    ) {
      paired_numerator <-
        axis_welfare_draws[[
          "Laplace"
        ]] -
        axis_welfare_draws[[
          "Logistic"
        ]]
      
      paired_ratio <- ratio_estimate_and_se(
        numerator_draws =
          paired_numerator,
        
        denominator_draws =
          first_best_draws
      )
      
      gap_estimate <- unname(
        paired_ratio["estimate"]
      )
      
      gap_standard_error <- unname(
        paired_ratio["standard_error"]
      )
      
      gap_counter <- gap_counter + 1L
      
      gap_rows[[
        gap_counter
      ]] <- data.frame(
        axis =
          axis_name,
        
        epsilon =
          epsilon_value,
        
        delta =
          delta_fixed,
        
        mu =
          implied_mu,
        
        privacy_x =
          if (
            axis_name == "axis3"
          ) {
            implied_mu
          } else {
            epsilon_value
          },
        
        gap =
          gap_estimate,
        
        gap_se =
          gap_standard_error,
        
        ci_lower =
          gap_estimate -
          1.96 *
          gap_standard_error,
        
        ci_upper =
          gap_estimate +
          1.96 *
          gap_standard_error,
        
        stringsAsFactors = FALSE
      )
    }
  }
}


welfare_results <- bind_rows(
  summary_rows
) %>%
  mutate(
    channel = factor(
      channel,
      levels = channels
    ),
    
    axis = factor(
      axis,
      levels = privacy_axes
    )
  )


paired_gap_results <- bind_rows(
  gap_rows
) %>%
  mutate(
    axis = factor(
      axis,
      levels = privacy_axes
    )
  )


overflow_results <- bind_rows(
  overflow_rows
) %>%
  mutate(
    channel = factor(
      channel,
      levels = channels
    ),
    
    axis = factor(
      axis,
      levels = privacy_axes
    )
  )


write.csv(
  welfare_results,
  file.path(
    data_dir,
    "exp9_welfare_three_axes.csv"
  ),
  row.names = FALSE
)


write.csv(
  paired_gap_results,
  file.path(
    data_dir,
    "exp9_laplace_minus_logistic_paired_gap.csv"
  ),
  row.names = FALSE
)


write.csv(
  overflow_results,
  file.path(
    data_dir,
    "exp9_overflow_diagnostics.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 27. WELFARE PLOTS
# ================================================================

make_welfare_plot <- function(
    axis_name,
    x_axis_label,
    filename_base
) {
  plot_data <- welfare_results %>%
    filter(
      axis == axis_name
    )
  
  plot_object <- ggplot(
    plot_data,
    aes(
      x = privacy_x,
      y = efficiency,
      color = channel,
      group = channel
    )
  ) +
    geom_ribbon(
      aes(
        ymin =
          efficiency -
          1.96 *
          efficiency_se,
        
        ymax =
          efficiency +
          1.96 *
          efficiency_se,
        
        fill =
          channel
      ),
      alpha = 0.08,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(
      linewidth = 1.05
    ) +
    geom_point(
      size = 2.2
    ) +
    scale_color_manual(
      values = channel_colors,
      name = "Channel",
      drop = TRUE
    ) +
    scale_fill_manual(
      values = channel_colors,
      guide = "none",
      drop = TRUE
    ) +
    labs(
      x = x_axis_label,
      y = expression(
        "Welfare efficiency  " *
          W /
          W[FB]
      )
    ) +
    theme_paper()
  
  if (SHOW_TITLES) {
    plot_object <- plot_object +
      labs(
        title = axis_labels[[
          axis_name
        ]]
      )
  }
  
  save_figure(
    plot_object = plot_object,
    filename_base = filename_base,
    width = 8.8,
    height = 5
  )
  
  plot_object
}


plot_axis1 <- make_welfare_plot(
  axis_name = "axis1",
  x_axis_label = expression(
    "Pure privacy budget  " *
      epsilon
  ),
  filename_base =
    "fig_exp9_welfare_axis1_equal_epsilon"
)


plot_axis2 <- make_welfare_plot(
  axis_name = "axis2",
  x_axis_label = expression(
    "Approximate-DP budget  " *
      epsilon *
      "  at fixed  " *
      delta
  ),
  filename_base =
    "fig_exp9_welfare_axis2_eps_delta"
)


plot_axis3 <- make_welfare_plot(
  axis_name = "axis3",
  x_axis_label = expression(
    "Gaussian privacy parameter  " *
      mu
  ),
  filename_base =
    "fig_exp9_welfare_axis3_common_mu"
)


# ================================================================
# 28. PAIRED GAP PLOTS
# ================================================================

make_gap_plot <- function(
    axis_name,
    x_axis_label,
    filename_base
) {
  plot_data <- paired_gap_results %>%
    filter(
      axis == axis_name
    )
  
  plot_object <- ggplot(
    plot_data,
    aes(
      x = privacy_x,
      y = gap
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
        ymin = ci_lower,
        ymax = ci_upper
      ),
      alpha = 0.12,
      fill = "#D55E00"
    ) +
    geom_line(
      linewidth = 1.1,
      color = "#D55E00"
    ) +
    geom_point(
      size = 2.4,
      color = "#D55E00"
    ) +
    labs(
      x = x_axis_label,
      y = expression(
        (
          W[Lap] -
            W[Log]
        ) /
          W[FB]
      )
    ) +
    theme_paper()
  
  if (SHOW_TITLES) {
    plot_object <- plot_object +
      labs(
        title = paste(
          axis_labels[[
            axis_name
          ]],
          ": paired continuous-type gap"
        )
      )
  }
  
  save_figure(
    plot_object = plot_object,
    filename_base = filename_base,
    width = 8.8,
    height = 5
  )
  
  plot_object
}


gap_plot_axis1 <- make_gap_plot(
  axis_name = "axis1",
  x_axis_label = expression(
    "Pure privacy budget  " *
      epsilon
  ),
  filename_base =
    "fig_exp9_gap_axis1_equal_epsilon"
)


gap_plot_axis2 <- make_gap_plot(
  axis_name = "axis2",
  x_axis_label = expression(
    "Approximate-DP budget  " *
      epsilon *
      "  at fixed  " *
      delta
  ),
  filename_base =
    "fig_exp9_gap_axis2_eps_delta"
)


gap_plot_axis3 <- make_gap_plot(
  axis_name = "axis3",
  x_axis_label = expression(
    "Gaussian privacy parameter  " *
      mu
  ),
  filename_base =
    "fig_exp9_gap_axis3_common_mu"
)


# ================================================================
# 29. REPRESENTATIVE CALIBRATION TABLE
# ================================================================

representative_epsilon <- if (
  1 %in% epsilon_grid
) {
  1
} else {
  epsilon_grid[
    ceiling(
      length(epsilon_grid) / 2
    )
  ]
}


calibration_table <- welfare_results %>%
  filter(
    abs(
      epsilon -
        representative_epsilon
    ) <
      1e-10
  ) %>%
  transmute(
    axis,
    channel,
    
    epsilon =
      epsilon,
    
    delta =
      delta,
    
    mu =
      round(
        mu,
        6
      ),
    
    scale =
      round(
        scale,
        6
      ),
    
    noise_sd =
      round(
        noise_sd,
        6
      ),
    
    pure_epsilon =
      ifelse(
        is.finite(pure_epsilon),
        round(
          pure_epsilon,
          6
        ),
        NA_real_
      ),
    
    mu_index =
      round(
        mu_index,
        6
      ),
    
    efficiency =
      round(
        efficiency,
        6
      ),
    
    efficiency_se =
      round(
        efficiency_se,
        6
      ),
    
    build_probability =
      round(
        build_probability,
        6
      ),
    
    overflow_rate =
      signif(
        overflow_rate,
        4
      )
  ) %>%
  arrange(
    axis,
    channel
  )


write.csv(
  calibration_table,
  file.path(
    data_dir,
    sprintf(
      "exp9_calibration_table_eps_%g.csv",
      representative_epsilon
    )
  ),
  row.names = FALSE
)


cat(
  sprintf(
    paste0(
      "\nCalibration and continuous welfare at ",
      "epsilon = %g, delta = %.0e:\n"
    ),
    representative_epsilon,
    delta_fixed
  )
)


print(
  as.data.frame(
    calibration_table
  ),
  row.names = FALSE
)


# ================================================================
# 30. DIAGNOSTIC REPORTS
# ================================================================

cat(
  "\nPaired continuous-type Laplace-minus-logistic gaps:\n"
)


print(
  paired_gap_results %>%
    transmute(
      axis,
      
      epsilon =
        epsilon,
      
      mu =
        round(
          mu,
          6
        ),
      
      gap =
        round(
          gap,
          7
        ),
      
      gap_se =
        round(
          gap_se,
          7
        ),
      
      ci_lower =
        round(
          ci_lower,
          7
        ),
      
      ci_upper =
        round(
          ci_upper,
          7
        )
    ) %>%
    as.data.frame(),
  row.names = FALSE
)


cat(
  "\nMaximum signal-grid overflow rates:\n"
)


overflow_summary <- overflow_results %>%
  group_by(
    axis,
    channel
  ) %>%
  summarize(
    maximum_overflow =
      max(overflow_rate),
    
    .groups = "drop"
  )


print(
  as.data.frame(
    overflow_summary
  ),
  digits = 10,
  row.names = FALSE
)


maximum_overflow <- max(
  overflow_results$overflow_rate
)


if (
  maximum_overflow >
  1e-5
) {
  warning(
    paste0(
      "Maximum signal-grid overflow rate is ",
      format(
        maximum_overflow,
        scientific = TRUE
      ),
      ". Increase the tail range or grid size."
    )
  )
}


cat(
  "\nEndpoint welfare at mu = 1 for selected priors:\n"
)


print(
  endpoint_results %>%
    filter(
      abs(
        mu - 1
      ) <
        1e-10,
      
      prior_high %in%
        c(
          0.1,
          0.3,
          0.5,
          0.7,
          0.9
        )
    ) %>%
    transmute(
      prior_high,
      
      efficiency_laplace =
        round(
          efficiency_laplace,
          6
        ),
      
      efficiency_logistic =
        round(
          efficiency_logistic,
          6
        ),
      
      welfare_gap =
        round(
          welfare_gap,
          6
        ),
      
      alpha_laplace =
        round(
          alpha_laplace,
          6
        ),
      
      alpha_logistic =
        round(
          alpha_logistic,
          6
        ),
      
      candidate_laplace,
      candidate_logistic
    ) %>%
    as.data.frame(),
  row.names = FALSE
)


# ================================================================
# 31. FINAL INTERPRETATION
# ================================================================

cat(
  paste0(
    "\nInterpretation:\n",
    "  * At equal scale, Laplace Blackwell-dominates logistic on\n",
    "    the complete state space.\n",
    "  * Equal pure-epsilon calibration gives equal scales and\n",
    "    therefore inherits the full Laplace Blackwell ordering.\n",
    "  * At common tight mu-GDP, the calibrated logistic endpoint\n",
    "    trade-off curve lies weakly below the Laplace endpoint\n",
    "    trade-off curve.\n",
    "  * Therefore, the common-mu logistic endpoint binary\n",
    "    experiment Blackwell-dominates the Laplace endpoint\n",
    "    experiment.\n",
    "  * Logistic consequently yields weakly higher optimized\n",
    "    welfare in every binary endpoint decision problem.\n",
    "  * This endpoint result does not establish a Blackwell\n",
    "    ordering of the full continuous-type experiments.\n",
    "  * Continuous-type welfare comparisons remain numerical and\n",
    "    specification-specific.\n",
    "  * Common underlying uniforms are used only as a Monte Carlo\n",
    "    variance-reduction device.\n",
    "  * Ratio and paired-gap standard errors account for the\n",
    "    randomness of the first-best denominator.\n",
    "  * Posterior integration is log-stable and trapezoidal.\n",
    "  * Signal-grid overflow rates should be negligible.\n"
  )
)


# ================================================================
# 32. COMPLETION MESSAGE
# ================================================================

cat(
  "\nExperiment 9 completed.\n"
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