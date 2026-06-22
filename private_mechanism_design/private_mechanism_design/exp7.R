# ================================================================
# EXPERIMENT 7
# Effect of three privacy channels on the planner's objective
#
# Optimal Binary Mechanism under Locally Private Signals
#
# ALLOCATION RULE
#
#   q_lambda(Y)
#     =
#     1{
#       sum_i [
#         E[X_i | Y_i]
#         +
#         lambda E[J_i(X_i) | Y_i]
#       ] >= 0
#     }.
#
# PRIORS
#
# Welfare experiment:
#
#   X ~ Uniform[-1,1],
#   lambda = 0.
#
# Knife-edge experiment:
#
#   X ~ Uniform[-0.5,1.5],
#   E[X] = 0.5,
#   E[J(X)] = -0.5,
#   lambda* = 1.
#
# WELFARE
#
#   W
#     =
#     E[
#       (sum_i X_i) q_lambda(Y)
#     ].
#
# First-best welfare:
#
#   W_FB
#     =
#     E[
#       (sum_i X_i)_+
#     ].
#
# Welfare efficiency:
#
#   W / W_FB.
#
# SIGNED REALIZED MYERSON REVENUE
#
#   R
#     =
#     E[
#       (sum_i E[J_i(X_i) | Y_i])
#       q_lambda(Y)
#     ].
#
# Revenue is not truncated at zero. It may be negative.
#
# PRIVACY CALIBRATION
#
# All three channels are calibrated to the same
#
#   (epsilon, delta)-LDP
#
# frontier with delta = 10^(-6).
#
# Gaussian and clamped Gaussian use the analytic Gaussian profile.
# Laplace and logistic use their hockey-stick privacy profiles.
#
# NUMERICAL IMPROVEMENTS
#
#   * common random numbers across channels;
#   * adaptive signal grids;
#   * stable log-scale posterior quadrature;
#   * trapezoidal integration;
#   * paired ratio standard errors for W/W_FB;
#   * vectorized Monte Carlo evaluation;
#   * tower-property diagnostics.
#
# AXES
#
# All plot axes use ordinary linear coordinates.
# There is no logarithmic or local-coordinate transformation.
#
# OUTPUT
#
#   figures/fig_exp7_efficiency_vs_eps.pdf
#   figures/fig_exp7_revenue_vs_eps.pdf
#   figures/fig_exp7_objective_vs_lambda.pdf
#
#   tables/exp7_privacy_calibration.csv
#   tables/exp7_welfare_vs_eps.csv
#   tables/exp7_revenue_vs_eps.csv
#   tables/exp7_vs_lambda.csv
#   tables/exp7_sanity_report.csv
#
# ================================================================


rm(list = ls())


# ================================================================
# 0. USER OPTIONS
# ================================================================

FAST        <- FALSE
SHOW_TITLES <- FALSE

set.seed(20260610)


out_dir <- paste0(
  "/Users/behroozmoosavi/Desktop/privacy/codes/",
  "R_codes/experiment/exp7"
)

fig_dir <- file.path(
  out_dir,
  "figures"
)

data_dir <- file.path(
  out_dir,
  "tables"
)


for (directory in c(
  out_dir,
  fig_dir,
  data_dir
)) {
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

theme_paper <- function(
    base_size = 13
) {
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
  
  invisible(
    output_file
  )
}


# ================================================================
# 3. NUMERICAL HELPERS
# ================================================================

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


safe_probability <- function(
    probability,
    tolerance = 1e-14
) {
  pmin(
    pmax(
      probability,
      tolerance
    ),
    1 -
      tolerance
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
    
    scale *
      log(
        2 *
          probability
      ),
    
    -scale *
      log(
        2 *
          (
            1 -
              probability
          )
      )
  )
}


ratio_delta_standard_error <- function(
    numerator_draws,
    denominator_draws
) {
  if (
    length(
      numerator_draws
    ) !=
    length(
      denominator_draws
    )
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
    !is.finite(
      denominator_mean
    ) ||
    abs(
      denominator_mean
    ) <=
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
# 4. COMMON PRIVACY PARAMETERS
# ================================================================

privacy_delta <- 1e-5

support_width <- 2

clamp_fraction <- 0.20


clamped_sensitivity <- (
  1 -
    2 *
    clamp_fraction
) *
  support_width


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


epsilon_grid <- if (FAST) {
  c(
    0.5,
    1,
    2,
    4
  )
} else {
  c(
    0.25,
    0.5,
    0.75,
    1,
    1.5,
    2,
    3,
    4,
    6
  )
}


# ================================================================
# 5. ANALYTIC GAUSSIAN PRIVACY PROFILE
#
# For
#
#   mu = Delta/sigma,
#
# the exact Gaussian hockey-stick profile is
#
#   delta_G(epsilon,mu)
#     =
#     Phi(mu/2 - epsilon/mu)
#     -
#     exp(epsilon)
#     Phi(-mu/2 - epsilon/mu).
#
# ================================================================

gaussian_delta_profile <- function(
    epsilon,
    mu
) {
  if (mu <= 0) {
    return(
      0
    )
  }
  
  value <- pnorm(
    mu /
      2 -
      epsilon /
      mu
  ) -
    exp(
      epsilon
    ) *
    pnorm(
      -mu /
        2 -
        epsilon /
        mu
    )
  
  min(
    max(
      value,
      0
    ),
    1
  )
}


solve_gaussian_mu <- function(
    epsilon,
    delta
) {
  objective <- function(mu) {
    gaussian_delta_profile(
      epsilon = epsilon,
      mu = mu
    ) -
      delta
  }
  
  lower <- 1e-8
  upper <- 1
  
  while (
    objective(
      upper
    ) <
    0
  ) {
    upper <- 2 *
      upper
    
    if (upper > 200) {
      stop(
        paste0(
          "Unable to bracket Gaussian privacy parameter at epsilon = ",
          epsilon,
          "."
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
    tol = 1e-10
  )$root
}


# ================================================================
# 6. LAPLACE AND LOGISTIC HOCKEY-STICK PROFILES
#
# For standardized density f and separation r,
#
#   delta(epsilon,r)
#     =
#     integral [
#       f(z) - exp(epsilon) f(z-r)
#     ]_+ dz.
#
# Here
#
#   r = Delta/scale.
#
# ================================================================

standard_log_density <- function(
    z,
    family
) {
  switch(
    family,
    
    "Laplace" =
      -log(
        2
      ) -
      abs(
        z
      ),
    
    "Logistic" =
      dlogis(
        z,
        location = 0,
        scale = 1,
        log = TRUE
      ),
    
    stop(
      paste(
        "Unknown location family:",
        family
      )
    )
  )
}


hockey_stick_integrand <- function(
    z,
    epsilon,
    separation,
    family
) {
  log_p <- standard_log_density(
    z,
    family
  )
  
  log_q <- standard_log_density(
    z -
      separation,
    family
  )
  
  maximum_log <- pmax(
    log_p,
    epsilon +
      log_q
  )
  
  scaled_difference <- exp(
    log_p -
      maximum_log
  ) -
    exp(
      epsilon +
        log_q -
        maximum_log
    )
  
  exp(
    maximum_log
  ) *
    pmax(
      scaled_difference,
      0
    )
}


standard_delta_profile <- function(
    epsilon,
    separation,
    family
) {
  if (separation <= 0) {
    return(
      0
    )
  }
  
  # For these standardized location families, separation <= epsilon
  # already satisfies pure epsilon-LDP.
  
  if (separation <= epsilon) {
    return(
      0
    )
  }
  
  integration_bound <- max(
    40,
    separation +
      30
  )
  
  integral_result <- integrate(
    f = function(z) {
      hockey_stick_integrand(
        z = z,
        epsilon = epsilon,
        separation = separation,
        family = family
      )
    },
    lower = -integration_bound,
    upper = integration_bound,
    subdivisions = 4000L,
    rel.tol = 1e-9,
    abs.tol = 1e-12,
    stop.on.error = TRUE
  )
  
  min(
    max(
      integral_result$value,
      0
    ),
    1
  )
}


solve_standardized_separation <- function(
    epsilon,
    delta,
    family
) {
  objective <- function(
    separation
  ) {
    standard_delta_profile(
      epsilon = epsilon,
      separation = separation,
      family = family
    ) -
      delta
  }
  
  lower <- epsilon
  
  upper <- max(
    epsilon +
      1,
    2 *
      epsilon +
      1
  )
  
  while (
    objective(
      upper
    ) <
    0
  ) {
    upper <- 2 *
      upper
    
    if (upper > 200) {
      stop(
        paste0(
          "Unable to bracket ",
          family,
          " privacy separation at epsilon = ",
          epsilon,
          "."
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
# 7. PRIVACY CALIBRATION TABLE
# ================================================================

calibration_rows <- lapply(
  epsilon_grid,
  function(epsilon_value) {
    cat(
      "Calibrating epsilon =",
      epsilon_value,
      "...\n"
    )
    
    gaussian_mu <- solve_gaussian_mu(
      epsilon = epsilon_value,
      delta = privacy_delta
    )
    
    laplace_separation <-
      solve_standardized_separation(
        epsilon = epsilon_value,
        delta = privacy_delta,
        family = "Laplace"
      )
    
    logistic_separation <-
      solve_standardized_separation(
        epsilon = epsilon_value,
        delta = privacy_delta,
        family = "Logistic"
      )
    
    data.frame(
      epsilon = epsilon_value,
      delta = privacy_delta,
      
      gaussian_mu =
        gaussian_mu,
      
      laplace_separation =
        laplace_separation,
      
      logistic_separation =
        logistic_separation,
      
      gaussian_sigma =
        support_width /
        gaussian_mu,
      
      clamped_gaussian_sigma =
        clamped_sensitivity /
        gaussian_mu,
      
      laplace_scale =
        support_width /
        laplace_separation,
      
      logistic_scale =
        support_width /
        logistic_separation,
      
      stringsAsFactors = FALSE
    )
  }
)


calibration_table <- bind_rows(
  calibration_rows
)


write.csv(
  calibration_table,
  file.path(
    data_dir,
    "exp7_privacy_calibration.csv"
  ),
  row.names = FALSE
)


cat(
  "\nPrivacy calibration table:\n"
)


print(
  calibration_table,
  digits = 6,
  row.names = FALSE
)


get_calibration_row <- function(
    epsilon_value
) {
  selected_row <- calibration_table[
    abs(
      calibration_table$epsilon -
        epsilon_value
    ) <
      1e-12,
    ,
    drop = FALSE
  ]
  
  if (nrow(selected_row) != 1L) {
    stop(
      paste(
        "Expected exactly one calibration row at epsilon =",
        epsilon_value
      )
    )
  }
  
  selected_row
}


# ================================================================
# 8. PRIOR REGIMES
# ================================================================

make_regime <- function(
    lower,
    upper
) {
  width <- upper -
    lower
  
  clamp_lower <- lower +
    clamp_fraction *
    width
  
  clamp_upper <- upper -
    clamp_fraction *
    width
  
  list(
    lower = lower,
    upper = upper,
    width = width,
    clamp_lower = clamp_lower,
    clamp_upper = clamp_upper
  )
}


welfare_regime <- make_regime(
  lower = -1,
  upper = 1
)


knife_edge_regime <- make_regime(
  lower = -0.5,
  upper = 1.5
)


# ================================================================
# 9. CHANNEL SCALE
# ================================================================

channel_scale <- function(
    channel_name,
    epsilon_value
) {
  calibration <- get_calibration_row(
    epsilon_value
  )
  
  switch(
    channel_name,
    
    "Gaussian" =
      calibration$gaussian_sigma[[1]],
    
    "Laplace" =
      calibration$laplace_scale[[1]],
    
    "Logistic" =
      calibration$logistic_scale[[1]],
    
    "Clamped Gaussian" =
      calibration$clamped_gaussian_sigma[[1]],
    
    stop(
      paste(
        "Unknown channel:",
        channel_name
      )
    )
  )
}


transformed_type_canonical <- function(
    z,
    channel_name
) {
  if (
    channel_name ==
    "Clamped Gaussian"
  ) {
    clamp_function(
      z,
      clamp_fraction *
        support_width,
      (
        1 -
          clamp_fraction
      ) *
        support_width
    )
  } else {
    z
  }
}


# ================================================================
# 10. CHANNEL SIMULATION USING COMMON UNIFORMS
# ================================================================

simulate_channel_from_uniforms <- function(
    x,
    uniform_noise,
    channel_name,
    epsilon_value,
    regime
) {
  uniform_noise <- safe_probability(
    uniform_noise
  )
  
  scale_value <- channel_scale(
    channel_name = channel_name,
    epsilon_value = epsilon_value
  )
  
  if (
    channel_name ==
    "Clamped Gaussian"
  ) {
    transformed_x <- clamp_function(
      x,
      regime$clamp_lower,
      regime$clamp_upper
    )
  } else {
    transformed_x <- x
  }
  
  switch(
    channel_name,
    
    "Gaussian" =
      transformed_x +
      qnorm(
        uniform_noise,
        mean = 0,
        sd = scale_value
      ),
    
    "Laplace" =
      transformed_x +
      laplace_quantile(
        probability = uniform_noise,
        scale = scale_value
      ),
    
    "Logistic" =
      transformed_x +
      qlogis(
        uniform_noise,
        location = 0,
        scale = scale_value
      ),
    
    "Clamped Gaussian" =
      transformed_x +
      qnorm(
        uniform_noise,
        mean = 0,
        sd = scale_value
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
# 11. ADAPTIVE CANONICAL SIGNAL GRIDS
# ================================================================

make_signal_grid <- function(
    channel_name,
    epsilon_value
) {
  scale_value <- channel_scale(
    channel_name = channel_name,
    epsilon_value = epsilon_value
  )
  
  transformed_lower <- transformed_type_canonical(
    z = 0,
    channel_name = channel_name
  )
  
  transformed_upper <- transformed_type_canonical(
    z = support_width,
    channel_name = channel_name
  )
  
  tail_multiplier <- switch(
    channel_name,
    
    "Gaussian" =
      9,
    
    "Clamped Gaussian" =
      9,
    
    "Laplace" =
      18,
    
    "Logistic" =
      20,
    
    stop(
      paste(
        "Unknown channel:",
        channel_name
      )
    )
  )
  
  seq(
    transformed_lower -
      tail_multiplier *
      scale_value,
    
    transformed_upper +
      tail_multiplier *
      scale_value,
    
    length.out = if (FAST) {
      2201
    } else {
      4201
    }
  )
}


# ================================================================
# 12. CONDITIONAL LOG-DENSITIES
# ================================================================

channel_log_density <- function(
    y,
    z,
    channel_name,
    epsilon_value
) {
  scale_value <- channel_scale(
    channel_name = channel_name,
    epsilon_value = epsilon_value
  )
  
  transformed_z <- transformed_type_canonical(
    z = z,
    channel_name = channel_name
  )
  
  switch(
    channel_name,
    
    "Gaussian" =
      dnorm(
        y,
        mean = transformed_z,
        sd = scale_value,
        log = TRUE
      ),
    
    "Laplace" =
      -log(
        2 *
          scale_value
      ) -
      abs(
        y -
          transformed_z
      ) /
      scale_value,
    
    "Logistic" =
      dlogis(
        y -
          transformed_z,
        location = 0,
        scale = scale_value,
        log = TRUE
      ),
    
    "Clamped Gaussian" =
      dnorm(
        y,
        mean = transformed_z,
        sd = scale_value,
        log = TRUE
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
# 13. STABLE CANONICAL POSTERIOR MOMENTS
#
# Let
#
#   Z ~ Uniform[0,2].
#
# For a regime X = a + Z,
#
#   E[X | Y]
#     =
#     a + E[Z | Y-a].
#
# For X ~ Uniform[a,a+2],
#
#   J_X(X)
#     =
#     2X-(a+2)
#     =
#     a + 2Z-2.
#
# Hence
#
#   E[J_X(X) | Y]
#     =
#     a + E[2Z-2 | Y-a].
#
# ================================================================

posterior_moment_table <- function(
    channel_name,
    epsilon_value
) {
  type_grid <- seq(
    0,
    support_width,
    length.out = if (FAST) {
      1001
    } else {
      2001
    }
  )
  
  type_step <- type_grid[2] -
    type_grid[1]
  
  trapezoid_weights <- rep(
    1,
    length(
      type_grid
    )
  )
  
  trapezoid_weights[
    c(
      1,
      length(
        trapezoid_weights
      )
    )
  ] <- 0.5
  
  signal_grid <- make_signal_grid(
    channel_name = channel_name,
    epsilon_value = epsilon_value
  )
  
  prior_density <- rep(
    1 /
      support_width,
    length(
      type_grid
    )
  )
  
  canonical_virtual_value <-
    2 *
    type_grid -
    support_width
  
  rows <- lapply(
    signal_grid,
    function(signal_value) {
      log_weights <- channel_log_density(
        y = signal_value,
        z = type_grid,
        channel_name = channel_name,
        epsilon_value = epsilon_value
      ) +
        log(
          prior_density
        ) +
        log(
          trapezoid_weights
        )
      
      maximum_log_weight <- max(
        log_weights
      )
      
      if (!is.finite(maximum_log_weight)) {
        return(
          c(
            zhat = NA_real_,
            Jhat0 = NA_real_
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
        type_step
      
      if (
        !is.finite(
          denominator
        ) ||
        denominator <=
        .Machine$double.eps
      ) {
        return(
          c(
            zhat = NA_real_,
            Jhat0 = NA_real_
          )
        )
      }
      
      posterior_z <- sum(
        type_grid *
          stabilized_weights
      ) *
        type_step /
        denominator
      
      posterior_J0 <- sum(
        canonical_virtual_value *
          stabilized_weights
      ) *
        type_step /
        denominator
      
      c(
        zhat = posterior_z,
        Jhat0 = posterior_J0
      )
    }
  )
  
  moment_matrix <- do.call(
    rbind,
    rows
  )
  
  list(
    signal_grid =
      signal_grid,
    
    zhat =
      moment_matrix[
        ,
        "zhat"
      ],
    
    Jhat0 =
      moment_matrix[
        ,
        "Jhat0"
      ]
  )
}


make_interpolator <- function(
    grid,
    values
) {
  force(
    grid
  )
  
  force(
    values
  )
  
  function(signal_values) {
    approx(
      x = grid,
      y = values,
      xout = signal_values,
      rule = 2,
      ties = "ordered"
    )$y
  }
}


# ================================================================
# 14. POSTERIOR CACHE
# ================================================================

posterior_cache <- new.env(
  parent = emptyenv()
)


posterior_cache_key <- function(
    channel_name,
    epsilon_value
) {
  paste(
    channel_name,
    format(
      epsilon_value,
      scientific = FALSE,
      trim = TRUE
    ),
    sep = "__"
  )
}


get_posterior_functions <- function(
    channel_name,
    epsilon_value
) {
  key <- posterior_cache_key(
    channel_name = channel_name,
    epsilon_value = epsilon_value
  )
  
  if (!exists(
    key,
    envir = posterior_cache,
    inherits = FALSE
  )) {
    cat(
      "Building posterior table:",
      channel_name,
      "| epsilon =",
      epsilon_value,
      "\n"
    )
    
    posterior_table <- posterior_moment_table(
      channel_name = channel_name,
      epsilon_value = epsilon_value
    )
    
    posterior_object <- list(
      zhat = make_interpolator(
        grid = posterior_table$signal_grid,
        values = posterior_table$zhat
      ),
      
      Jhat0 = make_interpolator(
        grid = posterior_table$signal_grid,
        values = posterior_table$Jhat0
      ),
      
      signal_lower = min(
        posterior_table$signal_grid
      ),
      
      signal_upper = max(
        posterior_table$signal_grid
      )
    )
    
    assign(
      key,
      posterior_object,
      envir = posterior_cache
    )
  }
  
  get(
    key,
    envir = posterior_cache,
    inherits = FALSE
  )
}


# ================================================================
# 15. SHARED MONTE CARLO DRAWS
# ================================================================

generate_regime_draws <- function(
    regime,
    replications,
    number_of_agents
) {
  type_matrix <- matrix(
    runif(
      replications *
        number_of_agents,
      min = regime$lower,
      max = regime$upper
    ),
    nrow = replications,
    ncol = number_of_agents
  )
  
  uniform_noise_matrix <- matrix(
    runif(
      replications *
        number_of_agents
    ),
    nrow = replications,
    ncol = number_of_agents
  )
  
  aggregate_type <- rowSums(
    type_matrix
  )
  
  first_best_draws <- pmax(
    aggregate_type,
    0
  )
  
  list(
    type_matrix =
      type_matrix,
    
    uniform_noise_matrix =
      uniform_noise_matrix,
    
    aggregate_type =
      aggregate_type,
    
    first_best_draws =
      first_best_draws,
    
    first_best_mean =
      mean(
        first_best_draws
      )
  )
}


# ================================================================
# 16. PRECOMPUTE ONE CHANNEL-EPSILON CELL
# ================================================================

compute_channel_statistics <- function(
    regime,
    channel_name,
    epsilon_value,
    shared_draws
) {
  posterior_functions <- get_posterior_functions(
    channel_name = channel_name,
    epsilon_value = epsilon_value
  )
  
  signal_matrix <- simulate_channel_from_uniforms(
    x = shared_draws$type_matrix,
    uniform_noise = shared_draws$uniform_noise_matrix,
    channel_name = channel_name,
    epsilon_value = epsilon_value,
    regime = regime
  )
  
  signal_vector <- as.vector(
    signal_matrix
  )
  
  canonical_signal_vector <-
    signal_vector -
    regime$lower
  
  zhat_vector <- posterior_functions$zhat(
    canonical_signal_vector
  )
  
  Jhat0_vector <- posterior_functions$Jhat0(
    canonical_signal_vector
  )
  
  posterior_x_vector <-
    regime$lower +
    zhat_vector
  
  posterior_J_vector <-
    regime$lower +
    Jhat0_vector
  
  if (
    any(
      !is.finite(
        posterior_x_vector
      )
    ) ||
    any(
      !is.finite(
        posterior_J_vector
      )
    )
  ) {
    stop(
      paste0(
        "Non-finite posterior moments for channel ",
        channel_name,
        " at epsilon = ",
        epsilon_value,
        "."
      )
    )
  }
  
  number_of_replications <- nrow(
    shared_draws$type_matrix
  )
  
  number_of_agents <- ncol(
    shared_draws$type_matrix
  )
  
  posterior_x_matrix <- matrix(
    posterior_x_vector,
    nrow = number_of_replications,
    ncol = number_of_agents
  )
  
  posterior_J_matrix <- matrix(
    posterior_J_vector,
    nrow = number_of_replications,
    ncol = number_of_agents
  )
  
  list(
    aggregate_xhat =
      rowSums(
        posterior_x_matrix
      ),
    
    aggregate_Jhat =
      rowSums(
        posterior_J_matrix
      ),
    
    overflow_rate =
      mean(
        canonical_signal_vector <
          posterior_functions$signal_lower |
          canonical_signal_vector >
          posterior_functions$signal_upper
      )
  )
}


# ================================================================
# 17. EVALUATE ALLOCATION RULE
# ================================================================

evaluate_rule <- function(
    channel_name,
    epsilon_value,
    lambda_value,
    shared_draws,
    channel_statistics
) {
  aggregate_score <-
    channel_statistics$aggregate_xhat +
    lambda_value *
    channel_statistics$aggregate_Jhat
  
  build_indicator <- as.numeric(
    aggregate_score >= 0
  )
  
  welfare_draws <-
    shared_draws$aggregate_type *
    build_indicator
  
  # Signed realized Myerson revenue.
  # Do not truncate aggregate_Jhat at zero.
  
  revenue_draws <-
    channel_statistics$aggregate_Jhat *
    build_indicator
  
  efficiency_estimate <- mean(
    welfare_draws
  ) /
    shared_draws$first_best_mean
  
  efficiency_standard_error <-
    ratio_delta_standard_error(
      numerator_draws =
        welfare_draws,
      
      denominator_draws =
        shared_draws$first_best_draws
    )
  
  revenue_estimate <- mean(
    revenue_draws
  )
  
  revenue_standard_error <- sd(
    revenue_draws
  ) /
    sqrt(
      length(
        revenue_draws
      )
    )
  
  data.frame(
    channel =
      channel_name,
    
    epsilon =
      epsilon_value,
    
    delta =
      privacy_delta,
    
    lambda =
      lambda_value,
    
    efficiency =
      efficiency_estimate,
    
    efficiency_se =
      efficiency_standard_error,
    
    efficiency_lower =
      efficiency_estimate -
      1.96 *
      efficiency_standard_error,
    
    efficiency_upper =
      efficiency_estimate +
      1.96 *
      efficiency_standard_error,
    
    revenue =
      revenue_estimate,
    
    revenue_se =
      revenue_standard_error,
    
    revenue_lower =
      revenue_estimate -
      1.96 *
      revenue_standard_error,
    
    revenue_upper =
      revenue_estimate +
      1.96 *
      revenue_standard_error,
    
    build_probability =
      mean(
        build_indicator
      ),
    
    mean_aggregate_xhat =
      mean(
        channel_statistics$aggregate_xhat
      ),
    
    mean_aggregate_Jhat =
      mean(
        channel_statistics$aggregate_Jhat
      ),
    
    signal_outside_grid_rate =
      channel_statistics$overflow_rate,
    
    stringsAsFactors = FALSE
  )
}


# ================================================================
# 18. MONTE CARLO SETTINGS
# ================================================================

number_of_agents <- if (FAST) {
  50
} else {
  100
}


monte_carlo_replications <- if (FAST) {
  4000
} else {
  20000
}


welfare_draws <- generate_regime_draws(
  regime = welfare_regime,
  replications = monte_carlo_replications,
  number_of_agents = number_of_agents
)


knife_edge_draws <- generate_regime_draws(
  regime = knife_edge_regime,
  replications = monte_carlo_replications,
  number_of_agents = number_of_agents
)


cat(
  sprintf(
    paste0(
      "\nFirst-best welfare means:\n",
      "  Uniform[-1,1]     : %.6f\n",
      "  Uniform[-0.5,1.5] : %.6f\n\n"
    ),
    welfare_draws$first_best_mean,
    knife_edge_draws$first_best_mean
  )
)


# ================================================================
# 19. PRECOMPUTE POSTERIOR AGGREGATES
# ================================================================

welfare_statistics <- list()

knife_edge_statistics <- list()


cell_key <- function(
    channel_name,
    epsilon_value
) {
  paste(
    channel_name,
    epsilon_value,
    sep = "__"
  )
}


for (epsilon_value in epsilon_grid) {
  for (channel_name in channels) {
    key <- cell_key(
      channel_name = channel_name,
      epsilon_value = epsilon_value
    )
    
    cat(
      "Computing welfare cell:",
      channel_name,
      "| epsilon =",
      epsilon_value,
      "\n"
    )
    
    welfare_statistics[[
      key
    ]] <- compute_channel_statistics(
      regime = welfare_regime,
      channel_name = channel_name,
      epsilon_value = epsilon_value,
      shared_draws = welfare_draws
    )
    
    cat(
      "Computing knife-edge cell:",
      channel_name,
      "| epsilon =",
      epsilon_value,
      "\n"
    )
    
    knife_edge_statistics[[
      key
    ]] <- compute_channel_statistics(
      regime = knife_edge_regime,
      channel_name = channel_name,
      epsilon_value = epsilon_value,
      shared_draws = knife_edge_draws
    )
  }
}


# ================================================================
# 20. WELFARE EFFICIENCY AGAINST epsilon
# ================================================================

welfare_rows <- list()

welfare_counter <- 0L


for (epsilon_value in epsilon_grid) {
  for (channel_name in channels) {
    welfare_counter <- welfare_counter +
      1L
    
    key <- cell_key(
      channel_name = channel_name,
      epsilon_value = epsilon_value
    )
    
    welfare_rows[[
      welfare_counter
    ]] <- evaluate_rule(
      channel_name = channel_name,
      epsilon_value = epsilon_value,
      lambda_value = 0,
      shared_draws = welfare_draws,
      channel_statistics =
        welfare_statistics[[
          key
        ]]
    )
  }
}


welfare_results <- bind_rows(
  welfare_rows
) %>%
  mutate(
    channel = factor(
      channel,
      levels = channels
    )
  )


write.csv(
  welfare_results,
  file.path(
    data_dir,
    "exp7_welfare_vs_eps.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 21. SIGNED REALIZED REVENUE AGAINST epsilon
#
# Knife-edge prior:
#
#   X ~ Uniform[-0.5,1.5]
#
# and
#
#   lambda = 1.
#
# Revenue can be negative.
# ================================================================

revenue_rows <- list()

revenue_counter <- 0L


for (epsilon_value in epsilon_grid) {
  for (channel_name in channels) {
    revenue_counter <- revenue_counter +
      1L
    
    key <- cell_key(
      channel_name = channel_name,
      epsilon_value = epsilon_value
    )
    
    revenue_rows[[
      revenue_counter
    ]] <- evaluate_rule(
      channel_name = channel_name,
      epsilon_value = epsilon_value,
      lambda_value = 1,
      shared_draws = knife_edge_draws,
      channel_statistics =
        knife_edge_statistics[[
          key
        ]]
    )
  }
}


revenue_results <- bind_rows(
  revenue_rows
) %>%
  mutate(
    channel = factor(
      channel,
      levels = channels
    )
  )


write.csv(
  revenue_results,
  file.path(
    data_dir,
    "exp7_revenue_vs_eps.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 22. WELFARE AND REVENUE AGAINST lambda AT epsilon = 1
# ================================================================

epsilon_fixed <- 1


lambda_grid <- if (FAST) {
  c(
    0,
    0.5,
    1,
    1.5,
    2
  )
} else {
  c(
    0,
    0.3,
    0.6,
    1,
    1.3,
    1.6,
    2
  )
}


lambda_rows <- list()

lambda_counter <- 0L


for (lambda_value in lambda_grid) {
  for (channel_name in channels) {
    lambda_counter <- lambda_counter +
      1L
    
    key <- cell_key(
      channel_name = channel_name,
      epsilon_value = epsilon_fixed
    )
    
    lambda_rows[[
      lambda_counter
    ]] <- evaluate_rule(
      channel_name = channel_name,
      epsilon_value = epsilon_fixed,
      lambda_value = lambda_value,
      shared_draws = knife_edge_draws,
      channel_statistics =
        knife_edge_statistics[[
          key
        ]]
    )
  }
}


lambda_results <- bind_rows(
  lambda_rows
) %>%
  mutate(
    channel = factor(
      channel,
      levels = channels
    )
  )


write.csv(
  lambda_results,
  file.path(
    data_dir,
    "exp7_vs_lambda.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 23. FIGURE 1
# WELFARE EFFICIENCY VERSUS epsilon
#
# Ordinary linear x- and y-axes.
# ================================================================

welfare_plot <- ggplot(
  welfare_results,
  aes(
    x = epsilon,
    y = efficiency,
    color = channel,
    group = channel
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray35",
    linewidth = 0.50
  ) +
  geom_ribbon(
    aes(
      ymin = efficiency_lower,
      ymax = efficiency_upper,
      fill = channel
    ),
    alpha = 0.08,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    size = 2.4
  ) +
  scale_color_manual(
    values = channel_colors,
    name = "Channel"
  ) +
  scale_fill_manual(
    values = channel_colors
  ) +
  scale_x_continuous(
    breaks = epsilon_grid,
    labels = label_number(
      accuracy = 0.01
    )
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(
      n = 8
    ),
    labels = label_number(
      accuracy = 0.01
    )
  ) +
  labs(
    x = expression(
      "privacy budget  " *
        epsilon
    ),
    y = expression(
      "welfare efficiency  " *
        W /
        W[FB]
    ),
    subtitle = paste0(
      "Common (epsilon, delta)-LDP frontier with ",
      "delta = 10^-6"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  welfare_plot <- welfare_plot +
    labs(
      title = "Welfare efficiency across privacy channels"
    )
}


save_figure(
  plot_object = welfare_plot,
  filename_base = "fig_exp7_efficiency_vs_eps",
  width = 7.6,
  height = 5.4
)


# ================================================================
# 24. FIGURE 2
# SIGNED REALIZED REVENUE VERSUS epsilon
#
# Ordinary linear x- and y-axes.
# Revenue is allowed to be negative.
# ================================================================

revenue_plot <- ggplot(
  revenue_results,
  aes(
    x = epsilon,
    y = revenue,
    color = channel,
    group = channel
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray35",
    linewidth = 0.45
  ) +
  geom_ribbon(
    aes(
      ymin = revenue_lower,
      ymax = revenue_upper,
      fill = channel
    ),
    alpha = 0.08,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    size = 2.4
  ) +
  scale_color_manual(
    values = channel_colors,
    name = "Channel"
  ) +
  scale_fill_manual(
    values = channel_colors
  ) +
  scale_x_continuous(
    breaks = epsilon_grid,
    labels = label_number(
      accuracy = 0.01
    )
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(
      n = 7
    ),
    labels = label_number(
      accuracy = 0.1
    )
  ) +
  labs(
    x = expression(
      "privacy budget  " *
        epsilon
    ),
    y = expression(
      "signed realized revenue  " *
        R
    ),
    subtitle = paste0(
      "lambda = 1; common (epsilon, delta)-LDP frontier"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  revenue_plot <- revenue_plot +
    labs(
      title = "Signed realized revenue at the knife edge"
    )
}


save_figure(
  plot_object = revenue_plot,
  filename_base = "fig_exp7_revenue_vs_eps",
  width = 7.6,
  height = 5.4
)


# ================================================================
# 25. FIGURE 3
# WELFARE AND SIGNED REVENUE VERSUS lambda
#
# Lambda uses an ordinary linear axis.
#
# The two panels use separate linear y-axis ranges because welfare
# efficiency and revenue have different units.
# ================================================================

lambda_long <- bind_rows(
  lambda_results %>%
    transmute(
      channel,
      lambda,
      outcome =
        "Welfare efficiency",
      value =
        efficiency,
      lower =
        efficiency_lower,
      upper =
        efficiency_upper
    ),
  
  lambda_results %>%
    transmute(
      channel,
      lambda,
      outcome =
        "Signed realized revenue",
      value =
        revenue,
      lower =
        revenue_lower,
      upper =
        revenue_upper
    )
) %>%
  mutate(
    channel = factor(
      channel,
      levels = channels
    ),
    
    outcome = factor(
      outcome,
      levels = c(
        "Welfare efficiency",
        "Signed realized revenue"
      )
    )
  )


lambda_plot <- ggplot(
  lambda_long,
  aes(
    x = lambda,
    y = value,
    color = channel,
    group = channel
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dotted",
    color = "gray35",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray75",
    linewidth = 0.35
  ) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper,
      fill = channel
    ),
    alpha = 0.07,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    size = 2.3
  ) +
  facet_wrap(
    ~ outcome,
    scales = "free_y",
    ncol = 2
  ) +
  scale_color_manual(
    values = channel_colors,
    name = "Channel"
  ) +
  scale_fill_manual(
    values = channel_colors
  ) +
  scale_x_continuous(
    breaks = lambda_grid,
    labels = label_number(
      accuracy = 0.1
    )
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(
      n = 7
    ),
    labels = label_number(
      accuracy = 0.01
    )
  ) +
  labs(
    x = expression(
      "revenue multiplier  " *
        lambda
    ),
    y = NULL,
    subtitle = paste0(
      "epsilon = 1, delta = 10^-6; ",
      "dotted line: lambda* = 1"
    )
  ) +
  theme_paper()


if (SHOW_TITLES) {
  lambda_plot <- lambda_plot +
    labs(
      title = "Welfare and signed revenue around the knife edge"
    )
}


save_figure(
  plot_object = lambda_plot,
  filename_base = "fig_exp7_objective_vs_lambda",
  width = 8.4,
  height = 5.0
)


# ================================================================
# 26. SANITY CHECKS
# ================================================================

tower_check <- bind_rows(
  welfare_results %>%
    transmute(
      panel =
        "welfare",
      
      channel,
      epsilon,
      lambda,
      
      mean_aggregate_xhat,
      mean_aggregate_Jhat,
      
      expected_total_x =
        number_of_agents *
        (
          welfare_regime$lower +
            welfare_regime$upper
        ) /
        2,
      
      expected_total_J =
        number_of_agents *
        welfare_regime$lower,
      
      signal_outside_grid_rate
    ),
  
  revenue_results %>%
    transmute(
      panel =
        "revenue",
      
      channel,
      epsilon,
      lambda,
      
      mean_aggregate_xhat,
      mean_aggregate_Jhat,
      
      expected_total_x =
        number_of_agents *
        (
          knife_edge_regime$lower +
            knife_edge_regime$upper
        ) /
        2,
      
      expected_total_J =
        number_of_agents *
        knife_edge_regime$lower,
      
      signal_outside_grid_rate
    )
) %>%
  mutate(
    aggregate_xhat_error =
      mean_aggregate_xhat -
      expected_total_x,
    
    aggregate_Jhat_error =
      mean_aggregate_Jhat -
      expected_total_J
  )


write.csv(
  tower_check,
  file.path(
    data_dir,
    "exp7_sanity_report.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT 7 SANITY REPORT\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Welfare ranking at epsilon = 1, lambda = 0:\n"
)


print(
  welfare_results %>%
    filter(
      abs(
        epsilon -
          1
      ) <
        1e-12
    ) %>%
    select(
      channel,
      efficiency,
      efficiency_se,
      build_probability
    ) %>%
    arrange(
      desc(
        efficiency
      )
    ) %>%
    as.data.frame(),
  digits = 5,
  row.names = FALSE
)


cat(
  "\nSigned realized-revenue ranking at epsilon = 1, lambda = 1:\n"
)


print(
  revenue_results %>%
    filter(
      abs(
        epsilon -
          1
      ) <
        1e-12
    ) %>%
    select(
      channel,
      revenue,
      revenue_se,
      build_probability
    ) %>%
    arrange(
      desc(
        revenue
      )
    ) %>%
    as.data.frame(),
  digits = 5,
  row.names = FALSE
)


cat(
  "\nMaximum signal-grid overflow rates:\n"
)


overflow_summary <- tower_check %>%
  group_by(
    panel,
    channel
  ) %>%
  summarize(
    maximum_overflow =
      max(
        signal_outside_grid_rate
      ),
    
    .groups = "drop"
  )


print(
  as.data.frame(
    overflow_summary
  ),
  digits = 8,
  row.names = FALSE
)


cat(
  "\nMaximum absolute tower-property errors:\n"
)


tower_summary <- tower_check %>%
  group_by(
    panel,
    channel
  ) %>%
  summarize(
    maximum_xhat_error =
      max(
        abs(
          aggregate_xhat_error
        )
      ),
    
    maximum_Jhat_error =
      max(
        abs(
          aggregate_Jhat_error
        )
      ),
    
    .groups = "drop"
  )


print(
  as.data.frame(
    tower_summary
  ),
  digits = 6,
  row.names = FALSE
)


cat(
  paste0(
    "\nInterpretation:\n",
    "  * All three channels use the same (epsilon,delta)-LDP frontier.\n",
    "  * Epsilon and lambda are displayed on ordinary linear axes.\n",
    "  * Signed revenue is not truncated at zero.\n",
    "  * Negative revenue is possible under the knife-edge prior.\n",
    "  * Welfare efficiency is W/W_FB by definition.\n",
    "  * Shared type and uniform-noise draws are used across channels.\n",
    "  * Signal-grid overflow rates should be negligible.\n",
    "  * Tower-property errors should be small relative to scale.\n"
  )
)


# ================================================================
# 27. COMPLETION REPORT
# ================================================================

cat(
  "\nFigures written to:\n"
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp7_efficiency_vs_eps.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp7_revenue_vs_eps.pdf"
  ),
  "\n",
  sep = ""
)


cat(
  "  ",
  file.path(
    fig_dir,
    "fig_exp7_objective_vs_lambda.pdf"
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