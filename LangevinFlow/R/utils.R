# Internal helpers — not exported

# Validate a numeric initial state vector.
.check_init <- function(init_x) {
  if (!is.numeric(init_x) || length(init_x) < 1L) {
    stop("'init_x' must be a non-empty numeric vector.", call. = FALSE)
  }
  if (any(!is.finite(init_x))) {
    stop("'init_x' must contain only finite values.", call. = FALSE)
  }
  as.numeric(init_x)
}

# Validate a positive scalar.
.check_pos_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be a single positive finite number.", name),
         call. = FALSE)
  }
  as.numeric(x)
}

# Validate a positive integer.
.check_pos_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x <= 0 || x != round(x)) {
    stop(sprintf("'%s' must be a single positive integer.", name),
         call. = FALSE)
  }
  as.integer(x)
}

# Validate a function and (optionally) probe its output shape at init_x.
.check_fun <- function(f, name, init_x = NULL, expect_scalar = FALSE) {
  if (!is.function(f)) {
    stop(sprintf("'%s' must be a function.", name), call. = FALSE)
  }
  if (!is.null(init_x)) {
    out <- tryCatch(f(init_x), error = function(e) {
      stop(sprintf("Calling '%s(init_x)' failed: %s", name, conditionMessage(e)),
           call. = FALSE)
    })
    if (expect_scalar) {
      if (!is.numeric(out) || length(out) != 1L || !is.finite(out)) {
        stop(sprintf("'%s' must return a single finite numeric value.", name),
             call. = FALSE)
      }
    } else {
      if (!is.numeric(out) || length(out) != length(init_x)) {
        stop(sprintf(
          "'%s' must return a numeric vector of length %d (got length %d).",
          name, length(init_x), length(out)),
          call. = FALSE)
      }
      if (any(!is.finite(out))) {
        stop(sprintf("'%s' returned non-finite values at 'init_x'.", name),
             call. = FALSE)
      }
    }
  }
  f
}

# Construct the S3 chain object returned by ula() / mala().
.new_langevin_chain <- function(samples,
                                algorithm,
                                step_size,
                                beta,
                                n_iter,
                                burn_in,
                                acceptance_rate = NA_real_,
                                accepted = NULL,
                                elapsed = NA_real_) {
  if (burn_in > 0 && burn_in < nrow(samples)) {
    samples_post <- samples[-seq_len(burn_in), , drop = FALSE]
    if (!is.null(accepted)) accepted <- accepted[-seq_len(burn_in)]
  } else {
    samples_post <- samples
  }

  out <- list(
    samples         = samples_post,
    algorithm       = algorithm,
    step_size       = step_size,
    beta            = beta,
    n_iter          = n_iter,
    burn_in         = burn_in,
    dimension       = ncol(samples_post),
    acceptance_rate = acceptance_rate,
    accepted        = accepted,
    elapsed_secs    = elapsed
  )
  class(out) <- "langevin_chain"
  out
}
