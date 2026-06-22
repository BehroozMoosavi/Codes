test_that("ULA is reproducible under set.seed()", {
  g <- function(x) x

  set.seed(123)
  a <- ula(c(1, -1), g, 0.05, 200)
  set.seed(123)
  b <- ula(c(1, -1), g, 0.05, 200)

  expect_identical(a$samples, b$samples)
})

test_that("MALA is reproducible under set.seed()", {
  U <- function(x) 0.5 * sum(x^2)
  g <- function(x) x

  set.seed(456)
  a <- mala(c(1, -1), U, g, 0.3, 200)
  set.seed(456)
  b <- mala(c(1, -1), U, g, 0.3, 200)

  expect_identical(a$samples, b$samples)
  expect_identical(a$accepted, b$accepted)
})

test_that("Different seeds produce different chains", {
  g <- function(x) x

  set.seed(1)
  a <- ula(c(0, 0), g, 0.05, 200)
  set.seed(2)
  b <- ula(c(0, 0), g, 0.05, 200)

  expect_false(isTRUE(all.equal(a$samples, b$samples)))
})
