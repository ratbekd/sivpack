#' Computes Two-sample Anderson-Darling statistic
#'
#' @param x Numeric vector (first sample).
#' @param y Numeric vector (second sample).
#' @return Two-sample Anderson-Darling statistic.
#' @export
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100, mean = 0.5)
#' ad2_stat(x, y)

ad2_stat <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y)) stop("Inputs must be numeric vectors")
  if (length(x) < 2 || length(y) < 2) stop("Samples must have at least two observations")

  # Sample sizes
  n <- length(x)
  m <- length(y)

  # Pooled sample, sorted
  z <- sort(c(x, y))
  z <- z[-length(z)]  # Remove the largest value robustly

  # Pooled ECDF (Explicitly referencing `stats::ecdf`)
  ecdf_x <- stats::ecdf(x)
  ecdf_y <- stats::ecdf(y)
  W <- rank(z) / (n + m)

  # Compute Anderson-Darling statistic
  eps <- .Machine$double.eps  # Avoid division by zero
  stat <- (n * m / (n + m)^2) * sum((ecdf_x(z) - ecdf_y(z))^2 / pmax((1 - W) * W, eps))

  return(stat)
}

