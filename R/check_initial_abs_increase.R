#' Check if Absolute Values Initially Increase
#'
#' @param x A numeric vector.
#' @return 1 if absolute values start increasing, 0 otherwise.
#' @export
#' @examples
#' check_initial_abs_increase(c(-1, -2, -3, -2, -1))  # Returns 0
#' check_initial_abs_increase(c(0.1, 0.5, 1.0, 2.0, 3.0))  # Returns 1

check_initial_abs_increase <- function(x) {
  # Ensure input is numeric
  if (!is.numeric(x)) stop("Input must be a numeric vector")

  # Handle small inputs
  if (length(x) < 2) return(0)

  # Compute absolute values
  abs_x <- abs(x)

  # Extract the first 20 (or fewer) elements
  initial_segment <- abs_x[1:min(20, length(abs_x))]

  # Check if there's an increasing trend
  increasing <- all(diff(initial_segment) >= 0) && any(diff(initial_segment) > 0)

  return(as.integer(increasing))
}
