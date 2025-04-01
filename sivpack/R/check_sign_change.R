#' Find the  Sign Change in a Series
#'
#' @param x A numeric vector.
#' @return The index of the first sign change, or NA if no sign change occurs.
#' @export
#' @examples
#' check_sign_change(c(-1, -2, -3, 3, 4))  # Returns 4
#' check_sign_change(c(5, 4, 3, 2, 1))  # Returns NA (no sign change)
#' check_sign_change(c(0, 0, 1, -1, 2))  # Returns 4

check_sign_change <- function(x) {
  # Ensure input is numeric
  if (!is.numeric(x)) stop("Input must be a numeric vector")

  # Handle cases where length is too short
  if (length(x) < 2) return(NA)

  # Adjust sign function to avoid ignoring zero transitions
  adjusted_sign <- sign(x) + (x == 0) * 1e-10

  # Detect first sign change
  sign_changes <- which(diff(adjusted_sign) != 0)

  # Return the first index where the sign changes
  return(ifelse(length(sign_changes) > 0, sign_changes[1] + 1, NA))
}
