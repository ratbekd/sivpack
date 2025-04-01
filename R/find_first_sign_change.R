#' Find the First Sign Change in a Series
#'
#' @param x A numeric vector.
#' @return The index of the first sign change, or NA if no sign change occurs.
#' @export
#' @examples
#' find_first_sign_change(c(-1, -2, -3, 3, 4))  # Returns 4
#' find_first_sign_change(c(5, 4, 3, 2, 1))  # Returns NA (no sign change)
#' #'
find_first_sign_change <- function(x) {
sign_changes <- which(diff(sign(x)) != 0)  # Find indices where sign changes
if (length(sign_changes) > 0) {
  return(sign_changes[1] + 1)  # Return the first occurrence (adjust for diff)
} else {
  return(NA)  # Return NA if no sign change
}
}
