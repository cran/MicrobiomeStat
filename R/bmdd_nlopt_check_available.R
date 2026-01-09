#' Check if NLopt is available
#'
#' This function checks whether the package was compiled with NLopt support.
#' NLopt is an optional dependency that enables the high-performance bmdd.nlopt() function.
#'
#' @return Logical value: TRUE if NLopt is available, FALSE otherwise
#' @export
#'
#' @examples
#' # Check if NLopt is available
#' if (bmdd_nlopt_check_available()) {
#'   message("NLopt is available - you can use bmdd.nlopt()")
#' } else {
#'   message("NLopt is not available - please use bmdd() instead")
#' }
bmdd_nlopt_check_available <- function() {
  test_nlopt_available()
}
