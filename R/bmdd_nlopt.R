#' Production-ready BMDD using NLopt L-BFGS-B optimizer
#'
#' This version uses NLopt's L-BFGS-B optimizer in pure C++ for maximum speed
#' while maintaining the same convergence quality as nlminb.
#'
#' @inheritParams bmdd
#' @export

bmdd.nlopt <- function(W, type = c('count', 'proportion'),
                      para.alp.init = NULL, para.pi.init = NULL, gam.init = NULL,
                      iterlim = 500, tol = 1e-6, trace = FALSE,
                      inner.loop = TRUE, inner.iterlim = 20, inner.tol = 1e-6,
                      alp.iterlim = 100, alp.tol = 1e-6,
                      alp.min = 1e-3, alp.max = 1e3) {

  # Check if NLopt is available
  if (!test_nlopt_available()) {
    stop('NLopt library is not available. Please use bmdd() instead, or install NLopt and reinstall the package.\n',
         'Installation instructions:\n',
         '  - Ubuntu/Debian: sudo apt-get install libnlopt-dev\n',
         '  - macOS: brew install nlopt\n',
         '  - Fedora/RHEL: sudo dnf install NLopt-devel\n')
  }

  if(any(is.na(W))) {
    stop('The OTU table contains NAs! Please remove!\n')
  }

  W <- as.matrix(W)

  # Call NLopt C++ function
  result <- bmdd_nlopt(
    W = W,
    type = type,
    para_alp_init = para.alp.init,
    para_pi_init = para.pi.init,
    gam_init = gam.init,
    iterlim = iterlim,
    tol = tol,
    trace = trace,
    inner_loop = inner.loop,
    inner_iterlim = inner.iterlim,
    inner_tol = inner.tol,
    alp_iterlim = alp.iterlim,
    alp_tol = alp.tol,
    alp_min = alp.min,
    alp_max = alp.max
  )

  return(result)
}
