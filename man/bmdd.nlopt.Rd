\name{bmdd.nlopt}
\alias{bmdd.nlopt}
\title{NLopt C++ Implementation of BMDD}
\description{
High-performance implementation using NLopt L-BFGS-B optimizer in C++.
Provides 50-90x speedup over R implementation for cases without covariates.
This is a low-level function; most users should use \code{bmdd()} instead.
}
\usage{
bmdd.nlopt(W, type = c("count", "proportion"),
           para.alp.init = NULL, para.pi.init = NULL, gam.init = NULL,
           iterlim = 500, tol = 1e-6, trace = FALSE,
           inner.loop = TRUE, inner.iterlim = 20, inner.tol = 1e-6,
           alp.iterlim = 100, alp.tol = 1e-6,
           alp.min = 1e-3, alp.max = 1e3)
}
\arguments{
  \item{W}{Matrix of observed data (m taxa × n samples)}
  \item{type}{Data type: "count" or "proportion"}
  \item{para.alp.init}{Initial alpha values}
  \item{para.pi.init}{Initial pi values}
  \item{gam.init}{Initial gamma values}
  \item{iterlim}{Maximum iterations (default 500)}
  \item{tol}{Convergence tolerance (default 1e-6)}
  \item{trace}{Print progress (default FALSE)}
  \item{inner.loop}{Use inner loop optimization (default TRUE)}
  \item{inner.iterlim}{Inner loop max iterations (default 20)}
  \item{inner.tol}{Inner loop tolerance (default 1e-6)}
  \item{alp.iterlim}{Alpha optimization iterations (default 100)}
  \item{alp.tol}{Alpha tolerance (default 1e-6)}
  \item{alp.min}{Minimum alpha (default 1e-3)}
  \item{alp.max}{Maximum alpha (default 1e3)}
}
\value{
  A list containing estimated parameters (same as \code{bmdd})
}
\details{
This function provides direct access to the NLopt C++ implementation.
Most users should use \code{bmdd()} instead, which automatically
selects the optimal implementation.

Does not support covariates. For covariate modeling, use \code{bmdd()}
with method="R".
}

\references{
Zhou, H., Chen, J., & Zhang, X. (2025). BMDD: A probabilistic framework for 
accurate imputation of zero-inflated microbiome sequencing data. PLoOS 
Computational Biology, 21(10), e1013124.
}

\examples{
\dontrun{
# Simulate data
m <- 100
n <- 50
W <- matrix(rpois(m*n, 100), m, n)

# Direct NLopt usage (advanced)
result <- bmdd.nlopt(W, type = "count")

# Recommended: use bmdd() instead
result <- bmdd(W, type = "count")  # auto-selects NLopt
}
}
\seealso{
  \code{\link{bmdd}} for the main user interface
}
