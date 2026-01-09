\name{bmdd}
\alias{bmdd}
\title{Bimodal Dirichlet Distribution Estimation}
\description{
Estimates parameters of the bimodal Dirichlet distribution using a variational
EM algorithm. Automatically selects the optimal implementation: high-performance
C++ (NLopt) when possible, or R when covariates are needed.
}
\usage{
bmdd(W, type = c("count", "proportion"),
     Z = NULL, formula.Z = NULL, U = NULL, formula.U = NULL,
     Z.standardizing = TRUE, U.standardizing = TRUE,
     alp.eta = FALSE, alp.kap = FALSE,
     pi.xi = FALSE, pi.zeta = FALSE,
     para.alp.init = NULL, para.pi.init = NULL, gam.init = NULL,
     iterlim = 500, tol = 1e-6, trace = FALSE,
     method = c("auto", "nlopt", "R"),
     inner.loop = TRUE, inner.iterlim = 20, inner.tol = 1e-6,
     alp.iterlim = 100, alp.tol = 1e-6,
     alp.min = 1e-3, alp.max = 1e3, ...)
}
\arguments{
  \item{W}{Matrix of observed data (m taxa × n samples)}
  \item{type}{Data type: "count" or "proportion"}
  \item{Z}{Optional matrix of covariates for alpha (forces R implementation)}
  \item{formula.Z}{Optional formula for Z covariates}
  \item{U}{Optional matrix of covariates for pi (forces R implementation)}
  \item{formula.U}{Optional formula for U covariates}
  \item{Z.standardizing}{Standardize Z covariates (default TRUE)}
  \item{U.standardizing}{Standardize U covariates (default TRUE)}
  \item{alp.eta}{Model alpha0 as function of Z (default FALSE)}
  \item{alp.kap}{Model alpha1 as function of Z (default FALSE)}
  \item{pi.xi}{Model pi as function of U (default FALSE)}
  \item{pi.zeta}{Model pi variance as function of U (default FALSE)}
  \item{para.alp.init}{Initial alpha values}
  \item{para.pi.init}{Initial pi values}
  \item{gam.init}{Initial gamma values}
  \item{iterlim}{Maximum iterations (default 500)}
  \item{tol}{Convergence tolerance (default 1e-6)}
  \item{trace}{Print progress (default FALSE)}
  \item{method}{Force method: "auto", "nlopt", or "R" (default "auto")}
  \item{inner.loop}{Use inner loop for NLopt (default TRUE)}
  \item{inner.iterlim}{Inner loop max iterations (default 20)}
  \item{inner.tol}{Inner loop tolerance (default 1e-6)}
  \item{alp.iterlim}{Alpha optimization iterations (default 100)}
  \item{alp.tol}{Alpha tolerance (default 1e-6)}
  \item{alp.min}{Minimum alpha (default 1e-3)}
  \item{alp.max}{Maximum alpha (default 1e3)}
  \item{...}{Additional arguments (ignored)}
}
\value{
  A list containing:
  \item{gamma}{Estimated gamma parameters (bimodality indicators)}
  \item{pi}{Estimated pi parameters (mixing proportions)}
  \item{beta}{Estimated posterior Dirichlet parameters}
  \item{alpha0}{Estimated alpha0 parameters (mode 0)}
  \item{alpha1}{Estimated alpha1 parameters (mode 1)}
  \item{converge}{Logical indicating convergence}
  \item{iter}{Number of iterations}
  \item{method}{Method used: "nlopt" or "R"}
}
\details{
Automatically chooses best implementation:
\itemize{
  \item \strong{NLopt C++}: When no covariates. 50-90x faster using L-BFGS-B
    with analytical gradients. Scientifically equivalent to R (r > 0.999).
  \item \strong{R}: When covariates needed or explicitly requested.
    Supports full covariate modeling.
}

Requires NLopt library for high-performance mode:
\itemize{
  \item macOS: \code{brew install nlopt}
  \item Linux: \code{sudo apt-get install libnlopt-dev}
}

Refer to \url{https://github.com/zhouhj1994/BMDD} for detailed examples about zero 
imputation and posterior sample generation. 
}
\references{
Zhou, H., Chen, J., & Zhang, X. (2025). BMDD: A probabilistic framework for 
accurate imputation of zero-inflated microbiome sequencing data. PLoOS 
Computational Biology, 21(10), e1013124.
}
\examples{
\dontrun{
# Simulate data
m <- 100  # taxa
n <- 50   # samples
W <- matrix(rpois(m*n, 100), m, n)

# Auto-select method (uses NLopt for speed)
result <- bmdd(W, type = "count")

# Access results
head(result$beta)    # Posterior parameters
head(result$gamma)   # Bimodality indicators
result$method        # Shows which method was used

# Force specific method
result_nlopt <- bmdd(W, method = "nlopt")  # Force NLopt
result_r <- bmdd(W, method = "R")          # Force R

# With covariates (automatically uses R)
Z <- matrix(rnorm(m * 2), m, 2)
result_cov <- bmdd(W, Z = Z, alp.eta = TRUE)
}
}
