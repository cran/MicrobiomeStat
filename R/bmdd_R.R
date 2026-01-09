beta.gam.fun <- function(gam, alp, pi, W, type, m,
                         inner.loop, inner.iterlim, inner.tol, inner.trace) {
  alp0 <- alp$alp0
  alp1 <- alp$alp1
  alp.diff <- alp0 - alp1
  lg.alp.diff <- lgamma(alp0) - lgamma(alp1)
  alp.gam <- alp0 - gam * alp.diff
  csum.alp.gam <- colSums(alp.gam)
  if(type == 'count') {
    beta <- W + alp.gam
    csum.beta <- colSums(beta)
  } else {
    r <- alp.diff * log(W) - lg.alp.diff
  }

  if(inner.loop) {
    alp0 <- as.matrix(alp0)
    alp1 <- as.matrix(alp1)
    pi <- as.matrix(pi)
    alp.diff <- as.matrix(alp.diff)
    lg.alp.diff <- as.matrix(lg.alp.diff)

    iter <- 1
    while (iter <= inner.iterlim) {
      gam.tmp <- gam
      for(i in 1 : m) {
        tmp <- csum.alp.gam - alp0[i,] + gam[i,] * alp.diff[i,]
        if(type == 'count') {
          h <- lgamma(alp0[i,] + tmp) - lgamma(alp1[i,] + tmp) +
            alp.diff[i,] * (digamma(beta[i,]) - digamma(csum.beta)) - lg.alp.diff[i,]
        } else {
          h <- lgamma(alp0[i,] + tmp) - lgamma(alp1[i,] + tmp) + r[i,]
        }
        gam[i,] <- 1 / (1 + (1 / pi[i,] - 1) * exp(h))
        alp.gam.i <- alp0[i,] - gam[i,] * alp.diff[i,]
        csum.alp.gam <- tmp + alp.gam.i
        if(type == 'count') {
          csum.beta <- csum.beta - beta[i,]
          beta[i,] <- W[i,] + alp.gam.i
          csum.beta <- csum.beta + beta[i,]
        }
      }
      if(inner.trace) {
        print(paste0('inner loop: ', iter))
      }
      if(sum((gam - gam.tmp) ^ 2) / sum(gam.tmp ^ 2) < inner.tol) {
        break
      }
      iter <- iter + 1
    }
  } else {
    alp.gam <- t(-t(alp.gam) + csum.alp.gam)
    if(type == 'count') {
      h <- lgamma(alp0 + alp.gam) - lgamma(alp1 + alp.gam) +
        alp.diff * t(t(digamma(beta)) - digamma(csum.beta)) - lg.alp.diff
    } else {
      h <- lgamma(alp0 + alp.gam) - lgamma(alp1 + alp.gam) + r
    }
    gam <- 1 / (1 + (1 / pi - 1) * exp(h))
  }

  if(type == 'count') {
    res <- list(beta = beta, gam = gam)
  } else {
    res <- gam
  }
  return(res)
}

alp.fun <- function(para, alp.eta = FALSE, alp.kap = FALSE,
                    Z = NULL, U = NULL,
                    m = NULL, n = NULL, mp = NULL, nq = NULL) {
  if(!alp.eta & !alp.kap) {
    alp0 <- para[1:m]
    alp1 <- para[(m+1):(2*m)]
  } else {
    if(alp.eta) {
      eta0 <- matrix(para[1:mp],m)
      eta1 <- matrix(para[(mp+1):(2*mp)],m)
      tmp0 <- eta0 %*% t(Z)
      tmp1 <- eta1 %*% t(Z)
      if(alp.kap) {
        kap0 <- matrix(para[(2*mp+1):(2*mp+nq)],n)
        kap1 <- matrix(para[(2*mp+nq+1):(2*mp+2*nq)],n)
        tmp0 <- tmp0 + U %*% t(kap0)
        tmp1 <- tmp1 + U %*% t(kap1)
      }
    } else if(alp.kap) {
      kap0 <- matrix(para[1:nq],n)
      kap1 <- matrix(para[(nq+1):(2*nq)],n)
      tmp0 <- U %*% t(kap0)
      tmp1 <- U %*% t(kap1)
    }
    alp0 <- exp(tmp0)
    alp1 <- exp(tmp1)
  }
  return(list(alp0 = alp0, alp1 = alp1))
}

pi.fun <- function(para, pi.xi = FALSE, pi.zeta = FALSE,
                   Z = NULL, U = NULL,
                   m = NULL, n = NULL, mp = NULL, nq = NULL) {
  if(!pi.xi & !pi.zeta) {
    pi <- para
  } else {
    if(pi.xi) {
      xi <- matrix(para[1:mp],m)
      tmp <- xi %*% t(Z)
      if(pi.zeta) {
        zeta <- matrix(para[(mp+1):(mp+nq)],n)
        tmp <- tmp + U %*% t(zeta)
      }
    } else if(pi.zeta) {
      zeta <- matrix(para[1:nq],n)
      tmp <- U %*% t(zeta)
    }
    pi <- 1 / (1 + exp(-tmp))
  }
  return(pi)
}

#' Bimodal Dirichlet distribution (BMDD) for microbiome data
#'
#' The function implements a variational EM algorithm for estimating the parameters in bimodal Dirichlet distribution.
#' For the microbiome data, the unobserved true compositions are assumed to follow the bimodal Dirichlet distribution,
#' thus the function also estimates the posterior distribution of the true compositions.
#'
#' Model and Notations: Assume there are m taxa and n samples in the microbiome count data.
#' Given the sequencing depth and true compositions of a sample, which is assumed to be from the bimodal Dirichlet distribution,
#' the counts of taxa in the sample can be considered following multinormial distribution. Let (alpha_1,...,alpha_m) be the
#' parameters of a Dirichlet distribution. The bimodal Dirichlet distribution allows two possibilities for each alpha_i (i.e, the two
#' modals of the distribution of the corresponding composition), alpha0_i and alpha1_i, and the probability of being alpha1_i is referred to as pi_i.
#' Let pi = (pi_1,pi_2,...,pi_m), alpha0 = (alpha0_1,alpha0_2,...,alpha0_m) and alpha1 = (alpha1_1,alpha1_2,...,alpha1_m).
#' We use variational EM algorithm to estimate the parameters based on the mean-field approximation for the posterior.
#' Use delta_ij~Bernoulli(pi_i) to indicate which modal the taxon i in sample j is from. The mean-field approximation for the posterior of delta_ij
#' is Bernoulli(gamma_ij). The mean-filed approximation for the posterior of true compositions of sample j is Dirichlet(beta_1j,beta_2j,...,beta_mj).
#' Let gamma and beta be m by n matrices with (i,j)th entry being gamma_ij and beta_ij respectively.
#'
#' Incorporating sample or taxon covariates: We incorporate the covariates for alpha through the log linear model, and pi the logistic model.
#' Suppose the number of regression variables (after considering the dummy variables) from the sample covariates is p, and from taxon covariates is q.
#' Use eta0 and eta1 to denote the regression coefficients on sample covariates for alpha0 and alpha1,
#' and use kappa0 and kappa1 to denote the regression coefficients on taxon covariates for alpha0 and alpha1.
#' Use xi to denote the regression coefficients on sample covariates for pi,
#' and use zeta to denote the regression coefficients on taxon covariates for pi.
#' Thus, the lengths of eta0, eta1 and xi are mp and the lengths of kappa0, kappa1 and zeta are nq.
#' In this case, alpha0, alpha1 and pi are matrices of m by n.
#'
#' @param W data frame or matrix representing count table or table of samples from bimodal Dirichlet distribution.
#' Row: taxa; column: samples. NAs are not allowed.
#' @param type either 'count' or 'proportion'. If \code{W} is a count table of microbiome data, then \code{type = 'count'}; and
#' if \code{W} is a table of samples from bimodal Dirichlet distribution, then \code{type = 'proportion'}.
#' @param Z data frame of sample covariates. The rows of \code{Z} correspond to the columns of \code{W}. If NULL, no sample covariates are incorporated.
#' @param formula.Z formula of regression models on sample covariates. For example: \code{formula = '~z1*z2+z3'}, where z1, z2, and z3 are some of the
#' columns of \code{Z}. If \code{Z} is non-NULL, and \code{formula.Z} is NULL, then the \code{formula.Z} is set to be all columns of \code{Z} automatically.
#' @param U data frame of taxon covariates. The rows of \code{U} correspond to the rows of \code{W}. If NULL, no taxon covariates are incorporated.
#' @param formula.U formula of regression models on taxon covariates. If \code{U} is non-NULL, and \code{formula.U} is NULL,
#' then the \code{formula.U} is set to be all columns of \code{U} automatically.
#' @param Z.standardizing a logical value indicating whether to standardize the data. If TRUE, the numerical columns in \code{formula.Z} will be standardized.
#' @param U.standardizing a logical value indicating whether to standardize the data. If TRUE, the numerical columns in \code{formula.U} will be standardized.
#' @param alp.eta a logical value indicating whether to incorporate sample covariates for alpha.
#' @param alp.kap a logical value indicating whether to incorporate taxon covariates for alpha.
#' @param pi.xi a logical value indicating whether to incorporate sample covariates for pi.
#' @param pi.zeta a logical value indicating whether to incorporate taxon covariates for pi.
#' @param para.alp.init a vector of initial values for (alpha0, alpha1) or (eta0, eta1, kappa0, kappa1).
#' If no covariates are considered, then the length of \code{para.alp.init} is 2m; otherwise,
#' the length of \code{para.alp.init} is 2(mpk+nql), where k=1 if \code{alp.eta} is TRUE and 0 otherwise, and l=1
#' if \code{alp.kap} is TRUE and 0 otherwise. If not given, the function will do the initialization using random numbers.
#' @param para.pi.init a vector of initial values for pi or (xi, zeta).
#' If no covariates are considered, then the length of \code{para.pi.init} is m; otherwise,
#' the length of \code{para.pi.init} is mpk+nql, where k=1 if \code{pi.xi} is TRUE and 0 otherwise, and l=1
#' if \code{pi.zeta} is TRUE and 0 otherwise. If not given, the function will do the initialization.
#' @param gam.init initial values of gamma. If not given, the function will do the initialization.
#' @param iterlim maximum number of iterations of the variational EM algorithm.
#' @param tol a numerical value giving the tolerance in the relative change in the estimates of alpha0 and alpha1 below which
#' the algorithm is considered to be converged.
#' @param trace a logical value indicating whether to print the iteration process.
#' @param inner.loop a logical value indicating whether to conduct iteration across j=1,...,m inside each E-step.
#' @param inner.iterlim maximum number of iterations of the loop inside each E-step.
#' @param inner.tol a numerical value giving the tolerance in the relative change in gamma inside the E-step.
#' @param inner.trace a logical value indicating whether to print the iteration process inside the E-step.
#' @param alp.iterlim maximum number of iterations in the optimization function (nlminb) for alpha0 and alpha1 in M-step.
#' @param alp.tol tolerance in the optimization function for alpha0 and alpha1 in M-step.
#' @param alp.min minimum alpha0 and alpha1 value.
#' @param alp.max maximum alpha0 and alpha1 value.
#' @param para.alp.min minimum value of the regression coefficients for alpha0 and alpha1.
#' @param para.alp.min maximum value of the regression coefficients for alpha0 and alpha1.
#' @param pi.iterlim maximum number of iterations in the optimization function (nlminb) for pi in M-step.
#' @param pi.tol tolerance in the optimization function for pi in M-step.
#' @param para.pi.min minimum value of the regression coefficients for pi.
#' @param para.pi.min maximum value of the regression coefficients for pi.
#'
#' @return A list with the elements
#' \item{gamma}{estimate of gamma}
#' \item{beta}{estimate of beta. If \code{type='porportion'}, no this item is returned.}
#' \item{alpha}{a list with first element being the estimate of alpha0 and second the estimate of alpha1.}
#' \item{pi}{estimate of pi}
#' \item{para.alpha}{if no covariates are incorporated, then \code{para.alpha} is (alpha0, alpha1) itself. if there are covariates,
#' then \code{para.alpha} is the estimate of the regression coefficients. The length of \code{para.alpha} is the same as \code{para.alp.init}}
#' \item{para.pi}{if no covariates are incorporated, then \code{para.pi} is pi itself. if there are covariates,
#' then \code{para.pi} is the estimate of the regression coefficients. The length of \code{para.pi} the is same as \code{para.pi.init}}
#' @author Huijuan Zhou \email{huijuanzhou2019@gmail.com}
#' Lu Yang \email{yang.lu@mayo.edu}
#' Jun Chen \email{chen.jun2@mayo.edu}
#' Xianyang Zhang \email{zhangxiany@stat.tamu.edu}
#' @references Huijuan Zhou, Lu Yang, Jun Chen, and Xianyang Zhang. Bimodal Dirichlet Distributions And Its Application to Microbiome Data Analysis.
#' @examples
#'
#' #This is an example of using BMDD as an zero-imputation method for LinDA to identify colorectal cancer
#' #associated bacterial species. Data "phy" is a phyloseq-class experiment-level object.
#' #LinDA is a differential abundance analysis method for microbiome data.
#' #The package is available at https://CRAN.R-project.org/package=MicrobiomeStat
#' #and https://github.com/zhouhj1994/LinDA.
#' #Reference for the dataset: Yu et al. (2017). Metagenomic analysis of faecal microbiome as a tool
#' #towards targeted non-invasive biomarkers for colorectal cancer.
#' #Reference for LinDA: Zhou et al. (2022). LinDA: linear models for differential abundance analysis
#' #of microbiome compositional data.
#'
#' #install packages "phyloseq", "MicrobiomeStat"
#' data(phy)
#' otu_filter <- function(feature.dat, prev = 0.1, dep = 1000){
#'   idx <- apply(feature.dat, 1, function(x) sum(x > 0) > (ncol(feature.dat) * prev))
#'   idx2 <- colSums(feature.dat) > dep
#'   return(feature.dat[idx, idx2])
#' }
#' otu.tab <- as.data.frame(as.matrix(phyloseq::otu_table(phy)))
#' meta.dat <- as.data.frame(as.matrix(phyloseq::sample_data(phy)))
#' meta.dat$grp <- as.factor(meta.dat$grp)
#' feature.dat <- otu_filter(otu.tab)
#' meta.dat <- meta.dat[colnames(feature.dat), ]

#' bmdd.fit <- bmdd(W = feature.dat, type = 'count')
#' prop.bmdd <- t(t(bmdd.fit$beta) / colSums(bmdd.fit$beta))
#' bmdd.obj  <- MicrobiomeStat::linda(feature.dat = prop.bmdd, meta.dat = meta.dat,
#'                                    formula = '~grp', feature.dat.type = 'proportion')
#' bmdd.res <- bmdd.obj$output[[1]][,'padj',drop = F]

#' linda.obj  <- MicrobiomeStat::linda(feature.dat = feature.dat, meta.dat = meta.dat,
#'                                     formula = '~grp', feature.dat.type = 'count')
#' linda.res <- linda.obj$output[[1]][,'padj',drop = F]

#' @keywords internal
bmdd.R <- function(W, type = c('count', 'proportion'),
                 Z = NULL, formula.Z = NULL, U = NULL, formula.U = NULL,
                 Z.standardizing = TRUE, U.standardizing = TRUE,
                 alp.eta = FALSE, alp.kap = FALSE,
                 pi.xi = FALSE, pi.zeta = FALSE,
                 para.alp.init = NULL, para.pi.init = NULL, gam.init = NULL,
                 iterlim = 500, tol = 1e-6, trace  = FALSE,
                 inner.loop = TRUE, inner.iterlim = 20, inner.tol = 1e-6, inner.trace = FALSE,
                 alp.iterlim = 20, alp.tol = 1e-6, alp.min = 1e-3, alp.max = 1e3,
                 para.alp.min = -10, para.alp.max = 10,
                 pi.iterlim = 20, pi.tol = 1e-6,
                 para.pi.min = -100, para.pi.max = 100) {
  if(any(is.na(W))) {
    stop('The OTU table contains NAs! Please remove!\n')
  }
  W <- as.matrix(W)
  m <- nrow(W)
  n <- ncol(W)
  if(is.null(gam.init))  {
    gam.init <- matrix(runif(m * n), m)
  }

  if(!is.null(Z)) {
    if(is.null(formula.Z)) {
      formula.Z <- paste0('~',paste(names(Z),collapse = '+'))
    }
    allvars <- all.vars(as.formula(formula.Z))
    Z <- as.data.frame(Z[, allvars])
    names(Z) <- allvars
    if(Z.standardizing) {
      ind <- sapply(1 : ncol(Z), function(i) is.numeric(Z[, i]))
      Z[, ind] <- scale(Z[, ind])
    }
    Z <- model.matrix(as.formula(formula.Z), Z)
    mp <- m * ncol(Z)
  } else {
    mp <- 0
  }
  if(!is.null(U)){
    if(is.null(formula.U)) {
      formula.U <- paste0('~',paste(names(U),collapse = '+'))
    }
    allvars <- all.vars(as.formula(formula.U))
    U <- as.data.frame(U[, allvars])
    names(U) <- allvars
    if(U.standardizing) {
      ind <- sapply(1 : ncol(U), function(i) is.numeric(U[, i]))
      U[, ind] <- scale(U[, ind])
    }
    U <- model.matrix(as.formula(formula.U), U)
    nq <- n * ncol(U)
  } else {
    nq <- 0
  }

  if(!alp.eta & !alp.kap) {
    l.alp <- 2 * m
    alp.lower <- rep(alp.min, l.alp)
    alp.upper <- rep(alp.max, l.alp)
    if(is.null(para.alp.init)) {
      para.alp.init <- c(runif(m), runif(m, 1, 2))
    }
  } else {
    if(alp.eta) {
      l.alp <- 2 * mp
      if(alp.kap) l.alp <- l.alp + 2 * nq
    } else if(alp.kap) {
      l.alp <- 2 * nq
    }
    alp.lower <- rep(para.alp.min, l.alp)
    alp.upper <- rep(para.alp.max, l.alp)
    if(is.null(para.alp.init)) {
      para.alp.init <- c(runif(l.alp / 2, -1, 0), runif(l.alp / 2))
    }
  }
  if(!pi.xi & !pi.zeta) {
    l.pi <- m
    if(is.null(para.pi.init)) {
      para.pi.init <- rowMeans(gam.init)
    }
  } else {
    if(pi.xi) {
      l.pi <- mp
      if(pi.zeta) l.pi <- l.pi + nq
    } else if(pi.zeta) {
      l.pi <- nq
    }
    pi.lower <- rep(para.pi.min, l.pi)
    pi.upper <- rep(para.pi.max, l.pi)
    if(is.null(para.pi.init)) {
      para.pi.init <- runif(l.pi, -1, 1)
    }
  }

  alp.obj.fun <- function(para) {
    alp <- alp.fun(para, alp.eta, alp.kap, Z, U, m, n, mp, nq)
    alp0 <- alp$alp0
    alp1 <- alp$alp1
    alp.gam <- alp0 - gam * (alp0 - alp1)
    lg.alp.gam <- (1 - gam) * lgamma(alp0) + gam * lgamma(alp1)
    res <- sum(lgamma(colSums(alp.gam))) + sum(alp.gam * A - lg.alp.gam)
    return(-res)
  }
  
  # We have a closed form solution, right?
  pi.obj.fun <- function(para) {
    pi <- pi.fun(para, pi.xi, pi.zeta, Z, U, m, n, mp, nq)
    res <- sum(gam * log(pi) + (1 - gam) * log(1 - pi))
    return(-res)
  }

  if(type == 'proportion') {
    if(any(W == 0)) {
      W <- t(apply(W, 1, function (x) {
        x[x == 0] <- 0.5 * min(x[x != 0])
        return(x)
      }))
    }
    A <- log(W)
  }

  para.alp <- para.alp.init
  para.pi <- para.pi.init
  gam <- gam.init
  alp <- alp.fun(para.alp, alp.eta, alp.kap, Z, U, m, n, mp, nq)
  pi <- pi.fun(para.pi, pi.xi, pi.zeta, Z, U, m, n, mp, nq)
  iter <- 1

  while (iter <= iterlim) {
    para.alp.tmp <- para.alp
    if(type == 'count') {
      res <- beta.gam.fun(gam, alp, pi, W, type, m,
                          inner.loop, inner.iterlim, inner.tol, inner.trace)
      gam <- res$gam
      beta <- res$beta
      A <- t(t(digamma(beta)) - digamma(colSums(beta)))
    } else {
      gam <- beta.gam.fun(gam, alp, pi, W, type, m,
                          inner.loop, inner.iterlim, inner.tol, inner.trace)
    }
    para.alp <- nlminb(para.alp, alp.obj.fun, lower = alp.lower, upper = alp.upper,
                  control = list(iter.max = alp.iterlim, rel.tol = alp.tol))$par
    if(!pi.xi & !pi.zeta) {
      para.pi <- rowMeans(gam)
    } else {
      para.pi <- nlminb(para.pi, pi.obj.fun, lower = pi.lower, upper = pi.upper,
                        control = list(iter.max = pi.iterlim, rel.tol = pi.tol))$par
    }
    alp <- alp.fun(para.alp, alp.eta, alp.kap, Z, U, m, n, mp, nq)
    pi <- pi.fun(para.pi, pi.xi, pi.zeta, Z, U, m, n, mp, nq)

    if(trace) {
      cat('Iteration', iter, '\n')
    }
    if(sum((para.alp - para.alp.tmp) ^ 2) / sum(para.alp.tmp ^ 2) < tol) {
      break
    }
    iter <- iter + 1
  }

  if(type == 'count') {
    res <- list(gamma = gam, beta = beta, alpha = alp, pi = pi,
                para.alpha = para.alp, para.pi = para.pi)
  } else {
    res <- list(gamma = gam, alpha = alp, pi = pi,
                para.alpha = para.alp, para.pi = para.pi)
  }
  return(res)
}
