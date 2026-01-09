#========================
# Integrated Algorithm
#========================
# View.type: 1. "C", Compositional data;
# View.type: 2. "O", Omics Data
# In this version, we consider multi-variate view

# Globalize lambda.eta.seq
utils::globalVariables(c( "lambda.eta.seq" ))

#' Compositional Sparse Canonical Correlation Analysis
#'
#' A compositional sparse canonical correlation analysis (csCCA) framework for integrating microbiome data with other high-dimensional omics data.
#' @param Y a n*(K*p) matrix representing the observations.
#' @param View.ind a (K*p) integer vector indicating the classes of features. The features with the same View.ind is in the same class.
#' @param lambda.seq a K vector consisting of hyper-parameters.
#' @param a.old Optional initial value for the coefficient vector \code{a.new}.
#' @param View.type a K vector encoding the structure type of each feature class. There are two choices: "O" (Omics Data),"C" (Compositional Data).
#' @param eps.stop a numerical value controlling the convergence.
#' @param max.step an integer controlling the maximum step for interaction.
#' @param eps a numerical value controlling the convergence.
#' @param T.step an integer controlling the maximum step for interaction.
#'
#' @return \code{a.new} the estimated coefficient vector.
#' @export
#'
#' @examples
#' \donttest{
#' library(dplyr)
#'
#' n <- 100
#' p <- q <- 50
#' sigma.nu <- 5
#' sigma.eps <- 1
#' omega_X <- 0.85*c(rep(1/10,9),-9/10,rep(0,p-10))
#' omega_Y <- 0.85*c(seq(0.08,0.12,length = 10),rep(0,q-10))
#' Data1 <- DGP_OC(seed=10,n,p,q,sigma.nu,sigma.eps,omega_X,omega_Y)
#'
#' library(mlrMBO)
#' Res.sCCA.CV <- cscca.CV(Y=Data1$Y,View.ind=Data1$View.ind,
#'                           View.type=c("O","O"),
#'                           show.info = TRUE)
#'
#'
#' Res.CsCCA.CV <- cscca.CV(Y=Data1$Y,View.ind=Data1$View.ind,
#'                                    View.type=c("O","C"),
#'                                    show.info = TRUE)
#'
#' Res.sCCA <- cscca(Y=Data1$Y,View.ind=Data1$View.ind,
#'                      lambda.seq=Res.sCCA.CV$lam.opt.trgt,
#'                      View.type=c("O","O"))
#' Res.CsCCA <- cscca(Y=Data1$Y,View.ind=Data1$View.ind,
#'                      lambda.seq=Res.CsCCA.CV$lam.opt.trgt,
#'                      View.type=c("O","C"))
#' }
#' @references
#'
#' 1. Deng, L., Tang, Y., Zhang, X., et al. (2024). Structure-adaptive canonical correlation analysis for microbiome multi-omics data. Frontiers in Genetics, 15, 1489694.
#'
#' 2. Chen, J., Bushman, F. D., Lewis, J. D., et al. (2013). Structure-constrained sparse canonical correlation analysis with an application to microbiome data analysis. Biostatistics, 14(2), 244–258.
cscca <- function(Y,
                  View.ind,
                  lambda.seq,
                  a.old = NULL,
                  View.type=NULL,
                  eps.stop=0.0001,max.step=30,
                  eps=0.0001,T.step=1){
  n <- nrow(Y)
  p <- ncol(Y)
  K.View <- length(unique(View.ind))

  #===== Initiate parameter ========
  if(is.null(View.type)){
    View.type = rep("O",K.View)
  }
  #####################################
  #== Initiation a
  #Y.cen <- Y - matrix(colMeans(Y),ncol=dim(Y)[2],nrow=dim(Y)[1],byrow=T)
  Y.cen <- scale(Y)
  cov.Y <- t(Y.cen)%*%Y.cen/n
  if(is.null(a.old)){
    a.old <- rep(0,p)
    for(k in 1:K.View){
      ind.k <- which(View.ind==k)

      a.old[ind.k] <- as.vector(svd(cov.Y[ind.k,-ind.k],nu=1,nv=0)$u) #first left singular vector with unit length
    }
  }
  a.new <- a.old
  zeta.new <- rep(1,p)
  error <- 1
  step <- 0
  Random.a <- rep(TRUE,K.View)

  while(error>eps.stop & step<=max.step){
    for(k in 1:K.View){
      ## Interactively update from each view
      ind.k <- which(View.ind==k)

      lambda.k <- lambda.seq[k]*rep(1,length(ind.k))
      zeta.new.k <- zeta.new[ind.k]
      c_a <- cov.Y[ind.k,-ind.k,drop=F] %*% a.new[-ind.k]

      if(View.type[k]=="O"){
        result_a <- Soft_threshold_cpp(c_a,lambda.k)
        a.new.k <- result_a$b
        a.new.k <- as.numeric(a.new.k)
      }

      if(View.type[k]=="C"){
        result_a <- ADMM_ab_cpp(c_a,lambda.k)
        a.new.k <- result_a$b
        a.new.k <- as.numeric(a.new.k)
      }

      if (sum(abs(a.new.k))==0){
        Random.a[k] <- T
        a.new.k=(2*rbinom(length(ind.k),1,0.5)-1)
        a.new.k=a.new.k/(sum(a.new.k^2))^0.5
        zeta.new[ind.k] <- rep(1,length(ind.k))
      }else{
        Random.a[k] <- F
        a.new.k=a.new.k/(sum(a.new.k^2))^0.5

      }
      a.new[ind.k] <- a.new.k
      zeta.new[ind.k] <- zeta.new.k
    }

    #calculate the error
    error <- sum(abs(a.new-a.old))
    #judge the convergence
    step <- step + 1
    a.old <- a.new
  }

  for(k in 1:K.View){
    ## To avoid random result
    ind.k <- which(View.ind==k)
    if(Random.a[k]){a.new[ind.k] <- 0}
  }

  results <- list(a.new=a.new)

  return(results)
}


cscca.CV.sub <- function(Y,
                     sam.ind=NULL,
                     View.ind,
                     lambda.seq,
                     View.type=NULL,
                     eps.stop=0.0001,max.step=30,eps=0.0001,T.step=1,
                     n_fold=5,
                     Criterion="cov",
                     is.refit=F){
  n <- nrow(Y)
  if(is.null(sam.ind)){
    sam.ind <- sample(n)
  }

  K.View <- length(unique(View.ind))
  fold_num <- floor(n/n_fold)
  CCA.trgt.seq <- numeric(n_fold)
  CCA.trgt.train.seq <- numeric(n_fold)
  CCA.trgt.cov.seq <- numeric(n_fold)
  CCA.trgt.train.cov.seq <- numeric(n_fold)

  for(i_fold in 1: n_fold){
    # To determine the hyperparameter
    train.ind <- sam.ind[-(((i_fold-1)*fold_num+1):(i_fold*fold_num))]
    test.ind <- sam.ind[((i_fold-1)*fold_num+1):(i_fold*fold_num)]

    Res <- cscca(Y=Y[train.ind,,drop=F],
                 View.ind=View.ind,
                 lambda.seq=lambda.seq,
                 View.type = View.type,
                 a.old= NULL,
                 eps.stop = eps.stop,
                 max.step = max.step,
                 eps = eps,
                 T.step = T.step)

    a.hat <- Res$a.new
    a.hat.p <- (abs(a.hat)>0.0001)*1

    a.hat[abs(a.hat)<=0.0001] <- 0

    ind.sel <- which(abs(a.hat)>0.0001)
    if((length(ind.sel)!=0 & is.refit == T)){

      Res.refit <- cscca(Y=Y[train.ind,ind.sel,drop=F],
                         View.ind=View.ind[ind.sel],
                         lambda.seq=rep(0,K.View),
                         View.type = View.type,
                         a.old = Res$a.new[ind.sel],
                         eps.stop = eps.stop,
                         max.step = max.step,
                         eps = eps,
                         T.step = T.step)
      a.hat.refit <- a.hat
      a.hat.refit[ind.sel] <- Res.refit$a.new
      a.hat <- a.hat.refit
    }

    # Calculate the value for maximization
    Y.test.cen <- Y[test.ind,,drop=F]-matrix(colMeans(Y[test.ind,,drop=F]),ncol=dim(Y)[2],nrow=length(test.ind),byrow=T)

    cov.Y.test <- t(Y.test.cen)%*%Y.test.cen/length(test.ind)

    Y.train.cen <- Y[train.ind,,drop=F]-matrix(colMeans(Y[train.ind,,drop=F]),ncol=dim(Y)[2],nrow=length(train.ind),byrow=T)

    cov.Y.train <- t(Y.train.cen)%*%Y.train.cen/length(train.ind)
    for(k in 1:(K.View-1)){
      ind.k <- which(View.ind==k)
      for(l in (k+1):(K.View)){
        ind.l <- which(View.ind==l)

        cor.tmp <- a.hat[ind.k] %*% cov.Y.test[ind.k,ind.l] %*% a.hat[ind.l]/
          sqrt((a.hat[ind.k] %*% cov.Y.test[ind.k,ind.k] %*% a.hat[ind.k])*
                 (a.hat[ind.l] %*% cov.Y.test[ind.l,ind.l] %*% a.hat[ind.l]))
        cor.tmp <- cor(Y.test.cen[,ind.k] %*% a.hat[ind.k],
                       Y.test.cen[,ind.l] %*% a.hat[ind.l])[1]

        cor.tmp <- ifelse(is.na(cor.tmp),0,cor.tmp)
        cov.tmp <-   a.hat[ind.k] %*%cov.Y.test[ind.k,ind.l] %*% a.hat[ind.l]
        cov.tmp <- cov(Y.test.cen[,ind.k] %*% a.hat[ind.k],
                       Y.test.cen[,ind.l] %*% a.hat[ind.l])[1]

        cor.train.tmp <- a.hat[ind.k] %*% cov.Y.train[ind.k,ind.l] %*% a.hat[ind.l]/
          sqrt((a.hat[ind.k] %*% cov.Y.train[ind.k,ind.k] %*% a.hat[ind.k])*
                 (a.hat[ind.l] %*% cov.Y.train[ind.l,ind.l] %*% a.hat[ind.l]))
        cor.train.tmp <- cor(Y.train.cen[,ind.k] %*% a.hat[ind.k],
                             Y.train.cen[,ind.l] %*% a.hat[ind.l])[1]

        cor.train.tmp <- ifelse(is.na(cor.train.tmp),0,cor.train.tmp)
        cov.train.tmp <- cov(Y.train.cen[,ind.k] %*% a.hat[ind.k],
                             Y.train.cen[,ind.l] %*% a.hat[ind.l])[1]

        # Correlation as the criterion
        CCA.trgt.seq[i_fold] <- CCA.trgt.seq[i_fold] + cor.tmp
        CCA.trgt.train.seq[i_fold] <- CCA.trgt.train.seq[i_fold] + cor.train.tmp
        # We choose CV

        CCA.trgt.cov.seq[i_fold] <- CCA.trgt.cov.seq[i_fold] + cov.tmp
        CCA.trgt.train.cov.seq[i_fold] <- CCA.trgt.train.cov.seq[i_fold] + cov.train.tmp
      }
    }

  }
  if(Criterion=="cor"){
    CCA.trgt <- mean(CCA.trgt.seq)
  }else if(Criterion=="cov"){
    CCA.trgt <- mean(CCA.trgt.cov.seq)
  }else{
    stop("Criterion for cv is not defined, please choose \"cor\" or \"cov\"")
  }

  return(list(CCA.trgt = CCA.trgt,
              cor.CV= mean(CCA.trgt.seq),
              cov.CV= mean(CCA.trgt.cov.seq),
              cor.seq.CV = CCA.trgt.seq,
              cov.seq.CV = CCA.trgt.cov.seq
  ))
}


cscca.CV.Opt.assign <- function(hp.lower.seq,hp.upper.seq,
                                Y,
                                sam.ind,
                                View.ind,
                                lambda.eta.seq,
                                eta.lower.seq,eta.upper.seq,
                                View.type,
                                eps.stop,max.step,eps,T.step,
                                Criterion,
                                n_fold,is.refit){
  K.View <- length(hp.lower.seq)

  cscca.CV.Opt <- makeSingleObjectiveFunction(
    name = "cscca.CV.Opt",
    fn = function(Y=Y,
                  sam.ind = sam.ind,
                  View.ind=View.ind,
                  lambda.eta.seq=lambda.eta.seq,
                  View.type=View.type,
                  eps.stop=eps.stop,max.step=max.step,eps=eps,T.step=T.step,
                  n_fold=n_fold,
                  Criterion=Criterion,
                  is.refit=is.refit) {

      K.View <- length(unique(View.ind))
      lambda.seq <- lambda.eta.seq[1:K.View]
      Res <- cscca.CV.sub(Y=Y,
                      sam.ind = sam.ind,
                      View.ind=View.ind,
                      lambda.seq=lambda.seq,
                      View.type=View.type,
                      eps.stop=eps.stop,max.step=max.step,eps=eps,T.step=T.step,
                      n_fold=n_fold,Criterion=Criterion,is.refit=is.refit
      )
      res <- Res$CCA.trgt

      extras <- list(cor.CV=Res$cor.CV,cov.CV=Res$cov.CV)

      res <- BBmisc::setAttribute(res, "extras", extras)

      return(res)
    },
    par.set = makeParamSet(
      makeNumericVectorParam("lambda.seq", len = K.View,
                             lower = hp.lower.seq,
                             upper = hp.upper.seq)
    ),
    #has.simple.signature = FALSE,
    minimize = FALSE
  )
  return(cscca.CV.Opt)
}


#' Compositional Sparse Canonical Correlation Analysis (Cross Valication Version)
#'
#' The cross validation version of a compositional sparse canonical correlation analysis (sCCA) framework for integrating microbiome data with other high-dimensional omics data.
#' @param Y a n*(K*p) matrix representing the observations.
#' @param View.ind a (K*p) integer vector indicating the classes of features. The features with the same View.ind is in the same class.
#' @param View.type a K vector encoding the structure type of each feature class. There are two choices: "O" (Omics Data),"C" (Compositional Data).
#' @param eps.stop a numerical value controlling the convergence.
#' @param max.step an integer controlling the maximum step for interaction.
#' @param eps a numerical value controlling the convergence.
#' @param T.step an integer controlling the maximum step for interaction.
#' @param n_fold an integer representing the number of folds for cross validation.
#' @param seed.sam.ind a vector of the seeds for sampling.
#' @param show.info a bool suggesting whether to show information through the hyperparameter optimization.
#' @param Criterion a character indicating the criterion we choose for cross validation.
#' @param hp.lower a numerical value or K vector specifying the lower bound of the hyper-parameter.
#' @param hp.upper a numerical value or K vector specifying the upper bound of the hyper-parameter.
#' @param hp.eta.lower a numerical value or K vector specifying the lower bound of the hyper-parameter for eta.
#' @param hp.eta.upper a numerical value or K vector specifying the upper bound of the hyper-parameter for eta.
#' @param eta.warm.stat.mat a matrix providing statistics for warm start of eta.
#' @param opt_n_design an integer controlling the number of design points in the hyperparameter optimization.
#' @param opt_n_iter an integer controlling the number of iterations in the hyperparameter optimization.
#' @param des.init an initial design for hyperparameter optimization.
#' @param is.refit a bool suggesting whether to refit the model using the optimal hyper-parameters.
#' @param is.refix.eta a bool suggesting whether eta is fixed during refitting.
#' @param opt_n_design.eta_warm an integer controlling the number of design points for eta warm-start optimization.
#' @param opt_n_iter.eta_warm an integer controlling the number of iterations for eta warm-start optimization.
#' @param is.opt.hyper a bool suggesting whether to optimize the hyper-parameters.
#' @param hyper_n_grid an integer controlling the grid size for hyperparameter search.
#' @param ... additional arguments passed to the internal optimization procedures.
#' @return A list containing the following elements: (1) \code{a.hat.opt.trgt}: The coefficient vector estimated with the optimal hyper-parameter vector; (2) \code{lam.opt.trgt}: The optimal hyper-parameter vector.
#' @export
#'
#' @references
#' 1. Deng, L., Tang, Y., Zhang, X., et al. (2024). Structure-adaptive canonical correlation analysis for microbiome multi-omics data. Frontiers in Genetics, 15, 1489694.
#'
#' 2. Chen, J., Bushman, F. D., Lewis, J. D., et al. (2013). Structure-constrained sparse canonical correlation analysis with an application to microbiome data analysis. Biostatistics, 14(2), 244–258.
#' @examples
#' \donttest{
#' library(dplyr)
#'
#' n <- 100
#' p <- q <- 50
#' sigma.nu <- 5
#' sigma.eps <- 1
#' omega_X <- 0.85*c(rep(1/10,9),-9/10,rep(0,p-10))
#' omega_Y <- 0.85*c(seq(0.08,0.12,length = 10),rep(0,q-10))
#' Data1 <- DGP_OC(seed=10,n,p,q,sigma.nu,sigma.eps,omega_X,omega_Y)
#'
#' library(mlrMBO)
#' Res.sCCA.CV <- cscca.CV(Y=Data1$Y,View.ind=Data1$View.ind,
#'                           View.type=c("O","O"),
#'                           show.info = TRUE)
#'
#'
#' Res.CsCCA.CV <- cscca.CV(Y=Data1$Y,View.ind=Data1$View.ind,
#'                                    View.type=c("O","C"),
#'                                    show.info = TRUE)
#'
#' Res.sCCA <- cscca(Y=Data1$Y,View.ind=Data1$View.ind,
#'                      lambda.seq=Res.sCCA.CV$lam.opt.trgt,
#'                      View.type=c("O","O"))
#' Res.CsCCA <- cscca(Y=Data1$Y,View.ind=Data1$View.ind,
#'                      lambda.seq=Res.CsCCA.CV$lam.opt.trgt,
#'                      View.type=c("O","C"))
#' print(Res.sCCA.CV$Cri.opt.trgt)
#' print(Res.CsCCA.CV$Cri.opt.trgt)
#' }
cscca.CV <- function(Y,
                        View.ind,
                        View.type=NULL,
                        eps.stop=0.0001,max.step=30,eps=0.0001,T.step=10,
                        n_fold=5,
                        seed.sam.ind =NULL,
                        show.info=FALSE,
                        hp.lower=NULL,
                        hp.upper= NULL,
                        hp.eta.lower=NULL,
                        hp.eta.upper= NULL,
                        eta.warm.stat.mat = NULL,
                        opt_n_design=30,
                        opt_n_iter=20,
                        Criterion="cov",
                        des.init = NULL,
                        is.refit=F,
                        is.refix.eta=T,
                        opt_n_design.eta_warm=30,
                        opt_n_iter.eta_warm=20,
                        is.opt.hyper = TRUE,
                        hyper_n_grid = 20,
                        ...){
  # lambda.mat should be a matrix with num_lambda*num_view
  n <- nrow(Y)
  p <- ncol(Y)
  K.View <- length(unique(View.ind))
  lam.row <- 20 # Can be delet later
  fold_num <- floor(n/n_fold)

  # Create function to optimize hyper-parameters: lambda
  if(is.null(hp.lower)|is.null(hp.upper)){
    Hyper.D <- Hyper_Para.determine(Y=Y,View.ind=View.ind,
                                    View.type=View.type,
                                    eps.stop = eps.stop,max.step = max.step,
                                    eps = eps,T.step = T.step)
    hp.lower.seq=Hyper.D$lambda.lower.seq
    hp.upper.seq=Hyper.D$lambda.upper.seq
  }else{
    if(length(hp.lower)==1){
      hp.lower.seq=rep(hp.lower,K.View)
    }else{
      hp.lower.seq=hp.lower
    }
    if(length(hp.upper)==1){
      hp.upper.seq=rep(hp.upper,K.View)
    }else{
      hp.upper.seq=hp.upper
    }
  }

  ## Optimize hyper-parameters: eta
  if(is.null(hp.eta.lower)|is.null(hp.eta.upper)){
    eta.lower.seq <- rep(0,K.View)
    eta.upper.seq <- ifelse(sapply(View.type,
                                   function(view.type){
                                     view.type %in% c("O","C")}),0,0.5)
    hp.eta.lower <- unique(eta.lower.seq)
    hp.eta.upper <- unique(eta.upper.seq)
  }else{
    if(length(hp.lower)==1){
      eta.lower.seq=rep(hp.eta.lower,K.View)
    }else{
      eta.lower.seq=hp.eta.lower
    }
    if(length(hp.upper)==1){
      eta.upper.seq=rep(hp.eta.upper,K.View)
    }else{
      eta.upper.seq=hp.eta.upper
    }
  }

  if(is.null(seed.sam.ind)){
    seed.sam.ind <- sample(10000001:20000000,100)
  }


  Cri.opt.trgt <- matrix(NA,nrow=K.View,ncol=4)

  if(!is.opt.hyper){
    run = NULL
    ## Determine the hyperparameters to be searched
    for(k in 1:K.View){
      assign(paste0("hp",k),seq(hp.lower.seq[k],hp.upper.seq[k],length.out = hyper_n_grid))
    }
    eta0 = unique(seq(hp.eta.lower,hp.eta.upper,length.out=hyper_n_grid))
    hp.all = eval(parse(text =paste0("expand.grid(",paste0("hp",1:K.View,collapse = ","),",eta0)")))
    hp.all = cbind(hp.all,Var4=matrix(hp.all[,K.View+1],ncol=K.View-1,byrow = F))
    CV.seq = numeric(dim(hp.all)[1])
    ## Now, we use CV to choose the optimal hyper parameters
    for(ii.check in 1:dim(hp.all)[1]){
      lambda.seq.tmp <- unlist(hp.all[ii.check,1:K.View])

      set.seed(seed.sam.ind[1])
      sam.ind <- sample(1:n)

      CV.Res <- cscca.CV.sub(Y=Y,
                         sam.ind = sam.ind,
                         View.ind=View.ind,
                         lambda.seq=lambda.seq.tmp,
                         View.type=View.type,
                         eps.stop=eps.stop,max.step=max.step,eps=eps,T.step=T.step,
                         n_fold=n_fold,Criterion=Criterion,is.refit=is.refit)
      CV.seq[ii.check] <- unlist(CV.Res$CCA.trgt)

    }
    CV.mat = cbind(hp.all,CV.seq)
    ii.opt <-  which.max(CV.seq)
    lam.opt <- unlist(hp.all[ii.opt,1:K.View])
  }else{
    CV.mat <- NULL
    convrg.flag <- F

    ## A possible reason for too many positive coeffcients is the setting for sample splitting, so we need to resample.
    iter.sam <- 1
    max.sam <- 1


    while(!convrg.flag ){
      set.seed(seed.sam.ind[iter.sam])
      sam.ind <- sample(1:n)

      surr.km = makeLearner("regr.km", predict.type = "se",covtype = "matern3_2",jitter=TRUE,
                            control = list(trace = FALSE))
      surr.km = makeLearner("regr.randomForest", predict.type="se")

      control = makeMBOControl()

      control = setMBOControlInfill(control, crit = makeMBOInfillCritCB())




      cscca.CV.Opt <- cscca.CV.Opt.assign(hp.lower.seq=hp.lower.seq,
                                          hp.upper.seq=hp.upper.seq,
                                          eta.lower.seq=eta.lower.seq,
                                          eta.upper.seq=eta.upper.seq,
                                          Y=Y,
                                          sam.ind = sam.ind,
                                          View.ind=View.ind,
                                          lambda.eta.seq=lambda.eta.seq,
                                          View.type=View.type,
                                          eps.stop=eps.stop,max.step=max.step,eps=eps,T.step=T.step,
                                          n_fold=n_fold,Criterion=Criterion,is.refit=is.refit)

      # Generate design
      des = generateDesign(n = opt_n_design, par.set = getParamSet(cscca.CV.Opt), fun = randomLHS)

      # We can provide hyperparameters
      if(!is.null(des.init)){
        colnames(des.init) <- colnames(des)
        des <- rbind(des,des.init)
      }

      cscca.CV.Opt <- cscca.CV.Opt.assign(hp.lower.seq=hp.lower.seq/2,
                                          hp.upper.seq=hp.upper.seq*2,
                                          eta.lower.seq=eta.lower.seq,
                                          eta.upper.seq=eta.upper.seq,
                                          Y=Y,
                                          sam.ind = sam.ind,
                                          View.ind=View.ind,
                                          lambda.eta.seq=lambda.eta.seq,
                                          View.type=View.type,
                                          eps.stop=eps.stop,max.step=max.step,eps=eps,T.step=T.step,
                                          n_fold=n_fold,Criterion=Criterion,is.refit=is.refit)

      control = setMBOControlTermination(control, iters = opt_n_iter)
      lam_eta.opt <- NULL
      run <- NULL
      for(opt.ind in 1){
        lam_eta.opt <- tryCatch({
          if(!is.null(run)){
            des <- run$opt.path$env$path
          }

          run = mbo(cscca.CV.Opt,design = des, learner = surr.km,
                    control = control, show.info = show.info,
                    more.args=list(Y=Y,
                                   sam.ind = sam.ind,
                                   View.ind=View.ind,
                                   View.type=View.type,
                                   eps.stop=eps.stop,max.step=max.step,eps=eps,T.step=T.step,
                                   n_fold=n_fold,Criterion=Criterion,is.refit=is.refit))
          lam.opt <- (run$x)$lambda.seq
          list(lam.opt=lam.opt)
        }, error = function(e){
          print(e)
          lam_eta.opt
        })
      }

      if(is.null(lam_eta.opt)){
        max.sam <- max.sam+1
        lam.opt <- rep(0,K.View)
      }else{
        lam.opt <- lam_eta.opt$lam.opt
      }

      # To avoid too many positive coefficients
      convrg.flag=T
      iter.sam <- iter.sam+1
      if(iter.sam>max.sam){
        break
      }
    }

  }
  # Obtain the result for the optimal lambda
  # We just maintain the results for CCA.trgt
  Res <- cscca(Y=Y,View.ind=View.ind,
               lambda.seq=lam.opt,
               View.type = View.type,
               eps.stop = eps.stop,
               max.step = max.step,
               eps = eps,
               T.step = T.step)

  ind.sel <- which(abs(Res$a.new)>0.0001)

  a.hat <- Res$a.new
  ## Refit Coefficient
  if((length(ind.sel)!=0 & is.refit == T)){
    Res.refit <- cscca(Y=Y[,ind.sel,drop=F],
                       View.ind=View.ind[ind.sel],
                       lambda.seq=rep(0,K.View),
                       View.type = View.type,
                       a.old = Res$a.new[ind.sel],
                       eps.stop = eps.stop,
                       max.step = max.step,
                       eps = eps,
                       T.step = T.step)
    a.hat[ind.sel] <- Res.refit$a.new
  }
  a.hat.p <- (abs(a.hat)>0.0001)*1


  #===== Find Opt Lambda ========
  lam.opt.trgt <- lam.opt
  a.hat.opt.trgt <- a.hat

  results <- list(a.hat.opt.trgt=a.hat.opt.trgt,
                  lam.opt.trgt=lam.opt.trgt,
                  run=run,
                  CV.mat=CV.mat
  )
  return(results)
}


Hyper_Para.determine <- function(Y,View.ind,
                                 lambda.start.seq=NULL,
                                 View.type=NULL,
                                 eps.stop = 1e-04,
                                 max.step = 30,
                                 eps = 1e-04,
                                 T.step = 1){
  K.View <- length(unique(View.ind))
  Opt.Path <- matrix(nrow=0,ncol=(K.View)*2)
  if(is.null(lambda.start.seq)){
    lambda.start.seq <- rep(4,K.View)
  }

  if(is.null(View.type)){
    View.type <- rep("O",K.View)
  }

  ### Initialize
  ## Find the hyperparameter inducing non-zero coefficient
  a.new.upper.init <- rep(0,length(View.ind))
  while(all(a.new.upper.init==0)){
    Res.start <- cscca(Y=Y,View.ind=View.ind,
                       lambda.seq=lambda.start.seq,
                       View.type=View.type)
    a.new.upper.init <- Res.start$a.new
    Opt.Path <- rbind(Opt.Path,
                      c(lambda.start.seq,
                        sapply(1:K.View,function(k){sum(a.new.upper.init[View.ind==k]!=0)}))
    )
    if(all(a.new.upper.init==0)){
      lambda.start.seq <- lambda.start.seq/2
    }
  }

  lambda.start.upper.seq <- lambda.start.seq

  ## Find the hyperparameter inducing coefficient without zero
  a.new.lower.init <- a.new.upper.init
  while(any(a.new.lower.init==0)){
    Res.start <- cscca(Y=Y,View.ind=View.ind,
                       lambda.seq=lambda.start.seq,
                       View.type=View.type)
    a.new.lower.init <- Res.start$a.new

    Opt.Path <- rbind(Opt.Path,
                      c(lambda.start.seq,
                        sapply(1:K.View,function(k){sum(a.new.lower.init[View.ind==k]!=0)}))
    )
    if(any(a.new.lower.init!=0)){
      lambda.start.seq <- lambda.start.seq/2
    }
  }

  lambda.start.lower.seq <- lambda.start.seq

  ### Precisely Tuning
  ## Tuning the lower bound of hyperparameter sequence
  for(k in 1:(K.View)){
    while(all(a.new.lower.init[View.ind==k]!=0)){
      lambda.start.lower.seq[k] <- lambda.start.lower.seq[k]*1.1
      Res.start <- cscca(Y=Y,View.ind=View.ind,
                         lambda.seq=lambda.start.lower.seq,
                         View.type=View.type)

      a.new.lower.init <- Res.start$a.new

      Opt.Path <- rbind(Opt.Path,
                        c(lambda.start.lower.seq,
                          sapply(1:K.View,function(k){sum(a.new.lower.init[View.ind==k]!=0)}))
      )
    }
  }

  ## Tuning the upper bound of hyperparameter sequence
  for(k in 1:(K.View)){
    while(any(a.new.upper.init[View.ind==k]!=0)){
      lambda.start.upper.seq[k] <- lambda.start.upper.seq[k]*1.1

      lambda.seq.in <- lambda.start.upper.seq
      Res.start <- cscca(Y=Y,View.ind=View.ind,
                         lambda.seq=lambda.seq.in,
                         View.type=View.type)
      a.new.upper.init <- Res.start$a.new
      Opt.Path <- rbind(Opt.Path,
                        c(lambda.start.upper.seq,
                          sapply(1:K.View,function(k){sum(a.new.upper.init[View.ind==k]!=0)}))
      )
    }
  }

  return(list(lambda.upper.seq = lambda.start.upper.seq,
              lambda.lower.seq = lambda.start.lower.seq,
              Opt.Path=Opt.Path))
}

#' Data Generating Process (Omics Data versus Compositional data)
#'
#' @param seed an integer for the initial seed.
#' @param n an integer representing the sample size.
#' @param p an integer representing the feature size of the omics dataset.
#' @param q an integer representing the feature size of the compositional dataset.
#' @param sigma.nu a numerial value representing the strength of correlation.
#' @param sigma.eps a numerical value representing the strength of noise.
#' @param omega_X a p vector representing the coefficient for the omics data.
#' @param omega_Y a q vector representing the coefficient for the compositional data.
#'
#' @return A list containing the following elements: (a) \code{Y}: a n*(2p) matrix representing the full observations; (b) \code{View.ind}: a 2p integer vector indicating the classes of features. The features with the same View.ind is in the same class; (c) \code{omega} a 2p vector representing the true coefficients.
#' @export
#'
#' @examples
#' library(dplyr)
#' n <- 200
#' p <- q <- 100
#' sigma.nu <- 5
#' sigma.eps <- 1
#' omega_X <- 0.85*c(rep(1/10,9),-9/10,rep(0,p-10))
#' omega_Y <- 0.85*c(seq(0.08,0.12,length = 10),rep(0,q-10))
#' Data1 <- DGP_OC(seed=10,n,p,q,sigma.nu,sigma.eps,omega_X,omega_Y)
DGP_OC <- function(seed=10,n,p,q,sigma.nu,sigma.eps,omega_X,omega_Y){
  set.seed(seed)
  # generate nu
  nu <- rnorm(n,sd=sigma.nu)

  # generate U: compositional
  logX <- lapply(as.list(nu),function(nui){nui*omega_X+rnorm(p,sd=sigma.eps)})
  U <- do.call("rbind",logX) #%>% exp()
  U.out <- t(apply(U,1,function(x){scale(x,scale=F)}))
  U.out <- apply(U,2,function(x){scale(x,scale=F)})

  # generate Y: non-compositional
  Y <- lapply(as.list(nu),function(nui){nui*omega_Y+rnorm(q,sd=sigma.eps)})
  Y <- do.call("rbind",Y)
  Y.out <- apply(Y,2,function(x){scale(x,scale=F)})


  #=============================================
  # We re-orgnize data into a multi-view form
  #=============================================

  Y.big <- cbind(Y.out,U.out)
  View.ind <- c(rep(1,p),rep(2,p))
  omega.out <- c(omega_Y,omega_X)

  return(list(Y=Y.big,View.ind=View.ind,omega=omega.out))
}
