
#from dePaula

recoverNetwork <- function(data,lambda,exoeffects=1,docv=0,timeweights=1,eta=0.05,Wder=NA,adaestmethod="cvxr",dopostols=0,Wfixed=-1){

  # settings ==================================================================

  param <- list()

  param$opt$outermeth='L-BFGS-B'
  param$opt$outeriter=100
  param$opt$en$estmethod='glmnet-coupled'
  param$opt$en$rsn=1
  param$opt$adaen$estmethod=adaestmethod
  param$opt$adaen$rsn=1
  param$r <- 0

  param$Wder=Wder
  param$Wfixed=Wfixed
  param$Wfixed0=Wfixed
  
  param$exoeffects=exoeffects
  param$standardise=1

  param$kappa=2.5
  param$eta=eta 
  param$verbose=1

  param$N <- length(unique(data$id))
  param$T <- length(unique(data$time))
  param$K <- length(grep("[x]",colnames(data)))

  param$lambdaL1 <- lambda[1]
  param$lambdaL2 <- lambda[2]
  param$lambdaL1Star <- lambda[3]

  param$fe <- c("id","time")
  
  param$timeweights <- timeweights

  # force main diagonal of Wfixed to be 0
  if (length(param$Wfixed)==1) {
    param$Wfixed <- matrix(NA,nrow=param$N,ncol=param$N)
    }
  diag(param$Wfixed) <- 0

  param$dopostpop <- 0
  param$dopostols <- dopostols
  param$doexternaliv <- 0
  
  # organise data ==============================================================
  
  datalist <- organiseData(data,param)

  data <- datalist$data
  datasq <- datalist$datasq

  result <- list()
  result$sdx <- datalist$sdx
  result$sdxnorm <- datalist$sdxnorm

  # initial conditions =========================================================

  if (param$verbose==1){message('Initial conditions')}

  # coefficients for the non-network model

  formula <- paste0("y~",paste(colnames(data)[grep("[x]",colnames(data))],collapse = "+"))
  betahat <- stats::lm(formula,data=data)
  betahat <- betahat$coefficients[-1]

  if (param$exoeffects==1) {
    p = c(0.5,betahat,rep(0,length(betahat)))
  } else{
    p = c(0.5,betahat)
  }

  # pre-filtering of elements of W
  # note: the function computeDerivatives is not computationally efficient

  if (length(param$Wder)==1) {
    param$Wder <- computeDerivatives(data,datasq,param)
    result$Wder <- param$Wder
  } else {
    message('skip computation of derivative')
  }

  WfixedWder<- (-param$Wder > param$eta)+0
  iX_Wfixed <- which(WfixedWder==0)
  param$Wfixed[iX_Wfixed] <- 0

  if (length(param$Wfixed0)>1) {
    iiX <- which(!is.na(param$Wfixed0))
    param$Wfixed[iiX] <- param$Wfixed0[iiX] 
  }
  
  result$ini$p <- p
  result$ini$Wfixed <- param$Wfixed

  # elastic net ================================================================

  if (param$verbose==1){message('Elastic Net')}

  # penalty

  penalty <- param$Wfixed
  penalty[is.na(penalty)] <- 1
  penalty[penalty==0] <- Inf
  result$en$penalty <- penalty

  # run optimisation

  o1 = outerOptim(p,lambdaL1=param$lambdaL1,lambdaL2=param$lambdaL2,penalty=penalty,estmethod=param$opt$en$estmethod,rsn=param$opt$en$rsn,datasq,param)
  W1 = f(as.matrix(o1[1:(2*param$K+1)]),lambdaL1=param$lambdaL1,lambdaL2=param$lambdaL2,penalty,estmethod=param$opt$en$estmethod,rsn=param$opt$en$rsn,datasq,param,returnW=1)
  p1 = o1[1:(2*param$K+1)]
  colnames(p1)<-NULL
  rownames(p1)<-NULL

  # save results

  result$en$optim$value <- o1$value
  result$en$W <- W1
  result$en$p <- p1

  result$en$rho = p1[1][[1]]
  result$en$beta = p1[2:(param$K+1)][[1]] * result$sdxnorm
  if (param$exoeffects==1){
    result$en$gamma = p1[(param$K+2):(2*param$K+1)][[1]] * result$sdxnorm
  } else {
    result$en$gamma = 0*result$en$beta 
  }

  # adaptive elastic net =======================================================
  
  if (param$verbose==1){message('Adaptive Elastic Net')}
  
  # penalty
  
  penalty <- result$en$W
  penalty[penalty <= 0.05] = 0.05
  penalty <- penalty^(-param$kappa)
  penalty[which(param$Wfixed==0)] <- Inf
  result$adaen$penalty <- penalty

  # run optimisation
  
  o2 = outerOptim(p,lambdaL1=param$lambdaL1Star,lambdaL2=param$lambdaL2,penalty=penalty,estmethod=param$opt$adaen$estmethod,rsn=param$opt$adaen$rsn,datasq,param)
  
  if (param$exoeffects==1){
    W2 = f(as.matrix(o2[1:(2*param$K+1)]),lambdaL1=param$lambdaL1Star,lambdaL2=param$lambdaL2,penalty,estmethod=param$opt$adaen$estmethod,rsn=param$opt$adaen$rsn,datasq,param,returnW=1)
    p2 = o2[1:(2*param$K+1)]
  } else {
    W2 = f(as.matrix(o2[1:(param$K+1)]),lambdaL1=param$lambdaL1Star,lambdaL2=param$lambdaL2,penalty,estmethod=param$opt$adaen$estmethod,rsn=param$opt$adaen$rsn,datasq,param,returnW=1)
    p2 = o2[1:(param$K+1)]
  }
  
  colnames(p2)<-NULL
  rownames(p2)<-NULL
  
  # save results
  
  result$adaen$optim$value <- o2$value
  result$adaen$W <- W2
  result$adaen$p <- p2
  
  result$adaen$rho = p2[1][[1]]
  result$adaen$beta = p2[2:(param$K+1)][[1]] * result$sdxnorm
  if (param$exoeffects==1){
    result$adaen$gamma = p2[(param$K+2):(2*param$K+1)][[1]] * result$sdxnorm
  } else {
    result$adaen$gamma = 0*result$adaen$beta
  }
  
  # unpenalised gmm ============================================================
  
  if (param$verbose==1){message('Unpenalised GMM')}
  
  WfixedSaved <- param$Wfixed
  param$Wfixed <- result$adaen$W
    
  # penalty
  penalty <- param$Wfixed
  penalty[is.na(penalty)] <- 1
  penalty[penalty==0] <- Inf
  
  # run optimisation
  
  o3 = outerOptim(unlist(p2),lambdaL1=0,lambdaL2=0,penalty=penalty,estmethod="givenw",rsn=param$opt$adaen$rsn,datasq,param)
  if (param$exoeffects==1){
    p3 = o3[1:(2*param$K+1)]
  } else {
    p3 = o3[1:(param$K+1)]
  }
  colnames(p3)<-NULL
  rownames(p3)<-NULL
  
  # save results
  result$unpenalisedgmm$optim$value <- o3$value
  result$unpenalisedgmm$p <- p3
  result$unpenalisedgmm$W <- result$adaen$W
  
  result$unpenalisedgmm$rho = p3[1][[1]]
  result$unpenalisedgmm$beta = p3[2:(param$K+1)][[1]] * result$sdxnorm
  if (param$exoeffects==1){
    result$unpenalisedgmm$gamma = p3[(param$K+2):(2*param$K+1)][[1]] * result$sdxnorm
  } else {
    result$unpenalisedgmm$gamma = 0*result$unpenalisedgmm$beta
  }
  
  # asymtotics (work in progress) ==============================================
  
  rho = result$unpenalisedgmm$rho
  beta = unlist(result$unpenalisedgmm$p[2:(param$K+1)])
  if (param$exoeffects==1){
    gamma = unlist(result$unpenalisedgmm$p[(param$K+2):(2*param$K+1)])
  } else {
    gamma = 0*beta
  }
  
  G_rho <- (momentDer_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param,s="rho",1))
  G_beta <- sapply(1:param$K,function(i) momentDer_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param,s="beta",i))
  G_gamma <- sapply(1:param$K,function(i) momentDer_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param,s="gamma",i))
  iX <- which(is.na(WfixedSaved))
  G_W <- sapply(iX,function(i) (momentDer_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param,s="W",iX)))
  
  if (param$exoeffects==1){ 
    
    G <- (cbind(G_rho,G_beta,G_gamma, G_W))
    
  } else {
    
    G <- (cbind(G_rho,G_beta, G_W))
    
  }

  Omega <- moment_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param)
  
  invGtG <- ginv(t(G)%*%G)
  GtG <- t(G)%*%G
  vtheta <- ginv(GtG) %*% t(G)%*%Omega%*%G %*% ginv(GtG)
  vtheta <- vtheta/(param$T)
  
  vtheta <- diag(vtheta)

  result$unpenalisedgmm$asyvar <- vtheta
  
  # post peers-of-peers IV =====================================================

  if (param$dopostpop==1) {
    
    dataIV = popIV(result$unpenalisedgmm$W,data,datasq,param)
    
    if (param$exoeffects==1){
      
      textformula <- paste0("y ~ -1+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"+", paste(names(dataIV)[grep("Wx",names(dataIV))],collapse="+"),"|id+time|Wy~", paste(names(dataIV)[grep("WSqx",names(dataIV))],collapse="+"))
      
    } else {
      
      textformula <- paste0("y ~ -1+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"|id+time|Wy~", paste(names(dataIV)[grep("WSqx",names(dataIV))],collapse="+"))
      
    }
    
    panelest <- fixest::feols(as.formula(textformula),data=dataIV,vcov="hetero")
    panelest_coeftable <- panelest$coeftable
    panelest_etable <- fixest::etable(panelest, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10))
    panelest_1ststage <- fixest::fitstat(panelest,"ivwald")
    
    result$pop$panelest <- panelest
    result$pop$panelest_coeftable <- panelest_coeftable
    result$pop$panelest_etable <- panelest_etable
    result$pop$panelest_1ststage <- panelest_1ststage
    
    result$pop$rho <- panelest$coefficients["fit_Wy"]
    result$pop$beta <- panelest$coefficients[names(dataIV)[grep("x",names(dataIV))][1:param$K]] * result$sdxnorm
    
    if (param$exoeffects==1){
      
      result$pop$gamma  <- panelest$coefficients[names(dataIV)[grep("Wx",names(dataIV))]] * result$sdxnorm
      
    } else {
      
      result$pop$gamma <- 0*result$pop$beta
      
    }
  
  }

  # external IV ================================================================

  if (param$doexternaliv==1) {
    
    dataIV = popIV(result$unpenalisedgmm$W,data,datasq,param)
    
    sqnz <- names(data)[grep("[z]",names(data))]
    
    if (length(sqnz)>0){
      
      if (param$exoeffects==1){
        
        textformula <- paste0("y ~ -1+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"+", paste(names(dataIV)[grep("Wx",names(dataIV))],collapse="+"),"|id+time|Wy~", paste(names(dataIV)[grep("Wz",names(dataIV))],collapse="+"))
        
      } else {
        
        textformula <- paste0("y ~ -1+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"|id+time|Wy~", paste(names(dataIV)[grep("Wz",names(dataIV))],collapse="+"))
        
      }
      
      panelest <- fixest::feols(as.formula(textformula),data=dataIV,vcov="hetero")
      panelest_coeftable <- panelest$coeftable
      panelest_etable <- fixest::etable(panelest, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10))
      panelest_1ststage <- fixest::fitstat(panelest,"ivwald")
      
      result$extiv$panelest_coeftable <- panelest_coeftable
      result$extiv$panelest_etable <- panelest_etable
      result$extiv$panelest_1ststage <- panelest_1ststage
      
      result$extiv$rho <- panelest$coefficients["fit_Wy"]
      result$extiv$beta <- panelest$coefficients[names(dataIV)[grep("x",names(dataIV))][1:param$K]] * result$sdxnorm
      
      if (param$exoeffects==1){
        
        result$extiv$gamma  <- panelest$coefficients[names(dataIV)[grep("Wx",names(dataIV))]] * result$sdxnorm
        
      } else {
        
        result$extiv$gamma <- 0*result$extiv$beta
        
      }
      
    }
    
  }

  # post OLS ===================================================================

  if (param$dopostols==1){
    
    dataIV = popIV(result$unpenalisedgmm$W,data,datasq,param)
    
    if (param$exoeffects==1){
      
      textformula <- paste0("y ~ -1+Wy+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"+", paste(names(dataIV)[grep("Wx",names(dataIV))],collapse="+"),"|id+time")
      
    } else {
      
      textformula <- paste0("y ~ -1+Wy+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"|id+time")
      
    }
    
    panelest <- fixest::feols(as.formula(textformula),data=dataIV,vcov="hetero")
    panelest_coeftable <- panelest$coeftable
    panelest_etable <- fixest::etable(panelest, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10))
    
    result$ols$panelest_coeftable <- panelest_coeftable
    result$ols$panelest_etable <- panelest_etable
    
    result$ols$rho <- panelest$coefficients["Wy"]
    result$ols$beta <- panelest$coefficients[names(dataIV)[grep("x",names(dataIV))][1:param$K]] * result$sdxnorm
    
    if (param$exoeffects==1){
      
      result$ols$gamma  <- panelest$coefficients[names(dataIV)[grep("Wx",names(dataIV))]] * result$sdxnorm
      
    } else {
      
      result$ols$gamma <- 0*result$ols$beta
      
    }
    
  }
  
  # BIC ========================================================================
  
  value <- result$adaen$optim$value
  W <- result$adaen$W
  penalty <- result$adaen$penalty
  value_nopen <- value - (param$lambdaL1Star*sum(abs(W)*penalty,na.rm=TRUE) + param$lambdaL2*sum(W^2))
  
  BIC <- log(value_nopen) + sum(abs(result$adaen$W)>1e-5)*log(param$T)/param$T
  
  result$BIC <- BIC
  
  # cross-validation ===========================================================
  
  if (docv==1) {

    crossvalidation <- cv(data,param)
    result$crossvalidation <- crossvalidation

  }
  
  # return =====================================================================

  return(result)
  
}

# optimisation =================================================================

outerOptim = function(p,lambdaL1,lambdaL2,penalty,estmethod,rsn=0,datasq,param){
  
  # bounds
  
  if (param$exoeffects==1) {
    
    lower <- c(0,rep(-Inf,param$K),rep(-Inf,param$K))
    upper <- c(.95,rep(Inf,param$K),rep(Inf,param$K))
    
  } else {
    
    lower <- c(0,rep(-Inf,param$K))
    upper <- c(.95,rep(Inf,param$K))
    
  }
  
  if (param$opt$outermeth == 'L-BFGS-B') {
    
    o <- optimx::optimx(
      p,
      f,
      lower=lower,
      upper=upper,
      control=list(trace=0,maxit=param$maxit,kkt=FALSE,save.failures=TRUE,pgtol=1e-5),
      hessian=FALSE,
      method=param$opt$outermeth,
      itnmax=param$opt$outeriter,
      lambdaL1=lambdaL1,
      lambdaL2=lambdaL2,
      penalty=penalty,
      estmethod=estmethod,
      rsn=rsn,
      datasq=datasq,
      param=param,
      returnW=0)
    
  }
  
  if (param$opt$outermeth == 'newuoa') {
    
    o <- optimx::optimx(
      p,
      f,
      control=list(trace=0,kkt=FALSE,save.failures=TRUE),
      hessian=FALSE,
      method=param$opt$outermeth,
      lambdaL1=lambdaL1,
      lambdaL2=lambdaL2,
      penalty=penalty,
      estmethod=estmethod,
      rsn=rsn,
      datasq=datasq,
      param=param,
      returnW=0)
    
  }
  
  if (param$opt$outermeth == 'testing') {
    
    o <- optimx::optimx(
      p,
      f,
      method=c('Nelder-Mead','BFGS','L-BFGS-B','CG','nlm','nlminb','spg','ucminf','newuoa','bobyqa','nmkb','hjkb','Rcgmin'),
      lambdaL1=lambdaL1,
      lambdaL2=lambdaL2,
      penalty=penalty,
      estmethod=estmethod,
      rsn=rsn,
      datasq=datasq,
      param=param,
      returnW=0)
    
  }
  
  return(o)
  
}

f <- function(p,lambdaL1,lambdaL2,penalty,estmethod,rsn=0,datasq,param,returnW=0){
  
  # download parameters
  
  rho = p[1]
  beta = p[2:(param$K+1)]
  
  if (param$exoeffects==1){
    gamma = p[(param$K+2):(2*param$K+1)]
  } else {
    gamma = 0*beta
  }
  
  # penalisation parameters
  
  lambda <- lambdaL1+lambdaL2
  if (lambda==0){lambda <- 1}
  alpha <- lambdaL1/lambda
  
  # prepare data
  
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  
  # Y <- y-x*beta
  texteval <- paste0("Y<-datasq$sqy-",paste0("datasq$",sqn,"*beta[",1:length(sqn),"]",collapse="-"))
  eval(parse(text=texteval))
  
  # X <- t(rho*datasq$Y + gamma*datasq$X)
  evaltext <- paste0("X<-t(rho*datasq$sqy+",paste0("datasq$",sqn,"*gamma[",1:length(sqn),"]",collapse="+"),")")
  eval(parse(text=evaltext))
  
  # Z <- t(datasq$X)
  evaltext <- paste0("Z<-cbind(",paste(paste0("t(datasq$",sqn,")"),collapse=","),")")
  eval(parse(text=evaltext))
  
  # allocating space for W and the error term
  
  W <- param$Wfixed
  err <- matrix(ncol=0,nrow=1)
  
  # GLMNET - DECOUPLED
  
  if (estmethod=='glmnet'){
  
    for (i in 1:param$N){
      
      Yi <- Y[i,]
      
      jX <- which(is.na(param$Wfixed[i,]))
      Xi <- X[,jX]
      peni <- as.matrix(penalty[i,jX])
      
      ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
      ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
      
      if (param$r>0) {
        
        ZYia <- rbind(ZYi,param$r)
        ZXia <- rbind(ZXi,param$r)
        
      } else {
        
        ZYia <- ZYi
        ZXia <- ZXi
        
      }
      
      if (length(jX)>=2){
        
        fit <- glmnet::glmnet(ZXia,ZYia,
                                lambda=lambda,
                                alpha=alpha,
                                penalty.factor=peni,
                                lower.limits=0,
                                upper.limits=1,
                                intercept=FALSE,
                                standardize=FALSE)
        
        fit <- as.numeric(fit$beta)
        
        if (rsn==1){fit <- fit / sum(fit)}

        W[i,jX] <- fit
        
        err = c(err,ZYi-ZXi%*%fit)
        
      }
      
      if (length(jX)==1){
        
        W[i,jX] <- 1
        
        err = c(err,ZYi-ZXi*1)
        
      }
      
      if (length(jX)==0){
        
        err = c(err,ZYi)
        
      }
      
    }
    
  }
  
  # GLMNET - COUPLED
  
  if (estmethod=='glmnet-coupled'){
    
    # allocate space for ZY and ZX
    
    ZY <- matrix(ncol=1,nrow=0)
    ZX <- matrix(ncol=0,nrow=0)
    pen <- matrix(ncol=1,nrow=0)
    
    for (i in 1:param$N){
      
      Yi <- Y[i,]
      
      jX <- which(is.na(param$Wfixed[i,]))
      Xi <- X[,jX]
      peni <- as.matrix(penalty[i,jX])
      
      ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
      ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
      
      if (param$r > 0) {
        ZYia <- rbind(ZYi,param$r)
        ZXia <- rbind(ZXi,param$r)
      } else {
        ZYia <- ZYi
        ZXia <- ZXi
      }

      ZY <- rbind(ZY,ZYia) ## (N^2x1)
      ZX <- Matrix::bdiag(ZX,ZXia) ## (N^{2}XN*)
      pen <- rbind(pen,peni)
      
    }
    
    fit <- glmnet::glmnet(ZX,ZY,lambda=lambda,alpha=alpha,penalty.factor=pen,lower.limits=0,upper.limits=1,intercept=FALSE)

    fit = as.numeric(fit$beta)
    
    fit = (1+param$lambdaL2/(param$T^2))*fit ## Debiasing Step
    
    # put coefficients back in W
    W <- t(param$Wfixed)
    W[which(is.na(W))] <- as.numeric(fit)
    W <- t(W)
    
    ## row-sum normalize
    
    if (rsn==1) {
      
      rsnW = rowSums(W)
      
      for (i in 1:param$N){
        
        if (rsnW[i] > 0){W[i,] = W[i,]/rsnW[i]}else{W[i,] = W[i,]/1}
        
      }
    
    }
    
    # compute error and prepare return

    err = ZY-ZX%*%fit
    
    f_obj = sum(err^2) + lambdaL1*sum(abs(W)*penalty,na.rm=TRUE) + lambdaL2*sum(W^2)
    
    # prepare return
    
    if (returnW==0){
      
      return(f_obj)
      
    } else {
      
      return(W)
      
    }
    
  }
  
  # CVXR

  if (estmethod=='cvxr'){
    
    pen <- matrix(ncol=1,nrow=0)    
    ZY <- matrix(ncol=1,nrow=0)
    ZX <- matrix(ncol=0,nrow=0)
    
    theta = CVXR::Variable(sum(is.na(param$Wfixed)))
    
    length_optim_Wfixed = rowSums(is.na(param$Wfixed))
    stop = c(1,cumsum(length_optim_Wfixed))

    constraints = list()
    
    # theta_pos = rowSums(is.na(param$Wfixed))
    # theta_pos = c(0,cumsum(theta_pos))

    for (i in 1:param$N){
      
      Yi <- as.matrix(Y[i,])
      
      jX <- which(is.na(param$Wfixed[i,]))
      Xi <- X[,jX]
      peni <- as.matrix(penalty[i,jX])
      pen <- rbind(pen,peni)
      
      ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
      ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
      
      ZY <- rbind(ZY,ZYi) ## (N^2x1)
      ZX <- Matrix::bdiag(ZX,ZXi) ## (N^{2}XN*)
      
      if (i == 1){constraints = append(constraints,(sum(theta[stop[i]:stop[i+1]]) == 1))}
      if (i > 1){constraints = append(constraints,(sum(theta[(stop[i]+1):stop[i+1]]) == 1))}
      
    }
      
    constraints = append(constraints,(theta >= rep(0,sum(is.na(param$Wfixed)))))
    
    loss <- sum((ZY - ZX %*% theta)^2) / (2 * dim(ZX)[2])
    
    obj <- loss + elastic_reg(theta,lambda=lambda,alpha=alpha,pen=pen)
    
    prob <- CVXR::Problem(CVXR::Minimize(obj),constraints=constraints)
    
    result <- CVXR::psolve(prob,solver="ECOS")
    
    fit = result$getValue(theta)
        
    # put coefficients back in W
    
    W <- t(param$Wfixed)
    W[which(is.na(W))] <- as.numeric(fit)
    W <- t(W)
    
    # objective function
    
    err = ZY-ZX%*%fit

  }
  
  # NO ESTIMATION METHOD - GIVEN W
  
  if (estmethod=='givenw'){
    
    for (i in 1:param$N){
      
      Yi <- Y[i,]
      
      jX <- setdiff(seq(1:param$N),i)
      
      Xi <- X[,jX]
      
      ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
      ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
      
      fit <- W[i,jX]
      
      err = c(err,ZYi-ZXi%*%fit)
      
    }
    
  }
  
  # compute the final objective function
  
  f_obj = sum(err^2) + lambdaL1*sum(abs(W)*penalty,na.rm=TRUE) + lambdaL2*sum(W^2)
  
  # prepare return
  
  if (returnW==0){
    
    return(f_obj)
    
  } else {
    
    return(W)
    
  }
  
}

elastic_reg <- function(theta,lambda,alpha,pen=pen) {
  
  ridge <- (1 - alpha) / 2 * sum(theta^2)
  lasso <- alpha * CVXR::p_norm(pen*theta, 1)
  l <- lambda * (lasso + ridge)
  return(l)
  
}

bic <- function(p,W,datasq,param){
  
  # download parameters
  
  rho = p[1]
  beta = p[2:(param$K+1)]
  
  if (param$exoeffects==1){
    gamma = p[(param$K+2):(2*param$K+1)]
  } else {
    gamma = 0*beta
  }
  
  # prepare data
  
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  
  # Y <- y-x*beta
  texteval <- paste0("Y<-datasq$sqy-",paste0("datasq$",sqn,"*beta[",1:length(sqn),"]",collapse="-"))
  eval(parse(text=texteval))
  
  # X <- t(rho*datasq$Y + gamma*datasq$X)
  evaltext <- paste0("X<-t(rho*datasq$sqy+",paste0("datasq$",sqn,"*gamma[",1:length(sqn),"]",collapse="+"),")")
  eval(parse(text=evaltext))
  
  # Z <- t(datasq$X)
  evaltext <- paste0("Z<-cbind(",paste(paste0("t(datasq$",sqn,")"),collapse=","),")")
  eval(parse(text=evaltext))
  
  # allocate space for ZY and ZX
  
  ZY <- matrix(ncol=1,nrow=0)
  ZX <- matrix(ncol=0,nrow=0)
  pen <- matrix(ncol=1,nrow=0)
  fit <- matrix(ncol=1,nrow=0)
  
  for (i in 1:param$N){
    
    Yi <- Y[i,]
    
    jX <- c(1:param$N)
    Xi <- X[,jX]

    ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
    ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
    
    ZY <- rbind(ZY,ZYi) ## (N^2x1)
    ZX <- Matrix::bdiag(ZX,ZXi) ## (N^{2}XN*)

    fiti <- W[i,jX]
    fit <- c(fit,fiti)
      
  }
  
  # compute error and prepare return
  
  err = ZY-ZX%*%fit
  
  f_obj = sum(err^2)
  
  return(f_obj)
  
}

# data management ==============================================================

organiseData <- function(data,param){
  
  # sort data by id and time 
  
  iX <- base::order(x<-data$time,y<-data$id)
  data <- data[iX,]
  
  uId <- base::sort(unique(data$id))
  uTime <- base::sort(unique(data$time))
  
  id <- NA * data$id
  time <- NA * data$time
  for (i in 1:length(uId)) {id[which(data$id==uId[i])] <- i}
  for (t in 1:length(uTime)) {time[which(data$time==uTime[t])] <- t}
  
  data$id <- id
  data$time <- time
  
  # remove fixed effects 
  
  iX <- grep("[xyz]",colnames(data))
  
  data_dm <- fixest:::demean(X=data[,iX],f=data[,c("id","time")])
  data[,iX] <- data_dm
  
  iXt <- which(data$time==min(data$time))
  data <- data[-iXt,]
  data$time <- data$time-1
  
  # standardise
  
  if (param$standardise==1){
    
    sdx <- matrix(NA,nrow=1,ncol=dim(data)[2])
    colnames(sdx) <- colnames(data)
    
    for (i in iX){
      
      x <- data[,i]
      sdx[,i] <- sqrt(var(data[,i]))
      data[,i] <- x/as.vector(sdx[,i])
      
      sdx_y <- sdx[which(colnames(sdx)=="y")]
      sdx_x <- sdx[grep("[x]",colnames(sdx))]
      sdxnorm <- sdx_y * (sdx_x^(-1))
      
    }
    
  } else {
    
    sdx <- 1
    sdxnorm <- 1
    
  }
  
  # add weights
  
  if (length(param$timeweights)>1) {
    
    for (i in iX){
      
      data[,i] <- data[,i]*sqrt(rep(param$timeweights,each=param$N))
      
    }
    
  }
  
  # square the data
  
  datalist <- list()
  
  datalist$datasq <- squareData(data,param)
  datalist$data <- data
  
  datalist$sdx <- sdx
  datalist$sdxnorm <- sdxnorm
    
  return(datalist)
  
}

squareData <- function(data,param){
  
  cnames <- colnames(data)[grep("[xyz]",colnames(data))]
  
  datasq <- list()
  
  for (i in 1:length(cnames)) {
    
    txteval <- paste0("datasq$sq",cnames[i], " = base::as.matrix(tidyr::pivot_wider(data,id_cols=id,names_from=time,values_from=",cnames[i],")[,-1])")
    eval(parse(text=txteval))
    
    
  }
  
  return(datasq)
  
}

popIV <- function(What,data,datasq,param) {
  
  # creates peer-of-peers instrument
  
  WhatSq <- What %*% What
  
  Wy_pop <- matrix(0,nrow=param$N*(param$T-1))
  sqnx <- names(data)[grep("[x]",names(data))]
  
  for (i in 1:length(sqnx)){
    texteval <- paste0("W",sqnx[i],"_pop <- matrix(0,nrow=param$N*(param$T-1))")
    eval(parse(text=texteval))
    texteval <- paste0("WSq",sqnx[i],"_pop <- matrix(0,nrow=param$N*(param$T-1))")
    eval(parse(text=texteval))
  }
  
  #Wz_pop <- matrix(0,nrow=length(data$x))
  sqnz <- names(data)[grep("[z]",names(data))]
  
  if (length(sqnz)>0) {
    
    for (i in 1:length(sqnz)){
      texteval <- paste0("W",sqnz[i],"_pop <- matrix(0,nrow=param$N*(param$T-1))")
      eval(parse(text=texteval))
    }
    
  }
  
  for (t in 1:(param$T-1)) {
    
    # select the observations for time period t
    Wy_pop[which(data$time==t),] <- What %*% datasq$sqy[,t]
    
    for (i in 1:length(sqnx)) {
      texteval <- paste0("W",sqnx[i],"_pop[which(data$time==t),] <- What %*% datasq$sq",sqnx[i],"[,t]")
      eval(parse(text=texteval))
      texteval <- paste0("WSq",sqnx[i],"_pop[which(data$time==t),] <- WhatSq %*% datasq$sq",sqnx[i],"[,t]")
      eval(parse(text=texteval))
    }
    

    if (length(sqnz)>0){
      for (i in 1:length(sqnz)) {
        texteval <- paste0("W",sqnz[i],"_pop[which(data$time==t),] <- What %*% datasq$sq",sqnz[i],"[,t]")
        eval(parse(text=texteval))
      }
    }
    
  }
  
  data$Wy <- Wy_pop
  for (i in 1:length(sqnx)) {
    texteval <- paste0("data$W",sqnx[i]," <- W",sqnx[i],"_pop")
    eval(parse(text=texteval))
    texteval <- paste0("data$WSq",sqnx[i]," <- WSq",sqnx[i],"_pop")
    eval(parse(text=texteval))
  }
  if (length(sqnz)>0){
  for (i in 1:length(sqnz)) {
    texteval <- paste0("data$W",sqnz[i]," <- W",sqnz[i],"_pop")
    eval(parse(text=texteval))
  }
  }
  
  return(data)
  
}

# moment cond derivative =======================================================

computeDerivatives <- function(data,datasq,param){
  
  # compute beta hat without considering the spatial effects
  formula <- paste0("y~",paste(colnames(data)[grep("[x]",colnames(data))],collapse = "+"),"|id+time")
  panelest <- fixest::feols(as.formula(formula),data)
  betahat <- base::as.vector(panelest$coefficients)
  
  # compute the derivative starting from the empty network, estimatate beta as above, rho = .5 and gamma = 0
  Wstart <- base::matrix(0,param$N,param$N)
  
  #if (param$ncores==0){
    
    # IF SINGLE-CORE THREAD
    g <- momentCond2(rho=.5,beta=betahat,gamma=rep(0,length(betahat)),W=Wstart,datasq,param)
    WDer <- base::lapply(1:(param$N^2),momentDer_i,rho=.5,
                         beta=betahat[1],gamma=0,W=Wstart,datasq=datasq,param=param,g=g)
    WDer <- base::unlist(WDer)
    dim(WDer) <- c(param$N,param$N)
    
  #}
  
  r <- unlist(WDer)
  
  # store
  param$WDer <- WDer
  return(WDer)
  
}

momentCond2 <- function(rho,beta,gamma,W,datasq,param) {
  
  tt <- 1:(dim(datasq$sqy)[2])
  g <- base::sapply(tt,function(t) momentCond2_aux(rho=rho,beta=beta,gamma=gamma,W=W,datasq=datasq,param=param,t=t) )
  g <- Matrix::rowSums(g)/(dim(datasq$sqy)[2])
  
  return(g)
  
}


momentCond2_aux <- function(rho,beta,gamma,W,datasq,param,t){
  
  # select the X
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  xt <- sapply(sqn, function(i) eval(parse(text=paste0("datasq$",i,"[,t]"))))
  sqn <- names(datasq)[grep("[z]",names(datasq))]
  zt <- xt
  dim(zt) <- c(dim(zt)[1]*dim(zt)[2],1)
  
  # calculate the residuals
  residT <- (diag(dim(W)[1]) - rho*W) %*% datasq$sqy[,t] - xt %*% beta - W %*% xt %*% gamma
  
  # compute the moment condition
  g <- base::kronecker(zt,residT)
  
  return(g)
  
}

momentDer_i <- function(rho,beta,gamma,W,datasq,param,g,i) {
  
  if (is.na(param$Wfixed[i])) {
    
    # create matrices used in the computations
    I <- diag(param$N)
    F1 <- solve(I-rho*W)
    F2 <- beta*I+gamma*W
    F3 <- F1*W
    
    dDer <- matrix(0,param$N^2,1)
    A <- matrix(0,param$N,param$N)
    A[i] <- 1
    
    xT <- as.matrix(datasq$sqx1)
    Gk <- (rho * F1 %*% A %*% F1 %*% F2 + gamma * F1 %*% A) %*% xT
    
    for (t in 1:(dim(datasq$sqy)[2])) {
      
      gt <- as.matrix(kronecker(xT[,t],Gk[,t]))
      
      dDer <- dDer + gt
      
    }
    
    dDer <- dDer / (dim(datasq$sqy)[2])
    
    g <- g[1:(length(g)/param$K)]
    WDer <- -2 * t(g) %*% dDer
    
  } else {
    
    WDer <- 0
    
  }
  
  return(WDer)
  
}


# asymptotic distribution ==================================================

momentDer_asyvar <- function(rho,beta,gamma,W,datasq,param,s,i) {
  
  sqn <- names(datasq)[grep("[x]",names(datasq))]

  dDer <- matrix(0,param$N^2*param$K,1)
  
  for (t in 1:(dim(datasq$sqy)[2])) {
    
    yT <- as.matrix(datasq$sqy)
    yT <- yT[,t]
    
    if (s=="rho"){
      
      Gk <- -1 * W %*% yT
      
      }
    
    if (s=="beta"){
      
      evaltext <- paste0("xT <- as.matrix(datasq$sqx",i,")")
      eval(parse(text=evaltext))
      xT <- xT[,t]
      
      Gk <- -xT
      
      }
    
    if (s=="gamma"){
      
      evaltext <- paste0("xT <- as.matrix(datasq$sqx",i,")")
      eval(parse(text=evaltext))
      xT <- xT[,t]
      
      Gk <- -1 * W %*% xT
      
      }
    
    if (s=="W"){
      
      A <- matrix(0,param$N,param$N)
      A[i] <- 1
      
      sqn <- names(datasq)[grep("[x]",names(datasq))]
      
      texteval <- paste0("xT<-cbind(",paste0("datasq$",sqn,"[,t]",collapse=","),")")
      eval(parse(text=texteval))
      
      Gk <- - rho * A %*% yT - A %*% xT %*% gamma
      
    }
    
    texteval <- paste0("xT<-cbind(",paste0("datasq$",sqn,"[,t]",collapse=","),")")
    eval(parse(text=texteval))
    
    dim(xT) <- c(dim(xT)[1]*dim(xT)[2],1)
    
    gt <- as.matrix(kronecker(xT,Gk))
    
    dDer <- dDer + gt
    
  }
  
  dDer <- dDer / dim(datasq$sqy)[2]
  
  return(dDer)
  
}

moment_asyvar <- function(rho,beta,gamma,W,datasq,param,i) {
  
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  
  d <- matrix(0,param$N^2*param$K,param$N^2*param$K)
  
  for (t in 1:(dim(datasq$sqy)[2])) {
    
    texteval <- paste0("xT<-cbind(",paste0("datasq$",sqn,"[,t]",collapse=","),")")
    eval(parse(text=texteval))
    
    yT <- as.matrix(datasq$sqy)
    yT <- yT[,t]
    
    Gk <- yT - rho*W%*%yT - xT%*%beta - W%*%xT%*%gamma
    
    dim(xT) <- c(dim(xT)[1]*dim(xT)[2],1)
    
    gt <- as.matrix(kronecker(xT,Gk))
    
    d <- d + gt %*% t(gt)
    
  }
  
  d <- d / dim(datasq$sqy)[2]
  
  return(d)
  
}


# cross-validation =============================================================

cv <- function(data,param){
  
  # compute to different folds
  
  l <- sapply(1:4,function(f) {cvinner(data,param,fold=f)})
  return(l)

}


cvinner <- function(data,param,fold) {
  
  # split points
  
  splitpoints <- seq(from=min(data$time),to=(max(data$time)+1),length.out=6)
  splitpoints <- splitpoints[-1]
  
  st <- splitpoints[fold]
  en <- splitpoints[fold+1]
  
  # slice the data
  
  data_estimation <- data %>% filter(time<st)
  
  data_evaluation <- data %>% filter(time>=st & time<en)
  data_estimation$time <- data_estimation$time-min(data_estimation$time)+1 
  
  # run procedure
  
  lambda <- c(param$lambdaL1,param$lambdaL2,param$lambdaL1Star)
  res <- recoverNetwork(data_estimation,lambda,docv=0)
  
  # compute msqpe
  
  W <- res$adaen$W
  rho <- res$adaen$rho
  beta <- res$adaen$beta
  gamma <- res$adaen$gamma
  
  # compute msqpe 
  
  msqpe <- 0
  
  for (t in unique(data_evaluation$time)) {
    
    yt <- data_evaluation %>% filter(time==t) %>% dplyr::select(y)
    xt <- data_evaluation %>% filter(time==t) %>% dplyr::select(starts_with("x"))
    
    sqn <- names(data)[grep("[x]",names(data))]
    
    fitt <- rho * W %*% as.matrix(yt) + as.matrix(xt) %*% as.matrix(beta) + W %*% as.matrix(xt) %*% as.matrix(gamma)
    errort <- sum((as.matrix(yt)-fitt)^2)
    
    msqpe <- msqpe + errort
    
  }
  
  return(msqpe)
  
}




