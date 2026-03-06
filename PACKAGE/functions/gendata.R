gendata <- function(setting,seed){

  # SET SEED
  base::set.seed(seed)
  
  # load all possible settings
  allsettings <- loadsettings()
  param <- allsettings[allsettings$setting==setting,]
  param$setting <- setting
  
  if (setting<=999){
    
    data <- gendatai(param) 
  
  } else if (setting>1000 & setting<=1999) {
    
    # map the param vectors
    param_start <- param
    param_start$W0 <- strsplit(param$W0,"-")[[1]][1]
    
    param_end <- param
    param_end$W0 <- strsplit(param$W0,"-")[[1]][2]
    
    # generate the data
    data_start <- gendatai(param_start)
    data_end <- gendatai(param_end)
    data_end$time <- data_end$time + max(data_start$time)
    
    # combine
    data <- rbind(data_start,data_end)
    
  } else if (setting>2000) {
    
    data <- gendatai(param) 
    
  }
  
  return(data)

}

# ---------- all settings ----------

gendatai <- function(param){
  
  setting <- param$setting
  
  # create fixed effects
  if (param$FEsetting == "standard") {
    FE <- stats::rnorm(param$N)+1
    TE <- stats::rnorm(param$N)+1
  }
  
  
  # create fixed effects
  if (param$FEsetting == "zero") {
    FE <- 0*stats::rnorm(param$N)+1
    TE <- 0*stats::rnorm(param$N)+1
  }
  
  # define variance-covariance matrix
  if (param$varcov == "orthogonal") {
    Sigma <- base::diag(param$N)
  }
  
  if (param$varcov == "q3") {
    q <- 0.3
    Sigma <- (1-q)*base::diag(param$N)+q
  }
  
  if (param$varcov == "q5") {
    q <- 0.5
    Sigma <- (1-q)*base::diag(param$N)+q
  }
  
  if (param$varcov == "q8") {
    q <- 0.8
    Sigma <- (1-q)*base::diag(param$N)+q
  }
  
  if (param$varcov == "q10") {
    q <- 1.0
    Sigma <- (1-q)*base::diag(param$N)+q
  }
  
  if (setting <= 2000){
    
    # create the reduced-form matrix
    base::eval(base::parse(text=paste('W0 <- W0$W0_', param$W0, '_N',param$N,sep='')))
    Pi0inv <- base::solve(diag(param$N)-param$rho*W0)
    
  } else {
    
    base::eval(base::parse(text=paste('W0erdos <- W0$W0_', 'erdosrenyi', '_N',param$N,sep='')))
    
    W0bilateral <- 0*W0erdos
    for (i in seq(1,29,by=2)){
      W0bilateral[i,i+1]<-1
      W0bilateral[i+1,i]<-1
    }
    
    mix <- (setting-2001)/5
    
    W0 <- (1-mix)*W0erdos + mix*(W0bilateral)
      
    Pi0inv <- base::solve(diag(param$N)-param$rho*W0)
    
  }
  

  # generate data for each time period
  data = dplyr::tibble(y=double(),x=double(),z=double(),id=integer(),time=integer())
  
  for (t in 1:param$T) {
    
    idt <- 1:param$N
    timet <- base::rep(t,param$N)
    
    et <- MASS::mvrnorm(n=1,mu=rep(0,param$N),Sigma)
    xt <- stats::rnorm(param$N*param$K)
    dim(xt) <- c(param$N,param$K)
    
    yt <- Pi0inv %*% xt %*% rep(param$beta,param$K) + Pi0inv %*% W0 %*% xt %*% rep(param$gamma,param$K) + Pi0inv %*% FE + Pi0inv %*% TE + Pi0inv %*% et
    xt <- as_tibble(as.data.frame(xt))
    names(xt) <- paste0('x',1:ncol(xt))
    
    datat <- dplyr::bind_cols(dplyr::tibble(id=idt,time=timet),dplyr::tibble(y=as.vector(yt)),xt)
    data <- base::rbind(data,datat)
    
  }

  # return data
  return(data)

}

# ---------- all settings ----------

loadsettings <- function(setting){
  
  K <- 1

  W0 <- 'erdosrenyi'
  N <- 15
  param <- data.frame(setting=1,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform',stringsAsFactors=FALSE)
  param <- base::rbind(param,data.frame(setting=0002,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0003,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0004,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0005,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0006,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0007,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0008,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0009,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  param <- base::rbind(param,data.frame(setting=0501,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0502,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0503,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0801,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0802,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0803,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 30
  param <- base::rbind(param,data.frame(setting=0011,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0012,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0013,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0014,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0015,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0016,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0017,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0018,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0019,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  param <- base::rbind(param,data.frame(setting=0511,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0512,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0513,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0811,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0812,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0813,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0601,W0=W0,N=N,T=150,K=K,rho=0.1,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0602,W0=W0,N=N,T=150,K=K,rho=0.5,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0603,W0=W0,N=N,T=150,K=K,rho=0.7,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0604,W0=W0,N=N,T=150,K=K,rho=0.9,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0611,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.2,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0612,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.6,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0621,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.3,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0622,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.7,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0631,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q3',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0632,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q5',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0633,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q8',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0634,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q10',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=2001,W0='100erdosrenyi-000bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2002,W0='080erdosrenyi-020bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2003,W0='060erdosrenyi-040bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2004,W0='040erdosrenyi-060bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2005,W0='020erdosrenyi-080bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2006,W0='000erdosrenyi-100bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  
  N <- 50
  param <- base::rbind(param,data.frame(setting=0021,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0022,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0023,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0024,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0025,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0026,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0027,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0028,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0029,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  param <- base::rbind(param,data.frame(setting=0521,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0522,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0523,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0821,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0822,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0823,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 70
  param <- base::rbind(param,data.frame(setting=0031,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0032,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0033,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0034,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0035,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0036,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0037,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0038,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0039,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  W0 <- 'politicalparty'
  N <- 15
  param <- base::rbind(param,data.frame(setting=0041,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0042,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0043,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0044,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0045,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0046,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0047,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0048,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0049,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  param <- base::rbind(param,data.frame(setting=0541,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0542,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0543,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0841,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0842,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0843,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 30
  param <- base::rbind(param,data.frame(setting=0051,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0052,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0053,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0054,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0055,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0056,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0057,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0058,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0059,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  param <- base::rbind(param,data.frame(setting=0551,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0552,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0553,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0851,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0852,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0853,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0651,W0=W0,N=N,T=150,K=K,rho=0.1,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0652,W0=W0,N=N,T=150,K=K,rho=0.5,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0653,W0=W0,N=N,T=150,K=K,rho=0.7,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0654,W0=W0,N=N,T=150,K=K,rho=0.9,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0661,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.2,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0662,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.6,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0671,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.3,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0672,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.7,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0681,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q3',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0682,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q5',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0683,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q8',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0684,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q10',Wknown='no',timeweights='uniform'))

  N <- 50
  param <- base::rbind(param,data.frame(setting=0061,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0062,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0063,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0064,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0065,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0066,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0067,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0068,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0069,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  param <- base::rbind(param,data.frame(setting=0561,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0562,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0563,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0861,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0862,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0863,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 70
  param <- base::rbind(param,data.frame(setting=0071,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0072,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0073,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0074,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0075,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0076,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0077,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0078,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0079,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  W0 <- 'highschool'
  N <- 15
  param <- base::rbind(param,data.frame(setting=0081,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0082,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0083,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0084,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0085,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0086,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0087,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0088,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0089,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  N <- 30
  param <- base::rbind(param,data.frame(setting=0091,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0092,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0093,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0094,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0095,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0096,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0097,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0098,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0099,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  N <- 50
  param <- base::rbind(param,data.frame(setting=0101,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0102,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0103,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0104,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0105,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0106,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0107,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0108,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0109,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  N <- 70
  param <- base::rbind(param,data.frame(setting=0111,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0112,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0113,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0114,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0115,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0116,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0117,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0118,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0119,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  W0 <- 'duflo'
  N <- 15
  param <- base::rbind(param,data.frame(setting=0121,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0122,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0123,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0124,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0125,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0126,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0127,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0128,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0129,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  N <- 30
  param <- base::rbind(param,data.frame(setting=0131,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0132,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0133,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0134,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0135,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0136,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0137,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0138,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0139,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  N <- 50
  param <- base::rbind(param,data.frame(setting=0141,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0142,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0143,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0144,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0145,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0146,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0147,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0148,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0149,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  N <- 70
  param <- base::rbind(param,data.frame(setting=0151,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0152,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0153,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0154,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0155,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0156,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0157,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0158,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0159,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))

  ##
  
  W0 <- 'erdosrenyi-politicalparty'
  N <- 30
  
  param <- base::rbind(param,data.frame(setting=1011,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1012,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1013,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1014,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1015,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1016,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1017,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1018,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1019,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=1011+10,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1012+10,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1013+10,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1014+10,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1015+10,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1016+10,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1017+10,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1018+10,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1019+10,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  
  param <- base::rbind(param,data.frame(setting=1011+20,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1012+20,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1013+20,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1014+20,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1015+20,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1016+20,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1017+20,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1018+20,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1019+20,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  
  W0 <- 'highschool-duflo'
  N <- 70
  
  param <- base::rbind(param,data.frame(setting=1111,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1112,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1113,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1114,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1115,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1116,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1117,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1118,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1119,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=1111+10,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1112+10,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1113+10,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1114+10,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1115+10,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1116+10,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1117+10,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1118+10,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1119+10,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))

  param <- base::rbind(param,data.frame(setting=1111+20,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1112+20,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1113+20,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1114+20,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1115+20,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1116+20,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1117+20,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1118+20,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1119+20,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  
  W0 <- 'geoneighuniform'
  param <- base::rbind(param,data.frame(setting=901,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighmigration'
  param <- base::rbind(param,data.frame(setting=902,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighmigration_trim1'
  param <- base::rbind(param,data.frame(setting=903,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighmigration_trim2'
  param <- base::rbind(param,data.frame(setting=904,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighlength'
  param <- base::rbind(param,data.frame(setting=905,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighdecay'
  param <- base::rbind(param,data.frame(setting=906,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  return(param)

}

# ---------- true coefficients ----------

truecoeffs <- function(param){

  # # load all possible settings
  # allsettings <- loadsettings()
  # param <- allsettings[allsettings$setting==setting,]

  setting <- param$setting
  
  # create the reduced-form matrix
  if (setting <= 2000){
    
    # create the reduced-form matrix
    base::eval(base::parse(text=paste('W0 <- W0$W0_', param$W0, '_N',param$N,sep='')))
    Pi0inv <- base::solve(diag(param$N)-param$rho*W0)
    
  } else {
    
    base::eval(base::parse(text=paste('W0erdos <- W0$W0_', 'erdosrenyi', '_N',param$N,sep='')))
    
    W0bilateral <- 0*W0erdos
    for (i in seq(1,29,by=2)){
      W0bilateral[i,i+1]<-1
      W0bilateral[i+1,i]<-1
    }
    
    mix <- (setting-2001)/5
    
    W0 <- (1-mix)*W0erdos + mix*(W0bilateral)
    
    Pi0inv <- base::solve(diag(param$N)-param$rho*W0)
    
  }
  
  # build list to return
  true <- list()
  true$N <- param$N
  true$T <- param$T
  true$rho <- param$rho
  true$beta <- param$beta
  true$gamma <- param$gamma
  true$W0 <- W0
  true$Pi0inv <- Pi0inv

  return(true)

}

