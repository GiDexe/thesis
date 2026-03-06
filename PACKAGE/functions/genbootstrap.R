genbootstrap <- function(data,W=NA,rho=NA,beta=NA,gamma=NA,seed,method="wild1") {
  
  if (method=="wild1") {
    
    set.seed(seed)
    
    # remove fixed effects
    
    param <- list()
    
    param$fe <- c("id","time")
    param$standardise <- 0
    param$timeweights <- 1
    
    param$N <- length(unique(data$id))
    param$T <- length(unique(data$time))
    param$K <- length(grep("[x]",colnames(data)))
    
    datalist <- organiseData(data,param)
    
    data <- datalist$data
    datasq <- datalist$datasq
  
    # find residuals
  
    Winv <- solve(diag(dim(W)[1]) - rho*W)
    
    for (t in 1:(param$T-1)) {
      
      # residuals from the original regression
      sqn <- names(datasq)[grep("[x]",names(datasq))]
      xt <- sapply(sqn, function(i) eval(parse(text=paste0("datasq$",i,"[,t]"))))
      residT <- (diag(dim(W)[1]) - rho*W) %*% datasq$sqy[,t] - xt %*% beta - W %*% xt %*% gamma
      
      # wild bootstrap residuals
      # sh <- -1+2*(runif(n=param$N)>0.5)
      sh <- -1+2*(runif(1)>0.5)
      residT <- sh*residT
      
      # put data back in place
      data$y[data$time==t] <- Winv %*% (xt %*% beta + W %*% xt %*% gamma + residT)
      
    }
    
    return(data)
    
    # data <- data %>% filter(time<=param$T-1) 
  
  }
  
  if (method=="block1") {
    
    set.seed(seed)
    
    T <- length(unique(data$time))
    N <- length(unique(data$id))
    Tboot <- sample(seq(1,T),size=T,replace=TRUE)
    
    databoot <- tibble()
    
    for (i in 1:T) {
      databoot <- rbind(databoot,data[which(data$time==Tboot[i]),])
    }
    
    databoot$time  <- rep(1:T,each=N)

    data <- databoot
    
  }
  
}