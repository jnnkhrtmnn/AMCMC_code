####################################################################################################################
################### XN ALGORITHM WITH FINE_TUNING OPTIONS ##########################################################


run_XNfine <- function(n=10000, x0=as.vector(c(0,0)), C0=diag(2), S=2000, h=100, c=.5, target.dens, ...){
  
  X <- matrix(NA, nrow=length(x0), ncol=n)
  X[,1] <- rowm <- x0
  
  k <- 1
  Ct <- C0
  
  U <- log(runif(n))
  
  accpt.binar <- rep(NA, times=n)
  cov.out <- matrix(NA, nrow=length(x0)^2, ncol=n)
  cov.out[,1:1000] = as.vector(C0)
  
  for (i in 2:n){
    
    # draw candidate y from adaptive proposal dist.    
    Y = mvrnorm(mu=as.vector(X[,i-1]), Sigma=if(i<1000){C0} else {k*Ct})
    
    if(i>1000){cov.out[,i] <- as.vector(k*Ct)}
    
    alpha <-  (target.dens(x=Y, ...) - target.dens(x=X[,i-1], ...))
    
    # set new value and update scale
    if(U[i] < alpha){
      X[,i] <- Y
      accpt.binar[i] =  1
      if(i > 1000 & i<=S){k = k + 3*k/(h*i^c)}
    } else {
      X[,i] <- X[,i-1]
      accpt.binar[i] =  0
      if(i > 1000 & i<=S){k = k - k/(h*i^c)}
    }
    
    # update shape
    rowm.minus.one <- rowm
    rowm <- ((i-1)/i) * rowm.minus.one + (1/i) * X[,i]          
    Ct <- ((i-1)/i)*Ct + (1/i)  * (i * tcrossprod(rowm.minus.one)
                                   - (i+1) * tcrossprod(rowm)
                                   + tcrossprod(X[,i])
    )
    
  }
  
  return(output = list(X=X 
                       ,accpt.binar=accpt.binar
                       ,cov.out=cov.out)
  )
}