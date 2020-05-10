
####################################################################################################################
################### RR ALGORITHM AS IN ROBERTS AND ROSENTHAL 2008 ##################################################


run_RR <- function(n=1000, x0=as.vector(c(0,0)), beta=0.1, target.dens, ...){
  
  X <- matrix(NA, nrow=length(x0), ncol=n)
  X[,1] <- x0
  
  sigma.beta <- (0.1)^2 * diag(length(x0)) / length(x0)
  scaling <- (2.38)^2 / length(x0)
  
  rowm <- x0
  
  U <- log(runif(n))
  
  accpt.binar <-  rep(NA, times=n)
  
  switch <- 1000*length(x0)
  cov.out <- matrix(NA, nrow=length(x0)^2, ncol=n)
  cov.out[,1:switch] = as.vector(sigma.beta)
  
  for (i in 2:n){
    
    # draw candidate y from adaptive proposal dist.
    if(i<=switch){
      Y = mvrnorm(mu=as.vector(X[,i-1]), Sigma=sigma.beta)
    }
    
    if(i>switch){
    Y = (1-beta)*mvrnorm(mu=as.vector(X[,i-1]), Sigma=scaling*sigma.n) + beta*mvrnorm(mu=as.vector(X[,i-1]), Sigma=sigma.beta)
    cov.out[,i] <- as.vector((1-beta)*scaling*sigma.n + beta*sigma.beta)
    }
    

    alpha <-  (target.dens(x=Y, ...) - target.dens(x=X[,i-1], ...))
    
    if(U[i] < alpha){
      X[,i] <- Y
      accpt.binar[i] =  1
    } else {
      X[,i] <- X[,i-1]
      accpt.binar[i] =  0
    }
    
    # update row means in recursive form
    rowm.minus.one <- rowm
    rowm <- ((i-1)/i) * rowm.minus.one + (1/i) * X[,i]          
    
    # update sigma.n as used in the next iteration
    if(i==2){sigma.n <- cov(t(X[,1:2]), use="na.or.complete")} 
    else {
      sigma.n <- ((i-1)/i)*sigma.n + (1/i)  * (i * tcrossprod(rowm.minus.one)
                                                    - (i+1) * tcrossprod(rowm)
                                                    + tcrossprod(X[,i])
                                                ) 
      }
  }
  
  return(output = list(X=X 
                       ,accpt.binar=accpt.binar
                       ,cov.out=cov.out)
  )
}


