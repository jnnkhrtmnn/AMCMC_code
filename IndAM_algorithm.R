
####################################################################################################################
################### AM LOCAL INDEPENDENCE SAMPLER ##################################################################

### Further adapt the recalc freq??

run_IndAM <- function(n=10000, x0=as.vector(c(0,0)), C0=diag(2), t0=2000, epsilon=0.1, recalc.freq=100, target.dens, ...){
  
  X <- matrix(NA, nrow=length(x0), ncol=n)
  X[,1] <- x0
  
  Ct <- C0
  Id <- diag(nrow(X))
  rowm <- x0
  
  U <- log(runif(n))
  
  accpt.binar <- rep(NA, times=n)
  
  cov.out <- matrix(NA, nrow=length(x0)^2, ncol=n)
  cov.out[,1:t0] = as.vector(C0)
  
  for (i in 2:n){
    
    if(i%%recalc.freq==0 && i>=recalc.freq){
      
      mean.crf = as.vector(rowMeans(X[,(i-recalc.freq):(i-1)], na.rm=TRUE))
      sigma.crf <- Ct
    }
    
    if (i <= t0){
      
      Y = mvrnorm(mu=as.vector(X[,i-1]), Sigma=C0)
      alpha <-  (target.dens(x=Y, ...) - target.dens(x=X[,i-1], ...))
    
    } else {
   
      Y = mvrnorm(mu=mean.crf, Sigma=sigma.crf)
    
      cov.out[,i] <- as.vector(sigma.crf)
      
      alpha <-  (target.dens(x=Y, ...) - target.dens(x=X[,i-1], ...)
                 + log(dmvnorm(x=X[,i-1], mean=mean.crf, sigma=sigma.crf)) - log(dmvnorm(x=Y, mean=mean.crf, sigma=sigma.crf)))
    
    }
    
    if(U[i] < alpha){
      X[,i] <- Y
      accpt.binar[i] = 1
    } else {
      X[,i] <- X[,i-1]
      accpt.binar[i] = 0
    }

    # update row means in recursive form
    rowm.minus.one <- rowm
    rowm <- ((i-1)/i) * rowm.minus.one + (1/i) * X[,i]          
    # update Ct as used in the next iteration
    Ct <- ((i-1)/i)*Ct + (1/i)  * (i * rowm.minus.one %*% t(rowm.minus.one)
                                                 - (i+1) * rowm %*% t(rowm)
                                                 + X[,i] %*% t(X[,i]) + epsilon*Id
                                                ) 
  }
  
  return(output = list(X=X 
                       ,accpt.binar=accpt.binar
                       ,cov.out=cov.out)
  )
}


####################################################################################################################
################### EXAMPLES #######################################################################################

#first = run_IndAM(n=100000, x0=as.vector(c(0,0)), C0=diag(2), t0=1000, epsilon=0.1, recalc.freq=200, target.dens=d_target, radius=3, sigma.sq=0.01)

#mean(first$accpt.binar, na.rm=T)
#ts.plot(first$cov.out[3,])

#thinned = thin.n.cut(first, nth=10, burnin=0)

#plot(thinned$X[1,], thinned$X[2,])
#plot(thinned$X[1,], thinned$X[2,], type="l")

#autocorrs <- apply(first$X, MARGIN=1, FUN=function(ac) c(acf(ac, plot=FALSE)$acf)[2])
#effsampsizes = round(ncol(first$X) * (1-autocorrs)/(1+autocorrs),0)
#print(effsampsizes)
