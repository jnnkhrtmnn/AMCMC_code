
####################################################################################################################
################### AM ALGORITHM WITH ACC.RATE RESCALING ###########################################################


run_AccAM <- function(n=10000, x0=as.vector(c(0,0)), C0=diag(2), t0=1000, epsilon=0.1, rescale.freq=100, target.dens, ...){
  
  X <- matrix(NA, nrow=length(x0), ncol=n)
  X[,1] <- x0
  
  Ctplusone <- C0
  Id <- diag(nrow(X))
  
  # bias! but what else? also use x0?
  rowm <- x0
  
  U <- log(runif(n))
  
  accpt.binar <- rep(NA, times=n)
  
  scaling <- 1
  scalings <- rep(NA, times=(n/rescale.freq))
  
  cov.out <- matrix(NA, nrow=length(x0)^2, ncol=n)
  cov.out[,1:t0] = scaling*as.vector(C0)
  
  for (i in 2:n){
    
    # draw candidate y from adaptive proposal dist.    
    Y = mvrnorm(mu=as.vector(X[,i-1]), Sigma=if(i<=t0) {C0} else {scaling*Ctplusone})
    
    if(i>t0){cov.out[,i] <- as.vector(Ctplusone)*scaling}

    alpha <-  (target.dens(x=Y, ...) - target.dens(x=X[,i-1], ...))
    
    if(U[i] < alpha){
      X[,i] <- Y
      accpt.binar[i] = 1
    } else {
      X[,i] <- X[,i-1]
      accpt.binar[i] = 0
    }
    
      if(i%%rescale.freq==0 && i>t0){
          temp.accrate <- mean(accpt.binar[(i-rescale.freq+1):i], na.rm=TRUE)
          scaling <- max(0.1, scaling + 0.5*(temp.accrate - 0.234))
          scalings[i/rescale.freq] <- scaling 
      }
    
    # update row means in recursive form
    rowm.minus.one <- rowm
    rowm <- ((i-1)/i) * rowm.minus.one + (1/i) * X[,i]          
    # update Ct as used in the next iteration
    Ctplusone <- ((i-1)/i)*Ctplusone + (1/i)  * (i * tcrossprod(rowm.minus.one)
                                                  - (i+1) * tcrossprod(rowm)
                                                  + tcrossprod(X[,i]) + epsilon*Id
                                                ) 
    
  }
  
  return(output = list(X=X 
                       ,accpt.binar=accpt.binar
                       ,scalings=scalings
                       ,cov.out=cov.out)
  )
}


####################################################################################################################
################### EXAMPLES #######################################################################################

# sd chosen as 2.48^2/d as proposed by Gelman 1996
#first = run_AccAM(n=100000, x0=as.vector(c(0.3,0.3)), C0=diag(2), t0=1000, sd=(2.38^2 / 2), epsilon=0.1, rescale.freq=100000, target.dens=b_target, A=20, B=20, C1=3, C2=3)
#first = run_AccAM(n=500000, x0=as.vector(c(0,0)), C0=diag(2), t0=2000, epsilon=0.1, rescale.freq=100, target.dens=t_target, v=3, Sigma=diag(2), mu=c(0,0))
#first2 = run_AM(n=500000, x0=as.vector(c(0,0)), C0=diag(2), t0=2000, epsilon=0.1, sd=2.38^2/2, target.dens=t_target, v=3, Sigma=diag(2), mu=c(0,0))

#first = run_AccAM(n=20000, x0=as.vector(c(.7,.7)), C0=diag(2), t0=10000, epsilon=0.1, rescale.freq=100, target.dens=d_target, radius=5, sigma.sq=.01)

#ts.plot(first$scalings)
#mean(first$accpt.binar, na.rm=T)
#ts.plot(first$cov.out[1,])
#ts.plot(first2$cov.out[1,])

#thinned = thin.n.cut(first2, nth=20, burnin=0)

#plot(thinned$X[1,], thinned$X[2,])
#plot(thinned$X[1,], thinned$X[2,], type="l")

#autocorrs <- apply(first$X, MARGIN=1, FUN=function(ac) c(acf(ac, plot=FALSE)$acf)[2])
#effsampsizes = round(ncol(first$X) * (1-autocorrs)/(1+autocorrs),0)
#print(effsampsizes)
