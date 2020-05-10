
####################################################################################################################
################### AP ALGORITHM AS IN HAARIO 1999 #################################################################


run_AP <- function(n=10000, x0=as.vector(c(0,0)), C0=diag(2), H=1000, sd=.5, target.dens, ...){
  
  X <- matrix(NA, nrow=length(x0), ncol=n)
  X[,1] <- x0
  
  Ctplusone <- C0*sd
  Id <- diag(nrow(X))
  
  # bias! but what else? also use x0?
  rowm <- x0
  
  U <- log(runif(n))
  
  accpt.binar <- rep(NA, times=n)
  cov.out <- matrix(NA, nrow=length(x0)^2, ncol=n)
  cov.out[,1:H] = as.vector(C0)
  
  for (i in 2:n){
    
    # draw candidate y from adaptive proposal dist.    
    Y = mvrnorm(mu=as.vector(X[,i-1]), Sigma=if(i<=H) {C0} else {Ctplusone})
    
    if(i>H){cov.out[,i] <- as.vector(Ctplusone)}

    
    alpha <-  (target.dens(x=Y, ...) - target.dens(x=X[,i-1], ...))
    
    if(U[i] < alpha){
      X[,i] <- Y
      accpt.binar[i] =  1
    } else {
      X[,i] <- X[,i-1]
      accpt.binar[i] =  0
    }
    
    
    # update row means in recursive form
    if(i==H){
      rowm = rowMeans(X[,1:H])
      Ctplusone = sd * cov(t(X[,1:H]))
      
    } else if(i>H){
        rowm.minus.one <- rowm
        rowm <- rowm.minus.one + (1/H) * X[,i] - (1/H) * X[,i-H]        
        # update Ct as used in the next iteration
    
        Ctplusone <- Ctplusone + sd * (tcrossprod(rowm.minus.one) - tcrossprod(rowm)
                                                + (1/H)* tcrossprod(X[,i]) - (1/H)* tcrossprod(X[,i-H])
                                      )
    }
    
  }
  
  return(output = list(X=X 
                       ,accpt.binar=accpt.binar
                       ,cov.out=cov.out
                       )
        )
}

####################################################################################################################
################### EXAMPLES #######################################################################################

#first = run_AP(n=1000000, x0=as.vector(c(0.3,0.3)), C0=diag(2), H=2000, sd=(2.38^2 / 2), target.dens=b_target, A=1, B=.05, C1=3, C2=3)
#first = run_AP(n=1000000, x0=as.vector(c(0.7,0.7)), C0=diag(2), H=5000, sd=(2.38^2 / 2), target.dens=d_target, r=1)

