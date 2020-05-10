
####################################################################################################################
################### AM ALGORITHM AS IN HAARIO 2001 #################################################################


run_AM <- function(n=10000, x0=as.vector(c(0,0)), C0=diag(2), t0=1000, sd=.5, epsilon=0.1, target.dens, ...){
  
  X <- matrix(NA, nrow=length(x0), ncol=n)
  X[,1] <- x0
  
  Ctplusone <- C0
  Id <- diag(nrow(X))
  
  # bias! but what else? also use x0?
  rowm <- x0
  
  U <- log(runif(n))
  
  accpt.binar <- rep(NA, times=n)
  cov.out <- matrix(NA, nrow=length(x0)^2, ncol=n)
  cov.out[,1:t0] = as.vector(C0)
  
  for (i in 2:n){
    
    # draw candidate y from adaptive proposal dist.    
    Y = mvrnorm(mu=as.vector(X[,i-1]), Sigma=if(i<=t0) {C0} else {Ctplusone})
    
    if(i>t0){cov.out[,i] <- as.vector(Ctplusone)}
    # is alpha actually needed? otherwise leave min away and leave as log, also transform U as log
    #alpha <- min(1, (target.dens(x=Y, ...) / target.dens(x=X[,i-1], ...)))
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
    # update Ct as used in the next iteration
    Ctplusone <- ((i-1)/i)*Ctplusone + (sd/i)  * (i * tcrossprod(rowm.minus.one)
                                                  - (i+1) * tcrossprod(rowm)
                                                  + tcrossprod(X[,i]) + epsilon*Id
                                                  ) 
    
    
    #if(i%%10000==0){print(i); flush.console()} # print iteration to show progress
    
  }
  
  return(output = list(X=X 
                      ,accpt.binar=accpt.binar
                      ,cov.out=cov.out)
        )
}


####################################################################################################################
################### EXAMPLES #######################################################################################

# sd chosen as 2.48^2/d as proposed by Gelman 1996
#first = run_AM(n=100000, x0=as.vector(c(0.3,0.3)), C0=diag(2), t0=1000, sd=(2.38^2 / 2), epsilon=0.1, target.dens=b_target, A=30, B=30, C1=3, C2=3)
#first = run_AM(n=100000, x0=as.vector(c(0,0)), C0=diag(2), t0=2000, sd=(2.38^2 / 2), epsilon=0.1, target.dens=t_target, v=3, Sigma=diag(2), mu=c(0,0))
#first = run_AM(n=100000, x0=as.vector(c(0,0)), C0=diag(2), t0=1000, sd=(2.38^2 / 2), epsilon=0.1, target.dens=d_target, radius=3, sigma.sq=0.01)
#mean(first$X[1,], na.rm=TRUE)

#A = 30
#B=30
#C1=C2=3
# mu for x2|x1 
#mean((B*first$X[1,] + C1 )/ (A*first$X[1,]^2 + 1), na.rm=T)

# error
# aux = E[x2|x1]
#aux = (B*first$X[1,] + C1 )/ (A*first$X[1,]^2 + 1)
#dif <- (first$X[2,]) - ((B*aux + C1 )/ (A*aux^2 + 1))
# var of dif corresponds to variance of x1|x2
#var(dif)
#1/(A*mean(first$X[2,], na.rm=T)^2 + 1)

#var(first$X[2,] - ((B*first$X[1,]+C2)/(A*first$X[1,]^2 + 1)))

#mean(first$X[2,] - ((B*first$X[1,]+C2)/(A*(var(first$X[1,]) + mean(first$X[1,])^2) + 1)))

# COV
#th.cov = mean(first$X[1,]*first$X[2,] - ((first$X[1,]*(B*first$X[1,] + C2))/(A*first$X[1,]^2 + 1)) - (((first$X[2,]*(B*first$X[2,] + C1))/(A*first$X[2,]^2 + 1))) + (((B*first$X[2,] + 1)*(B*first$X[1,] + 1))/((A*first$X[1,]^2 +1)*(A*first$X[2,]^2 + 1))))

