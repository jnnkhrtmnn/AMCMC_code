
####################################################################################################################
################### DIAGNOSE THE CHAIN #############################################################################


################### THINNING FUNCTION ##############################################################################

thin.n.cut <- function(obj, nth, burnin){
  # burn-in period?
  obj$X <- obj$X[,(burnin+1):ncol(obj$X)]
  indexes <- seq(1,ncol(obj$X), nth)
  obj$X <- obj$X[,indexes]
  return(obj)
}


################### DIAGNOSTIC STATISTICS ###########################################################################

diagnostic_stats <- function(obj, target.dens, ...){
#diagnostic_stats <- function(obj){
    
  
  meansx <- rowMeans(obj$X, na.rm=TRUE)[1:2]
  varx.end <- var(t(obj$X), na.rm=TRUE)[c(1,4,2)]
  prop.var.end <- c(obj$cov.out[1,ncol(obj$cov.out)], obj$cov.out[4,ncol(obj$cov.out)], obj$cov.out[2,ncol(obj$cov.out)])
  
  acceptrate <- mean(obj$accpt.binar, na.rm=TRUE)
  autocorrs <- apply(obj$X, MARGIN=1, FUN=function(ac) c(acf(ac, plot=FALSE)$acf)[2])
  # assumes first order AR process, for x1 and x2
  effsampsizes = round(ncol(obj$X) * (1-autocorrs)/(1+autocorrs),0)
  #if maxed then minimizes 1 order autcorr (Gelman et al 2010 ESJP paper)
  # MSEJD as in Fearnhead et al 2010

  dif.mat = rbind(obj$X[,-1], obj$X[,-dim(obj$X)[2]])
  msejd = 1/(ncol(obj$X)-1) * sum(apply(dif.mat, MARGIN=2, FUN=function(dif.mat){x = matrix(dif.mat, ncol=2); norm(x, type="2")}))
 
  #################### kl div calculation
  #kde = kde2d(obj$X[1,],y = obj$X[2,], n=200) 
  
  # find closest grid point
  #x.aux <- matrix(NA, nrow=nrow(obj$X),ncol=ncol(obj$X))
  #x.aux[1,] <- sapply(obj$X[1,], FUN=function(X.cord) kde$x[which.min(abs(X.cord-kde$x))])
  #x.aux[2,] <- sapply(obj$X[2,], FUN=function(X.cord) kde$y[which.min(abs(X.cord-kde$y))])
  
  # calculate density of observed (and interpolated) values
  #z.aux <- apply(x.aux, MARGIN=2, FUN=function(x.aux) exp(target.dens(x=x.aux, ...)))
  
  # rescale z.aux, as it is not normalized (rescale to scale of density estimation)
  #z.aux <- (z.aux - min(z.aux)) / (max(z.aux) - min(z.aux)) * (max(kde$z) - min(kde$z)) + min(kde$z)
  
  # find empirical z-values (densities) corresponding to the density estimates of the observed values
  #kde.z <- apply(x.aux, MARGIN=2, FUN=function(x.aux) kde$z[match(x.aux[1],kde$x) , match(x.aux[2],kde$y)])
  
  #com.aux = rbind(z.aux, kde.z)
  
  #com.aux[1,][com.aux[1,]==0] = min(com.aux[1,][com.aux[1,]>0])
  #com.aux[2,][com.aux[2,]==0] = min(com.aux[2,][com.aux[2,]>0])
  
  #kl.div <- (1/ncol(com.aux)) * sum(apply(com.aux , MARGIN=2, FUN=function(com.aux) com.aux[2] * log((com.aux[2])/com.aux[1])))
  
  kl.div=1
    
  return(diagnostic.stats = list(meansx=meansx
                                ,varx.end=varx.end
                                ,prop.var.end=prop.var.end
                                ,acceptrate = acceptrate 
                                ,autocorrs = autocorrs
                                ,effsampsizes = effsampsizes
                                ,msejd = msejd
                                ,kl.div=kl.div)
        )
}

# KL divergence as in Maire et al 2018