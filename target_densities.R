

####################################################################################################################
################### TARGET DENSITIES ###############################################################################

#Add check for boundaries of density, otherwise reject right away

# "doughnut"-shaped

# 2 dims only! include check for this?
# x[1]: x-coordinate
# x[2]: y-coordinate
# r: radius
d_target <- function(x, radius, sigma.sq){
  log.dens <- -(1/(2*sigma.sq))*(sqrt(x[1]^2 + x[2]^2) - radius)^2
  return(log.dens)
}


# "banana"-shaped
#-> see Gelman 1991 paper!
b_target <- function(x, A, B, C1, C2){
  #dens <- as.numeric(x[1] > 0 & x[2] > 0) * 
  #          (exp(-.5 * (A*x[1]^2 * x[2]^2 + x[1]^2 + x[2]^2
  #                      - B*x[1]*x[2] -2*C1*x[1] - 2*C2*x[2])))
  if(x[1] > 0 && x[2] > 0){
    log.dens <- -.5 * (A*x[1]^2 * x[2]^2 + x[1]^2 + x[2]^2
                       - 2*B*x[1]*x[2] -2*C1*x[1] - 2*C2*x[2])
  } else {
    log.dens <- -Inf  
  }
  return(log.dens)
}


# multivariate t-distribution

# mu: location (mean/median/mode)
# v: degrees of freedom
# Sigma: Shape marix (not Variance!!)
# x: p-dimensional input vector / support
t_target <- function(x, mu, v, Sigma){ # change to log density
  # check for Sigma nrow = ncol
  p = length(mu)
  # dens <- (1 + (1/v)*t(x-mu) %*% solve(Sigma) %*% (x-mu))^(-(v+p)/2)
  log.dens <-  (-(v+p)/2) * log((1 + (1/v)*t(x-mu) %*% solve(Sigma) %*% (x-mu)))
  return(log.dens)
}




