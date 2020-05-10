
###################################################################################################################
################### ADAPTIVE MCMC IN PRACTICE #####################################################################
###################################################################################################################

####################################################################################################################
################### LOAD LIBRARIES AND UDFS ########################################################################

library(mvtnorm)
library(MASS)
library(compiler)
enableJIT(3)

# find working directory
file_path <- getwd()

# load UDFs
source(paste(file_path, "AM_algorithm.R", sep="/"))
source(paste(file_path, "AP_algorithm.R", sep="/"))
source(paste(file_path, "AccAM_algorithm.R", sep="/"))
source(paste(file_path, "RR_algorithm.R", sep="/"))
source(paste(file_path, "XN_algorithm.R", sep="/"))
source(paste(file_path, "XNfine_algorithm.R", sep="/"))
source(paste(file_path, "IndAM_algorithm.R", sep="/"))
source(paste(file_path, "target_densities.R", sep="/"))
source(paste(file_path, "chain_diagnostics.R", sep="/"))
source(paste(file_path, "final_test_design.R", sep="/"))


####################################################################################################################
################### TEST EXECUTION #################################################################################

results51d <- run_experiment(R=50
                            ,radius=5
                            ,sigma.sq=.1
                            ,A=30, B=30, C1=3, C2=3
                            ,mu=c(0,0), v=3, Sigma=matrix(c(1,0,0,1), nrow=2)
                            ,n=105000, x0=c(0,0)
                            ,C0=diag(2), t0=1000, sd=(2.38^2 /2), epsilon=0.1
                            ,H=2000
                            ,beta=0.1
                            ,rescale.freq=100
                            ,S=2000
                            ,recalc.freq=100)

save(results51d, file="results51d.RData") #saves the files


results52  <- run_experiment(R=50 
                            ,radius=4 
                            ,sigma.sq=.01
                            ,A=30, B=30, C1=3, C2=3
                            ,mu=c(0,0), v=3, Sigma=matrix(c(3,0.5,0.5,3), nrow=2)
                            ,n=22000, x0=c(20,20)
                            ,C0=diag(2), t0=1000, sd=(2.38^2 /2), epsilon=0.1
                            ,H=2000
                            ,beta=0.1
                            ,rescale.freq=100
                            ,S=2000
                            ,recalc.freq=100)

save(results52, file="results52.RData") #saves the files


results53  <- run_experiment(R=50
                            ,radius=10
                            ,sigma.sq=.05
                            ,A=10, B=10, C1=1, C2=1
                            ,mu=c(1,1), v=3, Sigma=matrix(c(3,0.5,0.5,3), nrow=2)
                            ,n=55000, x0=c(20,20)
                            ,C0=diag(2), t0=5000, sd=(2.38^2 /2), epsilon=0.1
                            ,H=5000
                            ,beta=0.1
                            ,rescale.freq=200
                            ,S=5000
                            ,recalc.freq=200)

save(results53, file="results53.RData") #saves the files

####################################################################################################################
################### ANALYZE RESULTS ################################################################################

################## 4.1 #############################################################################################

load(paste(file_path, "results51d.RData", sep="/"))

results = results51d
results$algorithm = substring(results$algorithm, 5, 100L)
attach(results)


# table 4.1
aggregate(kl.div[target.density=="d_target"] ~ algorithm[target.density=="d_target"], data=results, FUN=mean)
aggregate(kl.div[target.density=="d_target"] ~ algorithm[target.density=="d_target"], data=results, FUN=var)


# Use MCMC to find "theortical " variance of doughnut
d_var <- function(x, radius, sigma.sq){
  if (x>0){
    dens = sqrt(pi/(2*sigma.sq)) * x^3 * exp((-1/(2*sigma.sq))*(x - radius)^2)
  }
  else{dens = 0}
  return(log(dens))
}

d.var.est = run_AM(n=200000, x0=as.vector(12), C0=1, t0=10000, sd=(2.38^2),
                   epsilon=0.1, target.dens=d_var, radius=5, sigma.sq=.1)
mean(d.var.est$X)
hist(d.var.est$X)


# figure 4.1
pdf(paste(file_path, "graphics", "res_51_dtarget_meanx1.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(mean.x1[target.density=="d_target"] ~ algorithm[target.density=="d_target"],
        main=expression( 'E[' * 'X'[1]* ']'* ' estimate' * ' doughnut distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=0, col="red")
dev.off()

# figure 4.1
pdf(paste(file_path, "graphics", "res_51_dtarget_meanx2.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(mean.x2[target.density=="d_target"] ~ algorithm[target.density=="d_target"],
        main=expression( 'E[' * 'X'[2]* ']'* ' estimate' * ' doughnut distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=0, col="red")
dev.off()






################## 4.2 #############################################################################################

load(paste(file_path, "results52.RData", sep="/"))

results = results52
results$algorithm = substring(results$algorithm, 5, 100L)
attach(results)


# table 4.2
aggregate(prop.covar.end[target.density=="b_target"] ~ algorithm[target.density=="b_target"], data=results, FUN=mean)
aggregate(prop.covar.end[target.density=="b_target"] ~ algorithm[target.density=="b_target"], data=results, FUN=var)
aggregate(covar.end[target.density=="b_target"] ~ algorithm[target.density=="b_target"], data=results, FUN=mean)

# table 4.3
aggregate(autocorr1[target.density=="b_target"] ~ algorithm[target.density=="b_target"], data=results, FUN=mean)
aggregate(autocorr2[target.density=="b_target"] ~ algorithm[target.density=="b_target"], data=results, FUN=mean)


# figure 4.2
pdf(paste(file_path, "graphics", "res_52_ttarget_meanx1.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(mean.x1[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'E[' * 'X'[1]* ']'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=0, col="red")
dev.off()

# figure 4.2
pdf(paste(file_path, "graphics", "res_52_ttarget_meanx2.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(mean.x2[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'E[' * 'X'[2]* ']'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=0, col="red")
dev.off()

# figure 4.3
pdf(paste(file_path, "graphics", "res_52_ttarget_varx1.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(var1.end[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'Var(' * 'X'[1]* ')'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=9, col="red")
dev.off()

# figure 4.3
pdf(paste(file_path, "graphics", "res_52_ttarget_varx2.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(var2.end[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'Var(' * 'X'[2]* ')'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=9, col="red")
dev.off()

# figure 4.3
pdf(paste(file_path, "graphics", "res_52_ttarget_covar.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(covar.end[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'Cov(' * 'X'[1]*','* 'X'[2]* ')'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=1.5, col="red")
dev.off()

# figure 4.4
pdf(paste(file_path, "graphics", "res_52_btarget_exchain.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
 first = run_AM(n=30000, x0=as.vector(c(30,30)), C0=diag(2), 
                t0=1000, sd=(2.38^2 / 2), target.dens=b_target, 
                A=30, B=30, C1=3, C2=3)

  thinned = thin.n.cut(first, nth=1, burnin=0)
  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example MCMC Chain",
       xlab=expression('X'[1]) , ylab=expression('X'[2]),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(thinned$X[1,], thinned$X[2,], pch=3)
  points(thinned$X[1,1], thinned$X[2,1], pch=1, col="red", cex=3)
dev.off()

# figure 4.4
pdf(paste(file_path, "graphics", "res_52_btarget_propcov.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  plot(thinned$cov.out[2,], type="l", main="Proposal Covariance",
       xlab="Iteration", ylab=expression('Cov(X'[1] * ',' * 'X'[2] *')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2); 
  abline(h=0, col="red", lty=2)
dev.off()

# figure 4.5
pdf(paste(file_path, "graphics", "res_52_btarget_APchain.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  first = run_AP(n=30000, x0=as.vector(c(30,30)), C0=diag(2), H=2000,
                 sd=(2.38^2 / 2), target.dens=b_target,
                 A=30, B=30, C1=3, C2=3)
  thinned = thin.n.cut(first, nth=1, burnin=0)
  plot(thinned$X[1,], thinned$X[2,], type="l", 
       main="Example AP-MCMC Chain", 
       xlab=expression('X'[1]) , ylab=expression('X'[2]),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(thinned$X[1,], thinned$X[2,], pch=3)
  points(thinned$X[1,1], thinned$X[2,1], pch=1, col="red", cex=3)
dev.off()

# figue 4.5
pdf(paste(file_path, "graphics", "res_52_btarget_APpropcov.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  plot(thinned$cov.out[2,], type="l", main="AP Proposal Covariance",
       xlab="Iteration", ylab=expression('Cov(X'[1] * ',' * 'X'[2] *')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2); 
  abline(h=0, col="red", lty=2)
dev.off()

# figure 4.6
pdf(paste(file_path, "graphics", "res_52_btarget_accrates.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  boxplot(acceptrate[target.density=="b_target" & algorithm==c("XN", "AccAM")] 
          ~ algorithm[target.density=="b_target"& algorithm==c("XN", "AccAM")],
          ylim=c(0,0.3), main="Acceptance rates AccAM vs. XN",
          cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2); 
  abline(h=0.234, col="red")
dev.off()



################## 4.3 #############################################################################################

load(paste(file_path, "results53.RData", sep="/"))

results = results53
results$algorithm = substring(results$algorithm, 5, 100L)
attach(results)

# table 4.4
aggregate(acceptrate[target.density=="d_target"] ~ algorithm[target.density=="d_target"], data=results, FUN=mean)
aggregate(msejd[target.density=="d_target"] ~ algorithm[target.density=="d_target"], data=results, FUN=mean)
aggregate(autocorr1[target.density=="d_target"] ~ algorithm[target.density=="d_target"], data=results, FUN=mean)
aggregate(autocorr2[target.density=="d_target"] ~ algorithm[target.density=="d_target"], data=results, FUN=mean)
aggregate(effsampsize1[target.density=="d_target"] ~ algorithm[target.density=="d_target"], data=results, FUN=mean)

# Example runs
ch1 = run_AM(n=55000, x0=as.vector(c(1,1)), C0=diag(2), t0=5000, epsilon=0.1, sd=(2.38^2 / 2), 
             target.dens=d_target, radius=10, sigma.sq=.05)
ch2 = run_XN(n=55000, x0=c(1,1), S=5000, target.dens=d_target, radius=10, sigma.sq=0.05)
ch3 = run_AP(n=55000, x0=as.vector(c(1,1)), C0=diag(2), H=5000, sd=(2.38^2 / 2),
             target.dens=d_target, radius=10, sigma.sq=.05)
ch4 = run_IndAM(n=55000, x0=as.vector(c(1,1)), C0=diag(2), t0=5000, epsilon=0.1,
                recalc.freq=100, target.dens=d_target, radius=10, sigma.sq=0.05)
ch5 = run_AccAM(n=55000, x0=as.vector(c(1,1)), C0=diag(2), t0=5000, epsilon=0.1,
                rescale.freq=500, target.dens=d_target, radius=10, sigma.sq=0.05)
ch6 = run_RR(n=55000, x0=as.vector(c(1,1)), beta=0.1, target.dens=d_target,
             radius=10, sigma.sq=0.05)

# figure 4.7
pdf(paste(file_path, "graphics", "res_53_dtarget_AccAM.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  thinned = thin.n.cut(ch5, nth=10, burnin=5000)
  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example AccAM Chain",
       xlab=expression('X'[1]) , ylab=expression('X'[2]),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()

# figure 4.7
pdf(paste(file_path, "graphics", "res_53_dtarget_AXN.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  thinned = thin.n.cut(ch2, nth=10, burnin=5000)
  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example XN Chain", 
       xlab=expression('X'[1]) , ylab=expression('X'[2]), 
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()

# figure 4.7
pdf(paste(file_path, "graphics", "res_53_dtarget_AM.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  thinned = thin.n.cut(ch1, nth=10, burnin=5000)
  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example AM Chain",
       xlab=expression('X'[1]) , ylab=expression('X'[2]),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()

# figure 4.8
pdf(paste(file_path, "graphics", "res_53_dtarget_propcov.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  plot(ch1$cov.out[2,], type="l", col="black", ylim=c(-60, 50), 
       ylab=expression('Cov(' * 'X'[1]*','* 'X'[2]* ')'), cex.lab=2, 
       cex.axis=2, cex.main=2, cex.sub=2,
       main=expression('Proposal Cov(' * 'X'[1]*','* 'X'[2]* ')'))
  points(ch2$cov.out[2,], type="l", col="purple")
  points(ch3$cov.out[2,], type="l", col="blue")
  points(ch4$cov.out[2,], type="l", col="red")
  points(ch5$cov.out[2,], type="l", col="green")
  points(ch6$cov.out[2,], type="l", col="orange")
  legend("bottomright", legend=c("AM", "XN", "AP", "IndAM", "AccAM", "RR"), 
         lty=c(1,1,1), col=c("black","purple", "blue", "red", "green", "orange"), cex=1.5)
dev.off()

# figure 4.8
pdf(paste(file_path, "graphics", "res_53_dtarget_propvar1.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  plot(ch1$cov.out[1,], type="l", col="black", ylim=c(0, 170),
       ylab=expression('Var(' * 'X'[1]*')'), cex.lab=2, 
       cex.axis=2, cex.main=2, cex.sub=2,
     main=expression('Proposal Var(' * 'X'[1]*')'))
  points(ch2$cov.out[1,], type="l", col="purple")
  points(ch3$cov.out[1,], type="l", col="blue")
  points(ch4$cov.out[1,], type="l", col="red")
  points(ch5$cov.out[1,], type="l", col="green")
  points(ch6$cov.out[1,], type="l", col="orange")
dev.off()

# figure 4.8
pdf(paste(file_path, "graphics", "res_53_dtarget_propvar2.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  plot(ch1$cov.out[4,], type="l", col="black", ylim=c(0, 170),
       ylab=expression('Var(' * 'X'[2]*')'), cex.lab=2,
       cex.axis=2, cex.main=2, cex.sub=2,
       main=expression('Proposal Var(' * 'X'[2]*')'))
  points(ch2$cov.out[4,], type="l", col="purple")
  points(ch3$cov.out[4,], type="l", col="blue")
  points(ch4$cov.out[4,], type="l", col="red")
  points(ch5$cov.out[4,], type="l", col="green")
  points(ch6$cov.out[4,], type="l", col="orange")
dev.off()

# figure 4.9
pdf(paste(file_path, "graphics", "res_53_dtarget_tracesx1.pdf", sep="/"))
  par(mfrow=c(3,1))  
  par(mar=c(5,5,3,2)+0.1)
  plot(ch1$X[1,], type="l", ylab=expression('AM' * ' X'[1]), xlab="Iteration",
       main=expression('Comparison of traces of '* 'X'[1]),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  plot(ch2$X[1,], type="l", ylab=expression('XN' * ' X'[1]), 
       xlab="Iteration", cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  plot(ch5$X[1,], type="l", ylab=expression('AccAM' * ' X'[1]), 
       xlab="Iteration", cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  par(mfrow=c(1,1))  
dev.off()



################## 4.4 #############################################################################################

ex.res = list()
for (i in 0:5){
  ex.res[[i+1]] = run_XNfine(n=10000, x0=c(1,1), C0=diag(2), S=1000, target.dens=b_target, h=100, c=0.5^i, A=30, B=30, C1=1, C2=1)
}

# figure 4.10
pdf(paste(file_path, "graphics", "res_54_btarget_varx1.pdf", sep="/"))
  par(mar=c(7,5,4,2)+0.1)
  plot(ex.res[[1]]$cov.out[1,], type="l", ylim=c(0,2), main=expression('Proposal Var(' * 'X'[1]* ')'), 
       xlab="Iteration", ylab=expression('Proposal Var(' * 'X'[1]* ')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(ex.res[[2]]$cov.out[1,], type="l", col=2)
  points(ex.res[[3]]$cov.out[1,], type="l", col=3)
  points(ex.res[[4]]$cov.out[1,], type="l", col=4)
  points(ex.res[[5]]$cov.out[1,], type="l", col=5)
  points(ex.res[[6]]$cov.out[1,], type="l", col=6)
  leg = 0.5^(0:5)
  legend("topright", legend=c(leg), cex=2, col=c(1:6), lty=1)
dev.off()

# figure 4.10
pdf(paste(file_path, "graphics", "res_54_btarget_varx2.pdf", sep="/"))
  par(mar=c(7,5,4,2)+0.1)
  plot(ex.res[[1]]$cov.out[4,], type="l", ylim=c(0,2), main=expression('Proposal Var(' * 'X'[2]* ')'), 
       xlab="Iteration", ylab=expression('Proposal Var(' * 'X'[2]* ')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(ex.res[[2]]$cov.out[4,], type="l", col=2)
  points(ex.res[[3]]$cov.out[4,], type="l", col=3)
  points(ex.res[[4]]$cov.out[4,], type="l", col=4)
  points(ex.res[[5]]$cov.out[4,], type="l", col=5)
  points(ex.res[[6]]$cov.out[4,], type="l", col=6)
dev.off()

# figure 4.10
pdf(paste(file_path, "graphics", "res_54_btarget_varx1x2.pdf", sep="/"))
  par(mar=c(7,5,4,2)+0.1)
  plot(ex.res[[1]]$cov.out[2,], type="l", ylim=c(-1,.5), main=expression('Proposal Cov(' * 'X'[1]* ','* 'X'[2]* ')'), 
       xlab="Iteration", ylab=expression('Proposal Cov(' * 'X'[1]* ','* 'X'[2]* ')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(ex.res[[2]]$cov.out[2,], type="l", col=2)
  points(ex.res[[3]]$cov.out[2,], type="l", col=3)
  points(ex.res[[4]]$cov.out[2,], type="l", col=4)
  points(ex.res[[5]]$cov.out[2,], type="l", col=5)
  points(ex.res[[6]]$cov.out[2,], type="l", col=6)
dev.off()


ex.res = list()
for (i in 0:5){
  ex.res[[i+1]] = run_XNfine(n=10000, x0=c(20,20), C0=diag(2), S=1000, target.dens=b_target, h=100, c=0.5^i, A=30, B=30, C1=1, C2=1)
}

# figure 4.12
pdf(paste(file_path, "graphics", "res_54_btarget_varx1off.pdf", sep="/"))
  par(mar=c(7,5,4,2)+0.1)
  plot(ex.res[[1]]$cov.out[1,], type="l", ylim=c(0,30), main=expression('Proposal Var(' * 'X'[1]* ')'), 
       xlab="Iteration", ylab=expression('Proposal Var(' * 'X'[1]* ')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(ex.res[[2]]$cov.out[1,], type="l", col=2)
  points(ex.res[[3]]$cov.out[1,], type="l", col=3)
  points(ex.res[[4]]$cov.out[1,], type="l", col=4)
  points(ex.res[[5]]$cov.out[1,], type="l", col=5)
  points(ex.res[[6]]$cov.out[1,], type="l", col=6)
  leg = 0.5^(0:5)
  legend("topright", legend=c(leg), cex=2, col=c(1:6), lty=1)
dev.off()





ex.res = list()
for (i in 0:5){
  ex.res[[i+1]] = run_XNfine(n=10000, x0=c(1,1), C0=diag(2), S=1000, target.dens=b_target, h=10^i, c=0.5, A=30, B=30, C1=1, C2=1)
}

# figure 4.11
pdf(paste(file_path, "graphics", "res_54_btarget_hvarx1.pdf", sep="/"))
  par(mar=c(7,5,4,2)+0.1)
  plot(ex.res[[1]]$cov.out[1,], type="l", ylim=c(0,2), main=expression('Proposal Var(' * 'X'[1]* ')'), 
       xlab="Iteration", ylab=expression('Proposal Var(' * 'X'[1]* ')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(ex.res[[2]]$cov.out[1,], type="l", col=2)
  points(ex.res[[3]]$cov.out[1,], type="l", col=3)
  points(ex.res[[4]]$cov.out[1,], type="l", col=4)
  points(ex.res[[5]]$cov.out[1,], type="l", col=5)
  points(ex.res[[6]]$cov.out[1,], type="l", col=6)
  leg = 10^(0:5)
  legend("topright", legend=c(leg), cex=2, col=c(1:6), lty=1)
dev.off()

# figure 4.11
pdf(paste(file_path, "graphics", "res_54_btarget_hvarx2.pdf", sep="/"))
  par(mar=c(7,5,4,2)+0.1)
  plot(ex.res[[1]]$cov.out[4,], type="l", ylim=c(0,2), main=expression('Proposal Var(' * 'X'[2]* ')'), 
       xlab="Iteration", ylab=expression('Proposal Var(' * 'X'[2]* ')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(ex.res[[2]]$cov.out[4,], type="l", col=2)
  points(ex.res[[3]]$cov.out[4,], type="l", col=3)
  points(ex.res[[4]]$cov.out[4,], type="l", col=4)
  points(ex.res[[5]]$cov.out[4,], type="l", col=5)
  points(ex.res[[6]]$cov.out[4,], type="l", col=6)
dev.off()

# figure 4.11
pdf(paste(file_path, "graphics", "res_54_btarget_hvarx1x2.pdf", sep="/"))
  par(mar=c(7,5,4,2)+0.1)
  plot(ex.res[[1]]$cov.out[2,], type="l", ylim=c(-1,.5), main=expression('Proposal Cov(' * 'X'[1]* ','* 'X'[2]* ')'), 
       xlab="Iteration", ylab=expression('Proposal Cov(' * 'X'[1]* ','* 'X'[2]* ')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(ex.res[[2]]$cov.out[2,], type="l", col=2)
  points(ex.res[[3]]$cov.out[2,], type="l", col=3)
  points(ex.res[[4]]$cov.out[2,], type="l", col=4)
  points(ex.res[[5]]$cov.out[2,], type="l", col=5)
  points(ex.res[[6]]$cov.out[2,], type="l", col=6)
dev.off()





ex.res = list()
for (i in 0:5){
  ex.res[[i+1]] = run_XNfine(n=10000, x0=c(20,20), C0=diag(2), S=1000, target.dens=b_target, h=10^i, c=0.5, A=30, B=30, C1=1, C2=1)
}

# figure 4.12
pdf(paste(file_path, "graphics", "res_54_btarget_hvarx1off.pdf", sep="/"))
  par(mar=c(7,5,4,2)+0.1)
  plot(ex.res[[1]]$cov.out[1,], type="l", ylim=c(0,50), main=expression('Proposal Var(' * 'X'[1]* ')'), 
       xlab="Iteration", ylab=expression('Proposal Var(' * 'X'[1]* ')'),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  points(ex.res[[2]]$cov.out[1,], type="l", col=2)
  points(ex.res[[3]]$cov.out[1,], type="l", col=3)
  points(ex.res[[4]]$cov.out[1,], type="l", col=4)
  points(ex.res[[5]]$cov.out[1,], type="l", col=5)
  points(ex.res[[6]]$cov.out[1,], type="l", col=6)
  leg = 10^(0:5)
  legend("topright", legend=c(leg), cex=2, col=c(1:6), lty=1)
dev.off()


# New experiment
A=B=30
C1=C2=1
S=2000
C0=diag(2)
x0=c(20,20)
n = 22000
c.param = c(0.5^(0:9))
h.param = c(10^(0:9))


C.p=length(c.param)
H.p=length(h.param)
cnt=1
R = 10
cl.names = c("c.p", "h.p", "mean.x1", "mean.x2","var1.end", "var2.end",
             "covar.end", "prop.var.end1", "prop.var.end2", "prop.covar.end",
             "acceptrate", "autocorr1", "autocorr2", "effsampsize1", "effsampsize2", "msejd")

results = as.data.frame(matrix(NA, nrow=C.p*H.p*R, ncol=length(cl.names)))

results[,1:16] = lapply(results[,1:16], as.numeric)

colnames(results) <- cl.names

# c different values of c
for (c.p in 1:C.p){
  # K different algorithms
  for (h.p in 1:H.p){
    # R replications
    for (r in (1:R)){
      
      param.list = list(n=n, x0=x0, target.dens=b_target, A=A,B=B, C1=C1,C2=C2)
      param.list$c = c.param[c.p]
      param.list$h = h.param[h.p]
      
      
      #### Run MCMC algorithm  
      alg.run = do.call(run_XNfine, param.list)
      thinned = thin.n.cut(alg.run, 5, 2000)
      diagnostics = diagnostic_stats(thinned, b_target, A=A, B=B, C1=C1, C2=C2)
      
      results[cnt,] = c(param.list$c
                        ,param.list$h
                        ,diagnostics$meansx
                        ,diagnostics$varx.end
                        ,diagnostics$prop.var.end
                        ,diagnostics$acceptrate
                        ,diagnostics$autocorrs
                        ,diagnostics$effsampsizes
                        ,diagnostics$msejd
      )
      
      print(cnt); flush.console();
      cnt = cnt + 1
    }
  }
}


save(results, file="results54ft.RData") #saves the files
load(paste(file_path, "results54ft.RData", sep="/")) 

# sort results according to c.p and h.p
results = results[order(results$c.p, results$h.p),]

# figure 4.13
pdf(paste(file_path, "graphics", "res_54_btarget_msejd.pdf", sep="/"))
  par(mar=c(5,5,3,2)+0.1)
  marker = 1
  aux.mat = rep(NA, length=100)
  for (i in 1:100){
    aux.mat[i] = mean(results$msejd[marker:(marker+9)])
    marker = marker+10
  }
  cont.z = matrix(aux.mat, nrow=10, ncol=10, byrow=FALSE)
  colnames(cont.z) = rev(c.param)
  rownames(cont.z) = h.param
  filled.contour(cont.z, x =log(rev(c.param)), y = log(h.param), 
                 main="MSEJD", xlab="log(c)", ylab="log(h)",
                 cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()

# figure 4.13
pdf(paste(file_path, "graphics", "res_54_btarget_ac1.pdf", sep="/"))
  par(mar=c(5,5,3,2)+0.1)
  marker = 1
  aux.mat = rep(NA, length=100)
  for (i in 1:100){
    aux.mat[i] = mean(results$autocorr1[marker:(marker+9)])
    marker = marker+10
  }

  cont.z = matrix(aux.mat, nrow=10, ncol=10, byrow=FALSE)
  colnames(cont.z) = rev(c.param)
  rownames(cont.z) = h.param
  filled.contour(cont.z, x =log(rev(c.param)), y = log(h.param), 
                 main=expression('Lag-1 autocorrelation ' * 'X'[1])
                 , xlab="log(c)", ylab="log(h)",
                 cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()

lmod1 <- lm((results$autocorr1) ~ (log(results$c.p) + log(results$h.p))^2, data=results)
summary(lmod1)
plot(lmod1)



################## END IN DOCUMENT #############################################################################################
################## END IN DOCUMENT #############################################################################################
################## END IN DOCUMENT #############################################################################################
################## END IN DOCUMENT #############################################################################################
################## END IN DOCUMENT #############################################################################################











