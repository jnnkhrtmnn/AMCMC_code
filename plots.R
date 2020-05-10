

############################################################################################################################
############# 5.1  #########################################################################################################


load(paste(getwd(), "results51d.RData", sep="/"))

results = results51d
results$algorithm = substring(results$algorithm, 5, 100L)
attach(results)



pdf(paste(getwd(), "graphics", "res_51_dtarget_meanx1.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(mean.x1[target.density=="d_target"] ~ algorithm[target.density=="d_target"],
        main=expression( 'E[' * 'X'[1]* ']'* ' estimate' * ' doughnut distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=0, col="red")
dev.off()

pdf(paste(getwd(), "graphics", "res_51_dtarget_meanx2.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(mean.x2[target.density=="d_target"] ~ algorithm[target.density=="d_target"],
        main=expression( 'E[' * 'X'[2]* ']'* ' estimate' * ' doughnut distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=0, col="red")
dev.off()







############################################################################################################################
###############  5.2  ######################################################################################################

load(paste(getwd(), "results52.RData", sep="/"))

results = results52
results$algorithm = substring(results$algorithm, 5, 100L)
attach(results)

pdf(paste(getwd(), "graphics", "res_52_ttarget_meanx1.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(mean.x1[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'E[' * 'X'[1]* ']'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=0, col="red")
dev.off()

pdf(paste(getwd(), "graphics", "res_52_ttarget_meanx2.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(mean.x2[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'E[' * 'X'[2]* ']'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=0, col="red")
dev.off()

pdf(paste(getwd(), "graphics", "res_52_ttarget_varx1.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(var1.end[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'Var(' * 'X'[1]* ')'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=9, col="red")
dev.off()

pdf(paste(getwd(), "graphics", "res_52_ttarget_varx2.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(var2.end[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'Var(' * 'X'[2]* ')'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=9, col="red")
dev.off()

pdf(paste(getwd(), "graphics", "res_52_ttarget_covar.pdf", sep="/"))
par(mar=c(7,5,4,2)+0.1)
boxplot(covar.end[target.density=="t_target"] ~ algorithm[target.density=="t_target"],
        main=expression( 'Cov(' * 'X'[1]*','* 'X'[2]* ')'* ' estimate' * ' Bivariate t-distribution'),
        cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, las=2);
abline(h=1.5, col="red")
dev.off()

#pdf(paste(getwd(), "graphics", "res_52_btarget_exchain.pdf", sep="/"))
#  par(mar=c(5,5,4,2)+0.1)
# first = run_AM(n=30000, x0=as.vector(c(30,30)), C0=diag(2), t0=1000, sd=(2.38^2 / 2), target.dens=b_target, A=30, B=30, C1=3, C2=3)

#  thinned = thin.n.cut(first, nth=1, burnin=0)
#  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example MCMC Chain", xlab=expression('X'[1]) , ylab=expression('X'[2]),
#       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
#  points(thinned$X[1,], thinned$X[2,], pch=3)
#  points(thinned$X[1,1], thinned$X[2,1], pch=1, col="red", cex=3)
#dev.off()

#pdf(paste(getwd(), "graphics", "res_52_btarget_propcov.pdf", sep="/"))
#  par(mar=c(5,5,4,2)+0.1)
#  plot(thinned$cov.out[2,], type="l", main="Proposal Covariance",
#       xlab="Iteration", ylab=expression('Cov(X'[1] * ',' * 'X'[2] *')'),
#       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2); 
#  abline(h=0, col="red", lty=2)
#dev.off()




#pdf(paste(getwd(), "graphics", "res_52_btarget_APchain.pdf", sep="/"))
#  par(mar=c(5,5,4,2)+0.1)
 #first = run_AP(n=30000, x0=as.vector(c(30,30)), C0=diag(2), H=2000, sd=(2.38^2 / 2), target.dens=b_target, A=30, B=30, C1=3, C2=3)
#  thinned = thin.n.cut(first, nth=1, burnin=0)
#  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example AP-MCMC Chain", xlab=expression('X'[1]) , ylab=expression('X'[2]),
#       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
#  points(thinned$X[1,], thinned$X[2,], pch=3)
#  points(thinned$X[1,1], thinned$X[2,1], pch=1, col="red", cex=3)
#dev.off()

#pdf(paste(getwd(), "graphics", "res_52_btarget_APpropcov.pdf", sep="/"))
#  par(mar=c(5,5,4,2)+0.1)
#  plot(thinned$cov.out[2,], type="l", main="AP Proposal Covariance",
#       xlab="Iteration", ylab=expression('Cov(X'[1] * ',' * 'X'[2] *')'),
#       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2); 
#  abline(h=0, col="red", lty=2)
#dev.off()

pdf(paste(getwd(), "graphics", "res_52_btarget_accrates.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  boxplot(acceptrate[target.density=="b_target" & algorithm==c("XN", "AccAM")] 
          ~ algorithm[target.density=="b_target"& algorithm==c("XN", "AccAM")],
          ylim=c(0,0.3), main="Acceptance rates AccAM vs. XN",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2); 
          abline(h=0.234, col="red")
dev.off()








############################################################################################################################
############  5.3  #########################################################################################################




ch1 = run_AM(n=55000, x0=as.vector(c(1,1)), C0=diag(2), t0=5000, epsilon=0.1, sd=(2.38^2 / 2), target.dens=d_target, radius=10, sigma.sq=.05)
ch2 = run_XN(n=55000, x0=c(1,1), S=5000, target.dens=d_target, radius=10, sigma.sq=0.05)
ch3 = run_AP(n=55000, x0=as.vector(c(1,1)), C0=diag(2), H=5000, sd=(2.38^2 / 2), target.dens=d_target, radius=10, sigma.sq=.05)
ch4 = run_IndAM(n=55000, x0=as.vector(c(1,1)), C0=diag(2), t0=5000, epsilon=0.1, recalc.freq=100, target.dens=d_target, radius=10, sigma.sq=0.05)
ch5 = run_AccAM(n=55000, x0=as.vector(c(1,1)), C0=diag(2), t0=5000, epsilon=0.1, rescale.freq=500, target.dens=d_target, radius=10, sigma.sq=0.05)
ch6 = run_RR(n=55000, x0=as.vector(c(1,1)), beta=0.1, target.dens=d_target, radius=10, sigma.sq=0.05)


pdf(paste(getwd(), "graphics", "res_53_dtarget_AccAM.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  thinned = thin.n.cut(ch5, nth=10, burnin=5000)
  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example AccAM Chain", xlab=expression('X'[1]) , ylab=expression('X'[2]),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()

pdf(paste(getwd(), "graphics", "res_53_dtarget_AM.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  thinned = thin.n.cut(ch1, nth=10, burnin=5000)
  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example AM Chain", xlab=expression('X'[1]) , ylab=expression('X'[2]),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()


pdf(paste(getwd(), "graphics", "res_53_dtarget_IndAM.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  thinned = thin.n.cut(ch4, nth=10, burnin=5000)
  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example IndAM Chain", xlab=expression('X'[1]) , ylab=expression('X'[2]),
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()


pdf(paste(getwd(), "graphics", "res_53_dtarget_AXN.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  thinned = thin.n.cut(ch2, nth=10, burnin=5000)
  plot(thinned$X[1,], thinned$X[2,], type="l", main="Example XN Chain", xlab=expression('X'[1]) , ylab=expression('X'[2]), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()










pdf(paste(getwd(), "graphics", "res_53_dtarget_propcov.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  plot(ch1$cov.out[2,], type="l", col="black", ylim=c(-60, 50), ylab=expression('Cov(' * 'X'[1]*','* 'X'[2]* ')'), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
       main=expression('Proposal Cov(' * 'X'[1]*','* 'X'[2]* ')'))
  points(ch2$cov.out[2,], type="l", col="purple")
  points(ch3$cov.out[2,], type="l", col="blue")
  points(ch4$cov.out[2,], type="l", col="red")
  points(ch5$cov.out[2,], type="l", col="green")
  points(ch6$cov.out[2,], type="l", col="orange")
  legend("bottomright", legend=c("AM", "XN", "AP", "IndAM", "AccAM", "RR"), 
         lty=c(1,1,1), col=c("black","purple", "blue", "red", "green", "orange"), cex=1.5)
dev.off()

pdf(paste(getwd(), "graphics", "res_53_dtarget_propvar1.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)

  plot(ch1$cov.out[1,], type="l", col="black", ylim=c(0, 170), ylab=expression('Var(' * 'X'[1]*')'), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
       main=expression('Proposal Var(' * 'X'[1]*')'))
  points(ch2$cov.out[1,], type="l", col="purple")
  points(ch3$cov.out[1,], type="l", col="blue")
  points(ch4$cov.out[1,], type="l", col="red")
  points(ch5$cov.out[1,], type="l", col="green")
  points(ch6$cov.out[1,], type="l", col="orange")
dev.off()

pdf(paste(getwd(), "graphics", "res_53_dtarget_propvar2.pdf", sep="/"))
  par(mar=c(5,5,4,2)+0.1)
  plot(ch1$cov.out[4,], type="l", col="black", ylim=c(0, 170), ylab=expression('Var(' * 'X'[2]*')'), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
       main=expression('Proposal Var(' * 'X'[2]*')'))
  points(ch2$cov.out[4,], type="l", col="purple")
  points(ch3$cov.out[4,], type="l", col="blue")
  points(ch4$cov.out[4,], type="l", col="red")
  points(ch5$cov.out[4,], type="l", col="green")
  points(ch6$cov.out[4,], type="l", col="orange")
dev.off()


