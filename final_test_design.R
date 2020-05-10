
####################################################################################################################
################### FINAL TEST DESIGN ##############################################################################

run_experiment <- function(
    #####  Set Test design parameters #####################################
    R=50,                                                                                                                                                                                                                  
    ##### Target density parameters 
    # doughnut
    radius=5, sigma.sq=.1,
    # banana
    A=30, B=30, C1=3, C2=3,
    # t-distribution
    mu=c(0,0), v=3, Sigma=matrix(c(1,0.5,0.5,1), nrow=2),

    ##### Algorithm specific parameters
    n=105000, x0=c(0,0),

    # AM
    C0=diag(2), t0=1000, sd=(2.38^2 /2), epsilon=0.1,  
    # AP
    H=2000,
    # RR
    beta=0.1,
    # AccAM
    rescale.freq=100,
    # XN
    S=2000,
    # Ind AM
    recalc.freq=100
    ##### End Test design parameters #####################################
){

  alg.list=list(`run_AM`, `run_AP`, `run_RR`, `run_AccAM`, `run_XN`, `run_IndAM`)
  names(alg.list) = c("run_AM", "run_AP", "run_RR", "run_AccAM", "run_XN", "run_IndAM")
  K=length(alg.list)

  target.list = list(`d_target`, `b_target`, `t_target`)
  names(target.list) = c("d_target", "b_target", "t_target")
  J=length(target.list)

  cnt=1

  #### Prepare empty results data frame
  cl.names <- c("target.density", "radius", "sigma.sq", "A", "B", "C1", "C2", "mu1", "mu2", "V"
                ,"algorithm", "n", "x01", "x02", "beta", "t0", "sd", "epsilon", "H", "rescale.freq", "recalc.freq", "S"
                ,"mean.x1", "mean.x2","var1.end", "var2.end", "covar.end", "prop.var.end1", "prop.var.end2", "prop.covar.end"
                ,"acceptrate", "autocorr1", "autocorr2", "effsampsize1", "effsampsize2", "msejd", "kl.div")

  results <- as.data.frame(matrix(NA, nrow=(J*K*R), ncol=length(cl.names)))
  results[,-c(1,11)] = lapply(results[,-c(1,11)], as.numeric)
  
  colnames(results) <- cl.names


  ##### Run the experiment

  # J different target densities
  for (j in 1:J){
      # K different algorithms
      for (k in 1:K){
          # R replications
          for (r in (1:R)){
            
              # store results in aux vector
              rs.ax <- as.vector(rep(NA, times=ncol(results)))
              rs.ax[1]=names(target.list)[j]
            
                # Specify which combination of algorithm and density to use
                target.dens.aux = target.list[[j]]
                # params valid for all settings
                rs.ax[11]=names(alg.list)[k]; rs.ax[12]=n; rs.ax[13:14]=x0;
                param.list = list(n=n, x0=x0, target.dens=target.dens.aux)
                diag.param.list = list(target.dens=target.dens.aux)
                
                  # add density specific parameters (and add those to results vector)
                  if(j==1){param.list$radius=rs.ax[2]=radius; param.list$sigma.sq=rs.ax[3]=sigma.sq;
                           diag.param.list$radius=radius; diag.param.list$sigma.sq=sigma.sq}
                  if(j==2){param.list$A=rs.ax[4]=A; param.list$B=rs.ax[5]=B; param.list$C1=rs.ax[6]=C1; param.list$C2=rs.ax[7]=C2; param.list$x0=rs.ax[13:14]=c(1,1);
                           diag.param.list$A=A; diag.param.list$B=B; diag.param.list$C1=C1; diag.param.list$C2=C2;}
                  if(j==3){param.list$mu=rs.ax[8:9]=mu; param.list$v=rs.ax[10]=v; param.list$Sigma=Sigma;
                            diag.param.list$mu=mu; diag.param.list$v=v; diag.param.list$Sigma=Sigma}
                
                  # add algorithm specific parameters (and add those to results vector)
                  if(k==1){param.list$C0=C0; param.list$t0=rs.ax[16]=t0; param.list$sd=rs.ax[17]=sd; param.list$epsilon=rs.ax[18]=epsilon}
                  if(k==2){param.list$C0=C0; param.list$H=rs.ax[19]=H; param.list$sd=rs.ax[17]=sd}
                  if(k==3){param.list$beta=rs.ax[15]=beta}
                  if(k==4){param.list$C0=C0; param.list$t0=rs.ax[16]=t0; param.list$epsilon=rs.ax[18]=epsilon; param.list$rescale.freq=rs.ax[20]=rescale.freq}
                  if(k==5){param.list$C0=C0; param.list$S=rs.ax[22]=S}
                  if(k==6){param.list$C0=C0; param.list$recalc.freq=rs.ax[21]=recalc.freq}
                
                  
                #### Run MCMC algorithm  
                alg.run = do.call(alg.list[[k]], param.list)
            
                #### diagnostic stats
                diag.param.list$obj = thin.n.cut(alg.run, nth=20, burnin=5000)
                #diagnostics= diagnostic_stats(thin.n.cut(alg.run, nth=20, burnin=5000))  
                diagnostics = do.call(diagnostic_stats, diag.param.list)
                
                    rs.ax[23:37]=c(diagnostics$meansx
                                   ,diagnostics$varx.end
                                   ,diagnostics$prop.var.end
                                   ,diagnostics$acceptrate
                                   ,diagnostics$autocorrs
                                   ,diagnostics$effsampsizes
                                   ,diagnostics$msejd
                                   ,diagnostics$kl.div
                                  )
                
            # add results as row to final data frame
            results[cnt,] <- rs.ax
            print(cnt); flush.console()
            cnt = cnt+1
          }
      }
  }
  results[,-c(1,11)] = lapply(results[,-c(1,11)], as.numeric)
  return(results)
}


