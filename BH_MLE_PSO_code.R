# Optimisation set up for Linux machines

require(pso)
require(foreach)
require(parallel)
require(doParallel)
require(doRNG)

source('HAMs.R')
source('CalibrationMethods.R')

# For optimisation, the output of ALW model is rescaled by a factor of 100

BHoptim_max <- function(calib_vec, preset_vec, calib_names, preset_names, obs, burnin, tru, h, bestlik, N) {
  
  cd <- rep(0, length(tru))
  
  calib_pars <- data.frame(t(calib_vec))
  colnames(calib_pars) <- calib_names
  
  preset_pars <- data.frame(t(preset_vec))
  colnames(preset_pars) <- preset_names
  
  for (i in 1:N) {
    
    sim <- 100*BH_model(calib_pars, preset_pars, obs, burnin)
    
    cd <- cd + cond_dens(sim, tru, h)
    
  }
  
  cd <- cd / N
  
  # Computing the goodness-of-fit value (i.e, model performance)
  gof <- fiterr_mle_max(cd, bestlik)
  
  if (is.na(gof) == TRUE | gof == Inf) {
    
    gof <- -10 * bestlik
    
  }
  
  return(gof)
  
}

calib_names <- c("g2", "b2", "g3", "b3")
calib_vec <- c(0.6, 0.2, 0.7, -0.2)

preset_names <- c("g1", "b1", "g4", "b4", "r", "beta", "sigma")
preset_vec <- c(0, 0, 1.01, 0, 0.01, 10, 0.04)

calib_pars <- data.frame(t(calib_vec))
colnames(calib_pars) <- calib_names

preset_pars <- data.frame(t(preset_vec))
colnames(preset_pars) <- preset_names

obs <- 1000
burnin <- 5

N <- 1

ncalib <- length(calib_names)

numcores <- detectCores()-1
cl <- makeCluster(numcores, type = "FORK")
registerDoParallel(cl)

starttime <- proc.time()
foreach(s = 1:30) %do% {
  set.seed(s)
  tru <- 100*BH_model(calib_pars, preset_pars, 1000, 5)
  h <- ((4/3)^(1/5)) * sd(tru)
  bestlik <- -sum(log(mgk(tru-tru, h))) 
  
  y <- foreach(cnt = 1:30, .packages= c("pso")) %dorng% { 
    
    psoptim(rep(NA,ncalib), BHoptim_max, calib_pars, preset_pars, obs, burnin, trumod, h, bestlik, N, control = list(maxf=3000, vectorize=FALSE))
    
  }
  
  res = data.frame(matrix(unlist(lapply(y, function(x) c(s, x$value, x$par))), nrow = 30, ncol = 2+ncalib, byrow=TRUE))
  colnames(res) = c('seed', 'gof', calib_names)
  
  saveRDS(res, paste('res',s, '.rds', sep=""))
  
}