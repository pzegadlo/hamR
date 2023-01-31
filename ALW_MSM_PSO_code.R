# Optimisation set up for Linux machines

require(pso)
require(foreach)
require(parallel)
require(doParallel)
require(doRNG)

source('HAMs.R')
source('CalibrationMethods.R')

ALWoptim_max <- function(calib_vec, preset_vec, calib_names, preset_names, obs, burnin, W, tm) {
  
  calib_pars <- data.frame(t(calib_vec))
  colnames(calib_pars) <- calib_names
  
  preset_pars <- data.frame(t(preset_vec))
  colnames(preset_pars) <- preset_names
  
  sim <- ALW_model(calib_pars, preset_pars, obs, burnin)
  tmsim <- moments_fx(sim)
  
  # Computing the goodness-of-fit value (i.e, model performance)
  gof <- fiterr_msm_max(tm, W, tmsim)
  
  if (is.na(gof) == TRUE) {
    
    gof <- -1000
    
  }
  
  return(gof)
  
}

calib_names <- c("a", "b", "sig0")
calib_vec <- c(0.0003, 0.0014, 0.03)

preset_names <- c("nf", "nx", "x0", "p0", "pf0")
preset_vec <- c(100, 100, 0, 100, 100)

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
  print(s)
  set.seed(s)
  tru <- ALW_model(calib_pars, preset_pars, 1000, 5)
  tm = moments_fx(tru)
  W = vcovNW_fx(tru, 7, tm)
  
  y <- foreach(i = 1:30, .packages= c("pso")) %dorng% {
    
    psoptim(rep(NA,ncalib), ALWoptim_max, calib_vec, preset_vec, calib_names, preset_names, obs, burnin, W, tm, lower = 0.0001, upper = 1, control = list(maxf=3000, vectorize=FALSE))
    
  }
  
  res = data.frame(matrix(unlist(lapply(y, function(x) c(s, x$value, x$par))), nrow = 30, ncol = 2+ncalib, byrow=TRUE))
  colnames(res) = c('seed', 'gof', calib_names)
  
  saveRDS(res, paste('res',s, '.rds', sep=""))
  
}