# Optimisation set up for Linux machines

require(foreach)
require(parallel)
require(doParallel)
require(doRNG)

source("HAMs.R")
source("CalibrationMethods.R")

BHoptim_min <- function(calib_vec, preset_vec, calib_names, preset_names, obs, burnin, W, tm) {
  
  calib_pars <- data.frame(t(calib_vec))
  colnames(calib_pars) <- calib_names
  
  preset_pars <- data.frame(t(preset_vec))
  colnames(preset_pars) <- preset_names
  
  sim <- BH_model(calib_pars, preset_pars, obs, burnin)
  tmsim <- moments_fx(sim)
  
  # Computing the goodness-of-fit value (i.e, model performance)
  gof <- fiterr_msm_min(tm, W, tmsim)
  
  if (is.na(gof) == TRUE) {
    
    gof <- 1000
    
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

iter <- 3#000

numcores <- detectCores()-1
cl <- makeCluster(numcores, type = "FORK")
registerDoParallel(cl)

starttime <- proc.time()
foreach(s = 1:30) %do% {
  print(s)
  set.seed(s)
  tru <- BH_model(calib_pars, preset_pars, 1000, 5)
  tm = moments_fx(tru)
  W = vcovNW_fx(tru, 7, tm)
  
  y <- foreach(i = 1:30) %dorng% {
    
    calib_m <- matrix(runif(iter*ncalib, min = 0.0001, max = 1), iter, ncalib)
    
    err <- rep(NA, times = iter)
    
    for (rpt in 1:iter) {
      
      calib_vec <- calib_m[rpt,]
      err[rpt] <- BHoptim_min(calib_vec, preset_vec, calib_names, preset_names, obs, burnin, W, tm)
      
    }
    
    opt <- which(err==min(err))
    
    out <- c(s, err[opt], calib_m[opt,])
    
  }
  
  res <- data.frame(matrix(unlist(y), nrow = 30, ncol = 2+ncalib, byrow=TRUE))
  colnames(res) = c("seed", "gof", calib_names)
  
  saveRDS(res, paste("res",s, ".rds", sep=""))
  
}