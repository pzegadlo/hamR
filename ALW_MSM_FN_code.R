require(dplyr)
require(randtoolbox) # for Sobol
require(RANN) # for nearest neighbours
require(foreach)
require(doParallel)
require(doRNG)
require(parallel)

source("HAMs.R")
source("CalibrationMethods.R")

calib_pars <- data.frame(cbind(0.0003, 0.0014, 0.03))
colnames(calib_pars) <- c("a", "b", "sig0")

preset_pars <- data.frame(cbind(100, 100, 0, 100, 100))
colnames(preset_pars) <- c("nf", "nx", "x0", "p0", "pf0")

obs <- 1000
burnin <- 5

trupar <- cbind(calib_pars, preset_pars)

parlist <- colnames(trupar)

npar <- length(parlist)

dimmom <- 7

### Optimisation parameters

# mult: error multiplication for better convergence = 10^mult
# was: 4, changed to 1 for testing
mult <- 1

# hyperparam1: number of the (initial) Sobol draws from the parameter space, def 500
nsob_mult <- 2 #should be 2
nsob_start <- 500 #should be 500

# hyperparam2: number of the best fit values (lowest error) to take for further analysis, def 250
# should this be percentage or number of observations?
bnum <- 250 # should be 250

# hyperparam3: window size for the estimation of mean, s.e., etc., def 30
swin <- 30 # should be 30

# hyperparam5: kappa for the LCB, higher kappa favours exploration, def 5
# maybe we can make it change with the phases of the optimisation
kappa <- 5

# hyperparam6: number of the top windows to take (or percentage), def 20
bwin <- 20L # should be 20

# hyperparam4: numbers of new draws in each top window for each acquisition functions (POI, EI, LCB), def 25, 0 , 0
wacq <- c(25L,0L,0L) # sould be 25

### bwin*wacq need to be nsob?

# Number of iterations in the algorithm
numrep <- 4

# choosing just the POI acquisition function
aq <- 1

# Parameter space definition - can be done outside the loop

parspace <- data.frame(matrix(0, 3, npar))

colnames(parspace) <- parlist
rownames(parspace) <- c("min","max","range")

parspace[1:2,] <- trupar

parspace[1:2,]$a <- c(0.0001,1) # full parameter space
parspace[1:2,]$b <- c(0.0001,1) # full parameter space
parspace[1:2,]$sig0 <- c(0.0001,1) # full parameter space

parspace[3,] <- parspace[2,] - parspace[1,]

numoptpar <- sum(parspace[3,] > 0)

# Calculate parameter values in the current parameter space

nsob <- nsob_mult * nsob_start

parset <- seq(1, nsob, 1)

parrand <- sobol(nsob, numoptpar, scrambling = 3, init = TRUE)

parval <- cbind(parrand,matrix(0,nsob,5))

for (i in 1:nsob) {
  parval[i,] <- unlist(parspace[1,] + t(parspace[3,]) * parval[i,])
}

parval <- data.frame(parval)
colnames(parval) <- parlist

numcores <- detectCores()-1
cl <- makeCluster(numcores, type = "FORK")
registerDoParallel(cl)

foreach(s = 1:30) %do% {
  set.seed(s)
  tru <- ALW_model(calib_pars, preset_pars, 1000, 5)
  tm = moments_fx(tru)
  W = vcovNW_fx(tru, dimmom, tm)
  
  y <- foreach(cnt = 1:30, .packages= c("dplyr", "randtoolbox", "RANN")) %dorng% { 
    
    # notebook setup to keep track of the effects of each step
    
    nbook <- data.frame(matrix(0,numrep+1,20))
    
    colnames(nbook) <- c(paste(rep(c("a","b","sig0"), each=2), rep(c("D","U"),3), sep="_"),
                         "a","b","sig0","MinErr",
                         paste(rep(c("a","b","sig0"), each=2), rep(c("wD","wU"),3), sep="_"),
                         "w_a","w_b","w_sig0","L2")
    rownames(nbook) <- 0:numrep
    
    
    # generate price sequences with different parameter sets
    
    fittable <- data.frame(cbind(parval[,1:3],0))
    
    for (i in 1:nsob) {
      
      sim <- ALW_model(parval[i,1:3], preset_pars, obs, burnin)
      sm <- moments_fx(sim)[1:dimmom]
      fittable[i,4] <- fiterr_msm_min(tm, W, sm)*(10^mult)
      
    }
    
    colnames(fittable) <- c("a","b","sig0","fit")
    
    fittable$fit[is.na(fittable$fit)] <- 2*max(na.omit(fittable$fit)) # replace NAs in the fit function with 2 * max error
    
    # Write to notebook
    
    pos <- "0"
    
    nbook[pos,1] <- parspace$a[1]
    nbook[pos,2] <- parspace$a[2]
    nbook[pos,3] <- parspace$b[1]
    nbook[pos,4] <- parspace$b[2]
    nbook[pos,5] <- parspace$sig0[1]
    nbook[pos,6] <- parspace$sig0[2]
    nbook[pos,7:10] <- fittable[fittable$fit == min(fittable$fit),]
    
    # Optimisation LOOP START
    
    nsob <- nsob_start
    
    parset <- seq(1, nsob, 1)
    
    for (rnd in 1:(numrep+1)) {
      
      # Choose the windows to pick the next group of points from  
      
      ## Points in multispace
      
      mintable <- fittable %>% arrange(fit)
      mintable <- mintable[1:bnum,]
      
      nearest <- nn2(mintable[,1:3], mintable[,1:3], k = swin)$nn.idx
      
      aplot <- data.frame(matrix(ncol = 11, nrow = 0))
      
      for (pp in 1:bnum) {
        
        ttab <- mintable[nearest[pp,],]
        
        tmetrics <- c(min(ttab[,1]), max(ttab[,1]), min(ttab[,2]), max(ttab[,2]), min(ttab[,3]), max(ttab[,3]), mean(ttab[,4]), se(ttab[,4]), min(ttab[,4]))
        
        aplot <- rbind(aplot, tmetrics)
        
      }
      
      colnames(aplot) <- c(paste(rep(c("a","b","sig0"), each=2), rep(c("D","U"),3), sep="_"),"mu","serr","minim")
      
      # minimum from window averages as objective
      objective <- min(aplot$mu)
      
      aplot$gamma <- (objective - aplot$mu) / aplot$serr
      
      
      # calculations for POI only at the moment
      
      aplot$poi <- pnorm(aplot$gamma)
      #plot(rownames(aplot), aplot$poi)
      
      #aplot$ei <- aplot$serr * (aplot$gamma * aplot$poi + dnorm(aplot$gamma))
      #plot(rownames(aplot), aplot$ei)
      
      #aplot$lcb <- aplot$mu - kappa*aplot$serr
      #plot(rownames(aplot), aplot$lcb)
      
      # This ranking should be done taking into account the weighting of the method
      
      aplot$rpoi <- rank(-aplot$poi)
      #aplot$rei <- rank(-aplot$ei)
      #aplot$rlcb <- rank(aplot$lcb)
      
      # Write the best window to notebook
      
      pos <- as.character(rnd-1)
      
      wplot <- aplot %>% arrange(rpoi) #+rei+rlcb)
      
      nbook[pos,11:12] <- wplot[1,1:2]
      nbook[pos,13:14] <- wplot[1,3:4]
      nbook[pos,15:16] <- wplot[1,5:6]
      nbook[pos,17] <- sum(wplot[1,1:2])/2
      nbook[pos,18] <- sum(wplot[1,3:4])/2
      nbook[pos,19] <- sum(wplot[1,5:6])/2
      nbook[pos,20] <- sqrt(sum((nbook[pos,17:19] - trupar[1,1:3])^2))
      
      # Show some current progress stats
      
      # Finish the algorithm when the computational budget is taken, but we have the last window
      
      if (rnd > numrep) { break }
      
      # Get the new points from each window
      
      # Points from POI - EI - LCB
      
      parrand <- sobol(nsob, numoptpar, scrambling = 3, init = FALSE)
      
      parval <- cbind(parrand,matrix(0,nsob,5))
      
      inc <- 0
      
      # choice of acq method
      
      aplot <- wplot
      
      for (candwin in 1:bwin) {
        
        parspace[1:2,]$a <- t(aplot[candwin,1:2])
        parspace[1:2,]$b <- t(aplot[candwin,3:4])
        parspace[1:2,]$sig0 <- t(aplot[candwin,5:6])
        
        parspace[3,] <- parspace[2,] - parspace[1,]
        
        for (i in ((candwin-1)*wacq[aq]+1+inc):(candwin*wacq[aq]+inc)) {
          
          parval[i,] <- unlist(parspace[1,] + t(parspace[3,]) * parval[i,])
        }
        
        parval <- data.frame(parval)
        colnames(parval) <- parlist
        
      }
      
      inc <- i
      
      # generate price sequences with different parameter sets
      
      # calculation of the goodness of fit
      
      # Defining moments in the "sim" datasets
      
      newtable <- data.frame(cbind(parval[,1:3],0))
      
      for (i in 1:nsob) {
        
        sim <- ALW_model(parval[i,1:3], preset_pars, obs, burnin)
        sm <- moments_fx(sim)[1:dimmom]
        
        newtable[i,4] <- fiterr_msm_min(tm, W, sm)*(10^mult)
        
      }
      
      colnames(newtable) <- c("a","b","sig0","fit")
      
      newtable$fit[is.na(newtable$fit)] <- 2*max(na.omit(fittable$fit)) # replace NAs in the fit function with 2 * max error
      
      fittable <- rbind(fittable,newtable)
      
      # Write to notebook
      
      pos <- as.character(rnd)
      
      nbook[pos,1] <- min(parval$a)
      nbook[pos,2] <- max(parval$a)
      nbook[pos,3] <- min(parval$b)
      nbook[pos,4] <- max(parval$b)
      nbook[pos,5] <- min(parval$sig0)
      nbook[pos,6] <- max(parval$sig0)
      nbook[pos,7:10] <- fittable[fittable$fit == min(fittable$fit),]
      
    }
    
    # save notebook for debugging and monitoring
    # saveRDS(nbook, paste("nbook",s,"_",cnt,".rds", sep=""))
    
    lrow <- nrow(nbook)
    
    finpar <- calib_pars
    finpar[1,] <- c(nbook$w_a[lrow], nbook$w_b[lrow], nbook$w_sig0[lrow])
    
    sim <- ALW_model(finpar, preset_pars, obs, burnin)
    sm <- moments_fx(sim)[1:dimmom]
    finfit <- fiterr_msm_min(tm, W, sm)*(10^mult)
    
    out <- c(s, finfit, as.numeric(finpar[1,]))
    
  } # foreach loop end
  
  res = data.frame(matrix(unlist(y), nrow = 30, ncol = 5, byrow=TRUE))
  colnames(res) = c("seed", "gof", "a", "b", "sig0")
  
  saveRDS(res, paste("res",s, ".rds", sep=""))
  
}