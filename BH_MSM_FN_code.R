require(dplyr)
require(randtoolbox) # for Sobol
require(RANN) # for nearest neighbours
require(foreach)
require(doParallel)
require(doRNG)
require(parallel)

source('HAMs.R')
source('CalibrationMethods.R')

calib_pars <- data.frame(cbind(0.6, 0.2, 0.7, -0.2))
colnames(calib_pars) <- c("g2", "b2", "g3", "b3")

preset_pars <- data.frame(cbind(0, 0, 1.01, 0, 0.01, 10, 0.04))
colnames(preset_pars) <- c("g1", "b1", "g4", "b4", "r", "beta", "sigma")

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

parspace[1:2,]$g2 <- c(-1,1) # full parameter space
parspace[1:2,]$b2 <- c(-1,1) # full parameter space
parspace[1:2,]$g3 <- c(-1,1) # full parameter space
parspace[1:2,]$b3 <- c(-1,1) # full parameter space

parspace[3,] <- parspace[2,] - parspace[1,]

numoptpar <- sum(parspace[3,] > 0)

# Calculate parameter values in the current parameter space

nsob <- nsob_mult * nsob_start

parset <- seq(1, nsob, 1)

parrand <- sobol(nsob, numoptpar, scrambling = 3, init = TRUE)

parval <- cbind(matrix(0,nsob,2),parrand,matrix(0,nsob,5))

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
  tru <- BH_model(calib_pars, preset_pars, 1000, 5)
  tm = moments_fx(tru)
  W = vcovNW_fx(tru, dimmom, tm)
  
  y <- foreach(cnt = 1:30, .packages= c("dplyr", "randtoolbox", "RANN")) %dorng% { 
    
    # notebook setup to keep track of the effects of each step
    
    nbook <- data.frame(matrix(0,numrep+1,26))
    
    colnames(nbook) <- c(paste(rep(c("g2","b2","g3","b3"), each=2), rep(c("D","U"),4), sep="_"),
                         "g2","b2","g3","b3","MinErr",
                         paste(rep(c("g2","b2","g3","b3"), each=2), rep(c("wD","wU"),4), sep="_"),
                         "w_g2","w_b2","w_g3","w_b3","L2")
    rownames(nbook) <- 0:numrep
    
    
    # generate price sequences with different parameter sets
    
    fittable <- data.frame(cbind(parval[,3:6],0))
    
    for (i in 1:nsob) {
      
      sim <- BH_model(parval[i,3:6], preset_pars, obs, burnin)
      sm <- moments_fx(sim)[1:dimmom]
      fittable[i,5] <- fiterr_msm_min(tm, W, sm)*(10^mult)
      
    }
    
    colnames(fittable) <- c("g2","b2","g3","b3","fit")
    
    fittable$fit[is.na(fittable$fit)] <- 2*max(na.omit(fittable$fit)) # replace NAs in the fit function with 2 * max error
    
    # Write to notebook
    
    pos <- "0"
    
    nbook[pos,1] <- parspace$g2[1]
    nbook[pos,2] <- parspace$g2[2]
    nbook[pos,3] <- parspace$b2[1]
    nbook[pos,4] <- parspace$b2[2]
    nbook[pos,5] <- parspace$g3[1]
    nbook[pos,6] <- parspace$g3[2]
    nbook[pos,7] <- parspace$b3[1]
    nbook[pos,8] <- parspace$b3[2]
    nbook[pos,9:13] <- fittable[fittable$fit == min(fittable$fit),]
    
    # Optimisation LOOP START
    
    nsob <- nsob_start
    
    parset <- seq(1, nsob, 1)
    
    for (rnd in 1:(numrep+1)) {
      
      # Choose the windows to pick the next group of points from  
      
      ## Points in multispace
      
      mintable <- fittable %>% arrange(fit)
      mintable <- mintable[1:bnum,]
      
      nearest <- nn2(mintable[,1:4], mintable[,1:4], k = swin)$nn.idx
      
      aplot <- data.frame(matrix(ncol = 11, nrow = 0))
      
      for (pp in 1:bnum) {
        
        ttab <- mintable[nearest[pp,],]
        
        tmetrics <- c(min(ttab[,1]), max(ttab[,1]), min(ttab[,2]), max(ttab[,2]), min(ttab[,3]), max(ttab[,3]), min(ttab[,4]), max(ttab[,4]), mean(ttab[,5]), se(ttab[,5]), min(ttab[,5]))
        
        aplot <- rbind(aplot, tmetrics)
        
      }
      
      colnames(aplot) <- c(paste(rep(c("g2","b2","g3","b3"), each=2), rep(c("D","U"),4), sep="_"),"mu","serr","minim")
      
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
      
      nbook[pos,14:15] <- wplot[1,1:2]
      nbook[pos,16:17] <- wplot[1,3:4]
      nbook[pos,18:19] <- wplot[1,5:6]
      nbook[pos,20:21] <- wplot[1,7:8]
      nbook[pos,22] <- sum(wplot[1,1:2])/2
      nbook[pos,23] <- sum(wplot[1,3:4])/2
      nbook[pos,24] <- sum(wplot[1,5:6])/2
      nbook[pos,25] <- sum(wplot[1,7:8])/2
      nbook[pos,26] <- sqrt(sum((nbook[pos,22:25] - trupar[1,3:6])^2))
      
      # Show some current progress stats
      
      # Finish the algorithm when the computational budget is taken, but we have the last window
      
      if (rnd > numrep) { break }
      
      # Get the new points from each window
      
      # Points from POI - EI - LCB
      
      parrand <- sobol(nsob, numoptpar, scrambling = 3, init = FALSE)
      
      parval <- cbind(matrix(0,nsob,2),parrand,matrix(0,nsob,5))
      
      inc <- 0
      
      # choice of acq method
      
      aplot <- wplot
      
      for (candwin in 1:bwin) {
        
        parspace[1:2,]$g2 <- t(aplot[candwin,1:2])
        parspace[1:2,]$b2 <- t(aplot[candwin,3:4])
        parspace[1:2,]$g3 <- t(aplot[candwin,5:6])
        parspace[1:2,]$b3 <- t(aplot[candwin,7:8])
        
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
      
      newtable <- data.frame(cbind(parval[,3:6],0))
      
      for (i in 1:nsob) {
        
        sim <- BH_model(parval[i,3:6], preset_pars, obs, burnin)
        sm <- moments_fx(sim)[1:dimmom]
        
        newtable[i,5] <- fiterr_msm_min(tm, W, sm)*(10^mult)
        
      }
      
      colnames(newtable) <- c("g2","b2","g3","b3","fit")
      
      newtable$fit[is.na(newtable$fit)] <- 2*max(na.omit(fittable$fit)) # replace NAs in the fit function with 2 * max error
      
      fittable <- rbind(fittable,newtable)
      
      # Write to notebook
      
      pos <- as.character(rnd)
      
      nbook[pos,1] <- min(parval$g2)
      nbook[pos,2] <- max(parval$g2)
      nbook[pos,3] <- min(parval$b2)
      nbook[pos,4] <- max(parval$b2)
      nbook[pos,5] <- min(parval$g3)
      nbook[pos,6] <- max(parval$g3)
      nbook[pos,7] <- min(parval$b3)
      nbook[pos,8] <- max(parval$b3)
      nbook[pos,9:13] <- fittable[fittable$fit == min(fittable$fit),]
      
    }
    
    # save notebook for debugging and monitoring
    # saveRDS(nbook, paste('nbook',s,'_',cnt,'.rds', sep=""))
    
    lrow <- nrow(nbook)
    
    finpar <- calib_pars
    finpar[1,] <- c(nbook$w_g2[lrow], nbook$w_b2[lrow], nbook$w_g3[lrow], nbook$w_b3[lrow])
    
    sim <- BH_model(finpar, preset_pars, obs, burnin)
    sm <- moments_fx(sim)[1:dimmom]
    finfit <- fiterr_msm_min(tm, W, sm)*(10^mult)
    
    out <- c(s, finfit, as.numeric(finpar[1,]))
    
  } # foreach loop end
  
  res = data.frame(matrix(unlist(y), nrow = 30, ncol = 6, byrow=TRUE))
  colnames(res) = c('seed', 'gof', 'g2', 'b2', 'g3', 'b3')
  
  saveRDS(res, paste('res',s, '.rds', sep=""))
  
}