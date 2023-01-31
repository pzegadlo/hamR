### MLE - fit error calculation

  # Gaussian Kernel
  
  gk <- function(x, h) {
    
    return (exp( (-x^2) / (2*h^2) ))
    
  }
  
  # modified Gaussian kernel
  
  mgk <- function(x, h) {
    
    return ( gk(x/h, h) / h )
    
  }
  
  # Fit error calculation
  # h argument is optional - it can be calculated once and supplied to the function as a given,
  # which speeds up the calculations
  
  cond_dens <- function(sim, tru, h) {
    
    if(missing(h)) { h <- ((4/3)^(1/5)) * sd(tru) }
    
    return( mgk(sim-tru, h) ) 
    
  }
  
  # Version to minimize target function
  
  fiterr_mle_min <- function(cd, bestlik) {
    
    return( -sum(log(cd)) - bestlik )
    
  }
  
  # Version to maximize target function
  
  fiterr_mle_max <- function(cd, bestlik) {
    
    return( bestlik + sum(log(cd)) )
    
  }
  

### MSM fit error calculation

  # Calculating 7  moments in the dataset
  
  moments_fx <- function(tru) {
    
    obs <- length(tru)
    
    sqtru <- tru^2
    abtru <- abs(tru)
    
    tm2 <- mean(tru^2)  # variance
    tm4 <- mean(tru^4)  # kurtosis is the 4th power of the standardised variable values, so this is directly proportional
    tma1 <- mean(tru)^2 - (mean(tru)/obs)*(tru[1]+tru[obs]+2*sum(tru[2:(obs-1)])) + sum(tru[2:(obs-1)]*tru[3:obs])/obs                # first-order autocovariance of raw returns
    tma2 <- mean(sqtru)^2 - (mean(sqtru)/obs)*(sqtru[1]+sqtru[obs]+2*sum(sqtru[2:(obs-1)])) + sum(sqtru[2:(obs-1)]*sqtru[3:obs])/obs  # first-order autocovariance of squared returns
    tmaa <- mean(abtru)^2 - (mean(abtru)/obs)*(abtru[1]+abtru[obs]+2*sum(abtru[2:(obs-1)])) + sum(abtru[2:(obs-1)]*abtru[3:obs])/obs  # first-order autocovariance of absolute returns
    tma25 <- mean(sqtru)^2 - (mean(sqtru)/obs)*(sum(sqtru[1:5])+sum(sqtru[(obs-4):obs])+2*sum(sqtru[6:(obs-5)])) + sum(sqtru[1:(obs-5)]*sqtru[6:obs])/obs  # 5-order autocovariance of squared returns
    tmaa5 <- mean(abtru)^2 - (mean(abtru)/obs)*(sum(abtru[1:5])+sum(abtru[(obs-4):obs])+2*sum(abtru[6:(obs-5)])) + sum(abtru[1:(obs-5)]*abtru[6:obs])/obs  # 5-order autocovariance of absolute returns
    
    tm <- c(tm2,tm4,tma1,tma2,tmaa,tma25,tmaa5)
    
    return(tm)
    
  }
  
  # Finding the Newey-West estimator of the vcov matrix of the empirical moments
  # Potentially need to take into account burnin too (correcting obs)
  
  vcovNW_fx <- function (tru, dimmom, tm) {
    
    obs <- length(tru)
    
    sqtru <- tru^2
    abtru <- abs(tru)
    
    p <- ceiling(obs^(1/4)) # rule of thumb for lag order choice, Franke 2008
    
    Gamma <- array(dim=c(p+1,dimmom,dimmom))
    
    for (j in 0:p) {        # calculating Gammas for the estimator
      
      # Loop not needed for these moments, but acov will be calculated in the loop
      
      m2t  <- tru[(j+1):obs]^2
      m2tj <- tru[1:(obs-j)]^2
      
      m4t  <- tru[(j+1):obs]^4
      m4tj <- tru[1:(obs-j)]^4
      
      ma1t  <- matrix(0,length(m2t),1)
      ma1tj <- matrix(0,length(m2t),1)        
      
      ma2t  <- matrix(0,length(m2t),1)
      ma2tj <- matrix(0,length(m2t),1)      
      
      maat  <- matrix(0,length(m2t),1)
      maatj <- matrix(0,length(m2t),1)  
      
      ma25t  <- matrix(0,length(m2t),1)
      ma25tj <- matrix(0,length(m2t),1)      
      
      maa5t  <- matrix(0,length(m2t),1)
      maa5tj <- matrix(0,length(m2t),1)         
      
      gam <- matrix(0,dimmom,dimmom)
      
      # potentially from j+1 to get more accurate gammas 
      for (t in (j+6):obs) {
        
        # single / momentary:
        
        ma1t[t-j]  <- (tru[t]-mean(tru[(j+1):t]))  * (tru[t-1]-mean(tru[(j+1):t]))
        ma1tj[t-j] <- (tru[t-j]-mean(tru[1:(t-j)])) * (tru[t-j-1]-mean(tru[1:(t-j)]))          
        
        ma2t[t-j]  <- (sqtru[t]-mean(sqtru[(j+1):t]))   * (sqtru[t-1]-mean(sqtru[(j+1):t]))
        ma2tj[t-j] <- (sqtru[t-j]-mean(sqtru[1:(t-j)])) * (sqtru[t-j-1]-mean(sqtru[1:(t-j)]))      
        
        maat[t-j]  <- (abtru[t]-mean(abtru[(j+1):t]))   * (abtru[t-1]-mean(abtru[(j+1):t]))
        maatj[t-j] <- (abtru[t-j]-mean(abtru[1:(t-j)])) * (abtru[t-j-1]-mean(abtru[1:(t-j)]))      
        
        ma25t[t-j]  <- (sqtru[t]-mean(sqtru[(j+5):t]))   * (sqtru[t-5]-mean(sqtru[(j+5):t]))
        ma25tj[t-j] <- (sqtru[t-j]-mean(sqtru[1:(t-j)])) * (sqtru[t-j-5]-mean(sqtru[1:(t-j)]))      
        
        maa5t[t-j]  <- (abtru[t]-mean(abtru[(j+5):t]))   * (abtru[t-5]-mean(abtru[(j+5):t]))
        maa5tj[t-j] <- (abtru[t-j]-mean(abtru[1:(t-j)])) * (abtru[t-j-5]-mean(abtru[1:(t-j)]))  
        
        mt <- c(m2t[t-j-1],m4t[t-j-1],ma1t[t-j-1],ma2t[t-j-1],maat[t-j-1],ma25t[t-j-1],maa5t[t-j-1]) 
        mtj <- c(m2tj[t-j-1],m4tj[t-j-1],ma1tj[t-j-1],ma2tj[t-j-1],maatj[t-j-1],ma25tj[t-j-1],maa5tj[t-j-1])
        
        mt <- mt[1:dimmom]
        mtj <- mtj[1:dimmom]
        
        gam <- gam + (mt-tm) %*% t(mtj-tm)
        
      }
      
      Gamma[j+1,,] <- gam
      
    }
    
    omegapart <- matrix(0,dimmom,dimmom)
    
    for (j in 1:p) { 
      
      omegapart <- omegapart + (1-j/(p+1))*(Gamma[j+1,,] + t(Gamma[j+1,,]))
      
    }
    
    omega <- Gamma[1,,] + omegapart
    
    W <- solve(omega)
    
    return(W)
    
  }
  
  
  # Version to minimize target function
  fiterr_msm_min <- function(tm, W, tmsim) { t(tmsim-tm) %*% W %*% (tmsim-tm) }
  
  # Version to maximize target function
  fiterr_msm_max <- function(tm, W, tmsim) { -1 * (t(tmsim-tm) %*% W %*% (tmsim-tm)) }
  
  
### Helper functions
  
  se <- function(x) {
    
    return(sd(x) / sqrt(length(x)))
    
  }