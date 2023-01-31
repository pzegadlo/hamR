### Alfarano Lux Wagner model

ALW_model <- function(calib_pars, preset_pars, obs, burnin, 
                      rescale = data.frame(
                      a = 1,
                      b = 1,
                      sig0 = 1,
                      nf = 1,
                      nx = 1,
                      x0 = 1,
                      p0 = 1,
                      pf0 = 1
                    ))

{
  
  if(is.data.frame(calib_pars) == FALSE) { stop('Argument calib_pars has to be a dataframe.')}
  if(is.data.frame(preset_pars) == FALSE) { stop('Argument preset_pars has to be a dataframe.')}
  
  pars <- cbind(calib_pars, preset_pars)
  
  if(setequal(colnames(pars), c('a', 'b', 'sig0', 'nf', 'nx', 'x0', 'p0', 'pf0')) == FALSE | length(colnames(pars)) != 8) {
    
    stop('Incorrect parameter names.') 
    
    }
  
  # rescaling maybe needed for some optimisation algorithms to work correctly
  
  a <- pars[['a']] * rescale[['a']]
  b <- pars[['b']] * rescale[['b']]
  sig0 <- pars[['sig0']] * rescale[['sig0']]
  nf <- pars[['nf']] * rescale[['nf']]
  nx <- pars[['nx']] * rescale[['nx']]
  x0 <- pars[['x0']] * rescale[['x0']]
  p0 <- pars[['p0']] * rescale[['p0']]
  pf0 <- pars[['pf0']] * rescale[['pf0']]
  
  out <- data.frame(matrix(0, ncol = 4, nrow = obs))
  colnames(out) <- c('x', 'pf', 'p', 'r')
  
  Nser = 0.5*((x0)*nx+nx) # number of optimist chartists
  
  #	Loop 1: Macro-Time and Micro-Time
  
  x = x0
  p = p0
  pf = pf0
  trun = 0.0
  jj = 1
  
  while ((trun < obs + burnin) || (jj <= obs + burnin)) {
    
    #	Loop: Agents decide whether to switch
    
    jjcomp = jj
    
    omeg1 = a + b*nx*(1 + x)/2
    omeg2 = a + b*nx*(1 - x)/2
    
    if (trun < jjcomp) {
      
      v1 = runif(1)
      lamda = 0.5*(1+x)*nx*omeg2+0.5*(1-x)*nx*omeg1
      tchange = -log(v1)/lamda
      trun = trun + tchange
      
      v2 = runif(1)
      if (v2 <= 0.5*(1+x)*nx*omeg2/lamda) { 
        Nser = Nser - 1 
      } else {
        Nser = Nser +1
      }
      xold = x
      x = 2*Nser/nx -1 # disposition of chartists (-1: fully pessimistic, +1: fully optimistic)
      
    }
    
    if (trun >= jjcomp) {
      
      pold = p
      
      r1 = rnorm(1) # standard normal random number
      pf = pf*exp(sig0*r1) # random evolution of the fundamental price
      p = pf*exp((nx/nf)*xold) # market price update equation
      
      if (jj > burnin) {
        
        out[jj-burnin,1] = xold
        out[jj-burnin,2] = pf	
        out[jj-burnin,3] = p
        out[jj-burnin,4] = log(p) - log(pold)
        
        
      }
      
      jj = jj + 1
      
    }
    
  }
  
  return(out$r)
  
}

### Brock Hommes model

BH_model <- function(calib_pars, preset_pars, obs, burnin) {
  
  if(is.data.frame(calib_pars) == FALSE) { stop('Argument calib_pars has to be a dataframe.')}
  if(is.data.frame(preset_pars) == FALSE) { stop('Argument preset_pars has to be a dataframe.')}
  
  pars <- cbind(calib_pars, preset_pars)
  
  if(setequal(colnames(pars), c("g1", "b1", "g2", "b2", "g3", "b3", "g4", "b4", "r", "beta", "sigma")) == FALSE) { stop('Incorrect parameter names.') }
  
  x <- matrix(0, obs+burnin, 1)
  U <- matrix(0, obs+burnin, 4)
  n <- matrix(0, obs+burnin, 4)
  
  R <- 1 + pars[['r']]
  beta <- pars[['beta']]
  sigma <- pars[['sigma']]
  
  g <- pars[1, grepl("^g" , names(pars))]
  g <- g[ , order(names(g))]
  b <- pars[1, grepl("^b[0-9]" , names(pars))]
  b <- b[ , order(names(b))]
  
  e <- rnorm(obs+burnin, 0, sigma)
  
  x[1:3] <- 0.01
  
  # HAM loop
  
  for (t in 3:(obs+burnin-1)) {
    
    U[t,] <- (x[t]-R*x[t-1])*(as.matrix(g)*x[t-2]+as.matrix(b) - R*x[t-1])
    
    n[t+1,] <- exp(beta*U[t,]) / sum(exp(beta*U[t,]))
    
    x[t+1] <- (1/R) * (sum(n[t+1,]*(g*x[t]+b)) + e[t+1])
    
  }
  
  output <- x[(burnin+1):(obs+burnin)]
  
}
