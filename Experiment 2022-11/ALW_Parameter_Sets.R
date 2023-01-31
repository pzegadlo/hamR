# 1: reference parameter set (following Ghonghadze and Lux, 2016)

alw_set1 <- data.frame(
                        a = 0.0003,
                        b = 0.0014,
                        sig0 = 0.03,
                        nf = 100,
                        nx = 100,
                        x0 = 0,
                        p0 = 100,
                        pf0 = 100
                      )

# 2: decrease noise

alw_set2 <- data.frame(
                        a = 0.0003,
                        b = 0.0014,
                        sig0 = 0.01,
                        nf = 100,
                        nx = 100,
                        x0 = 0,
                        p0 = 100,
                        pf0 = 100
                      )
 
# 3: increase noise
 
alw_set3 <- data.frame(
                        a = 0.0003,
                        b = 0.0014,
                        sig0 = 0.06,
                        nf = 100,
                        nx = 100,
                        x0 = 0,
                        p0 = 100,
                        pf0 = 100
                      )
   
# 4. swap a and b values
  
alw_set4 <- data.frame(
                        a = 0.0014,
                        b = 0.0003,
                        sig0 = 0.03,
                        nf = 100,
                        nx = 100,
                        x0 = 0,
                        p0 = 100,
                        pf0 = 100
                      )
  
# 5. swap a and b values, increase noise
  
alw_set5 <- data.frame(
                        a = 0.0014,
                        b = 0.0003,
                        sig0 = 0.06,
                        nf = 100,
                        nx = 100,
                        x0 = 0,
                        p0 = 100,
                        pf0 = 100
                      )