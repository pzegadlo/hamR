# 1: reference parameter set (following Platt, 2020)

bh_set1 <- data.frame(
                      g1 = 0,
                      b1 = 0,
                      g2 = 0.6,
                      b2 = 0.2,
                      g3 = 0.7,
                      b3 = -0.2,
                      g4 = 1.01,
                      b4 = 0,
                      r = 0.01,
                      beta = 10,
                      sigma = 0.04
                      )

# 2: decrease noise

bh_set2 <- data.frame(
                      g1 = 0,
                      b1 = 0,
                      g2 = 0.6,
                      b2 = 0.2,
                      g3 = 0.7,
                      b3 = -0.2,
                      g4 = 1.01,
                      b4 = 0,
                      r = 0.01,
                      beta = 10,
                      sigma = 0.01
                    )

# 3: decrease propensity to switch

bh_set3 <- data.frame(
                      g1 = 0,
                      b1 = 0,
                      g2 = 0.6,
                      b2 = 0.2,
                      g3 = 0.7,
                      b3 = -0.2,
                      g4 = 1.01,
                      b4 = 0,
                      r = 0.01,
                      beta = 5,
                      sigma = 0.04
                    )

# 4. Make agent 3 contrarian

bh_set4 <- data.frame(
                      g1 = 0,
                      b1 = 0,
                      g2 = 0.6,
                      b2 = 0.2,
                      g3 = -0.7,
                      b3 = -0.2,
                      g4 = 1.01,
                      b4 = 0,
                      r = 0.01,
                      beta = 10,
                      sigma = 0.04
                    )

# 5. Make agent 2 unbiased

bh_set5 <- data.frame(
                      g1 = 0,
                      b1 = 0,
                      g2 = 0.6,
                      b2 = 0.0,
                      g3 = 0.7,
                      b3 = -0.2,
                      g4 = 1.01,
                      b4 = 0,
                      r = 0.01,
                      beta = 10,
                      sigma = 0.04
                    )