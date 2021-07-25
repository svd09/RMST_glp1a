#######################################################################
##  WALD TEST 2 TAILED p-VALVE to compare sub group results from a   ##
##                          meta-analysis.                           ##
#######################################################################

# 2021/07/24
# Salil Deo 
# Created as part of paper on RMST for Glp1A drugs.
# create a function so that can run that for the other values.
# this function can be used for the other tests.

wald_test_MA <- function(Estimate, SE, Estimate2, SE2){

  stopifnot(is.numeric(Estimate))
  stopifnot(is.numeric(SE))
  stopifnot(is.numeric(Estimate2))
  stopifnot(is.numeric(SE2))
 
  a = (Estimate*Estimate) - (Estimate2*Estimate2)
  
  b = sqrt((SE*SE) + (SE2*SE2))
  
  z = a/b
  
  p = 2*pnorm(-abs(z))
  
  return(p)
}

# example
# wald_test_MA(1.2, 3, 1, 2)

# example of error
# wald_test_MA(1.3, 0.4, w, 0.03)


