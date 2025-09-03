#### Penalized Generalized Least Squares function

# Import the std_residual function
source("std_res.R")
source("calc_sim_mat.R")

pgls <- function(y,X,cov_inv,int_inv){
  tryCatch(
    {
      r <- pnorm((sort(y) - mean(y))/sd(y)) - (1:20)/21
      beta_hat <- int_inv %*% t(X) %*% cov_inv %*% r
      return(beta_hat)
      
    },
    error=function(e) {
      message('An error ocurred')
      print(e)
    },
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
}