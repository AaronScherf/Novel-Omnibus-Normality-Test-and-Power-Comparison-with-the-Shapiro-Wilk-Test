#### Residual Simulation Matrix Function 

source("std_res.R")

calc_sim_mat <- function(len,mean,sd,iter){
  tryCatch(
    {
      sim_mat <- matrix(0,nrow=iter,ncol=len)
      for (i in 1:iter) {
        y <- rnorm(n=len,mean=mean,sd=sd)
        sim_mat[i,] <- std_res(y=y) # each sample vector is stored as a row
      } 
      
      return(sim_mat)
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

