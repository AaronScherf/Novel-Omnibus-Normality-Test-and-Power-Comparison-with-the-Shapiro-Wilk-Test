# simulation of given distributions and application of test
source("pgls.R")

dist_iter <- function(X,cov_inv,int_inv,g_alpha,type,iter,len,df,l1,l2,l3,l4){
  tryCatch(
    {
      beta_mat <- matrix(0,nrow=iter,ncol=len)
      D_result <-  numeric(iter)
      SW_result <- numeric(iter)
      for (i in 1:iter) {
        if(type=="exp"){
          y <- rexp(n=len, rate = 1)
        } else if (type=="uni") {
          y <- runif(n=len,min=-1,max=1)
        } else if (type=="t"){
          y <- rt(n=len, df=df) 
        } else if (type=="chi") {
          y <- rchisq(n=len, df=df)
        } else if (type=="gamma"){
          y <- rgamma(n=len, shape=df)
        } else if (type=="glambda"){
          y <- rlambda(n=len,l1=l1,l2=l2,l3=l3,l4=l4)
        } else if (type=="gpareto"){
          y <- rpareto(n=len,alpha=l1,beta=l2,delta=l3,lambda=l4)
        }
        
        beta_hat <- pgls(y=y,X=X,cov_inv=cov_inv,int_inv=int_inv)
        beta_mat[i,] <- beta_hat
        g_result[i] <- t(beta_hat) %*% cov_beta_inv %*% beta_hat
        SW_result[i] <- shapiro.test(y)$p.value <= 0.05
      }
      power <- sum(g_result >= g_alpha) / length(g_result)
      print(paste(type,df," power: ",power))
      
      SW_power <- mean(SW_result)
      print(paste(type,df," SW power: ",SW_power))
      
      result_mat <- matrix(0,nrow=1,ncol=2)
      result_mat[1,1] <- power
      result_mat[1,2] <- SW_power
      
      return(result_mat)
      
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