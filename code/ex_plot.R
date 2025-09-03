### Example Plots describing s(i) for various distributions
source("std_res.R")

ex_plot <- function(len,type,iter,df,l1,l2,l3,l4,color){
  tryCatch(
    {
      plot_mat <- matrix(0,nrow=iter,ncol=len)
    
      for (i in 1:iter) {
        if(type=="Exponential"){
          y <- rexp(n=len, rate = 1)
        } else if (type=="Uniform") {
          y <- runif(n=len,min=-1,max=1)
        } else if (type=="t"){
          y <- rt(n=len, df=df) 
        } else if (type=="Chi Square") {
          y <- rchisq(n=len, df=df)
        } else if (type=="Gamma"){
          y <- rgamma(n=len, shape=df)
        } else if (type=="Generalized Lambda"){
          y <- rlambda(n=len,l1=l1,l2=l2,l3=l3,l4=l4)
        } else if (type=="Generalized Pareto"){
          y <- rpareto(n=len,alpha=l1,beta=l2,delta=l3,lambda=l4)
        } else if (type=="Normal"){
          y <- rnorm(n=len,mean=0,sd=1)
        }
        
        plot_mat[i,] <- std_res(y=y)
      }
      plot_ex <- colMeans(plot_mat)
        
      #h<-hist(plot_ex, breaks=len, xlab="s(i) value",
       #         main=paste("Empirical Process for",type,"Distribution"))
      
      ex_dens <- density(plot_ex,na.rm = T) # returns the density data
      
      plot(ex_dens, type = "b", pch=19,cex=0.5, main="Comparing Distributions after Empirical Process",
           col = color, xlab = "Expected s*(i) Value", ylab = "Density",xlim=c(-0.07,0.07),ylim = c(0,40))
      
      
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