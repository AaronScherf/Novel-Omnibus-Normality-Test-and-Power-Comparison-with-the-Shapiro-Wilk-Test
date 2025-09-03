### Generalized Pareto Distribution preparation calculations

g_pareto_prep <- function(g_alpha,g_beta,g_delta,g_lambda,type){
# delta skewness parameter -1 to 1; lambda held at 

  if(type=="skew"){
    # n=len,alpha=l1,beta=l2,delta=l3,lambda=l4
    skew_tab <- matrix(0,nrow=14,ncol=5)
    skew_tab[1,] <- c(0,   g_alpha, g_beta,  -0.1,   g_lambda)
    skew_tab[2,] <- c(0,g_alpha, g_beta, 0, g_lambda)
    skew_tab[3,] <- c(0, g_alpha, g_beta, 0.1,  g_lambda)
    skew_tab[4,] <- c(0,g_alpha, g_beta,   0.2,  g_lambda)
    skew_tab[5,] <- c(0,    g_alpha, g_beta,   0.3,   g_lambda)
    skew_tab[6,] <- c(0, g_alpha, g_beta,  0.4, g_lambda)
    skew_tab[7,] <- c(0,  g_alpha, g_beta,  0.5,  g_lambda) # normal distribution
    skew_tab[8,] <- c(0, g_alpha, g_beta,  0.6, g_lambda)
    skew_tab[9,] <- c(0,    g_alpha, g_beta,   0.7,   g_lambda)
    skew_tab[10,] <- c(0,    g_alpha, g_beta,   0.8,   g_lambda)
    skew_tab[11,] <- c(0,    g_alpha, g_beta,   0.9,   g_lambda) 
    skew_tab[12,] <- c(0,    g_alpha, g_beta,   1,   g_lambda)
    skew_tab[13,] <- c(0,    g_alpha, g_beta,   1.1,   g_lambda)
    skew_tab[14,] <- c(0,    g_alpha, g_beta,   1.2,   g_lambda)
    
    # scaling the beta and alpha parameters to fix L1 at 0 and L2 at 1
    for (i in 1:dim(skew_tab)[1]){
      skew_tab[i,3] <- ((skew_tab[i,5]+1)*(skew_tab[i,5]+2))
      skew_tab[i,2] <- (skew_tab[i,3]*(1-2*skew_tab[i,4])/(skew_tab[i,5]+1))
    }
    
    l_moments <- matrix(0,nrow=dim(skew_tab)[1],ncol=4)
    for (i in 1:dim(skew_tab)[1]){
      print(i)
      a<-skew_tab[i,2]
      b<-skew_tab[i,3]
      d<-skew_tab[i,4]
      l<-skew_tab[i,5]
      l_moments[i,1] <- a - (b*(1-2*d)/(l+1))
      l_moments[i,2] <- b / ((l+1)*(l+2))
      l_moments[i,3] <- (((l-1)*(1-2*d))/(l+3))
      l_moments[i,4] <- (((l-1)*(l-2))/((l+3)*(l+4)))
      skew_tab[i,1] <- l_moments[i,3]
    }
    print("L-Moments Ratio, Skew")
    
    print(l_moments)
    tab <- skew_tab
    
  } else if (type=="kurt"){
    # lambda kurt parameter -1 to 1; delta held at 0.5
    kurt_tab <- matrix(0,nrow=9,ncol=5)
    kurt_tab[1,] <- c(-0.9,  g_alpha, g_beta,g_delta, -0.9)
    kurt_tab[2,] <- c(-0.7,  g_alpha, g_beta,g_delta, -0.7)
    kurt_tab[3,] <- c(-0.5,  g_alpha, g_beta,g_delta, -0.5)
    kurt_tab[4,] <- c(-0.3,  g_alpha, g_beta,g_delta, -0.3)
    kurt_tab[5,] <- c(0,     g_alpha, g_beta,g_delta, 0.01)
    kurt_tab[6,] <- c(0.3,   g_alpha, g_beta,g_delta, 0.3)
    kurt_tab[7,] <- c(0.5,   g_alpha, g_beta,g_delta, 0.5)
    kurt_tab[8,] <- c(0.7,   g_alpha, g_beta,g_delta, 0.7)
    kurt_tab[9,] <- c(0.9,   g_alpha, g_beta,g_delta, 0.9)
    
    for (i in 1:dim(kurt_tab)[1]){
      kurt_tab[i,3] <- ((kurt_tab[i,5]+1)*(kurt_tab[i,5]+2))
    }
    
    l_moments <- matrix(0,nrow=dim(kurt_tab)[1],ncol=4)
    for (i in 1:dim(kurt_tab)[1]){
      print(i)
      a<-kurt_tab[i,2]
      l<-kurt_tab[i,5]
      b<-kurt_tab[i,3]
      d<-kurt_tab[i,4]
      
      l_moments[i,1] <- a - (b*(1-2*d)/(l+1))
      l_moments[i,2] <- b / ((l+1)*(l+2))
      l_moments[i,3] <- (((l-1)*(1-2*d))/(l+3))
      l_moments[i,4] <- (((l-1)*(l-2))/((l+3)*(l+4)))
      kurt_tab[i,1] <- l_moments[i,4]
    }
    print("L-Moments Ratio, Kurtosis")
    print(l_moments)
    tab <- kurt_tab
  }
  
  return(tab)
}
