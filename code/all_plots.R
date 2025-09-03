### Plots for Reproducing Paper
source("std_res.R")

all_plots <- function(len,t_results,chi_results,g_results,
                      skew_tab,gp_results_s,kurt_tab,gp_results_k){
  
  ### Example Plots Describing std_residual process
  
  n_plots <- 10
  for (i in 1:n_plots){
    y <- rnorm(n=len,mean=0,sd=1)
    s_y <- std_res(y=y)
    plot(1:len,s_y, type = "l", pch=19, main="Empirical Process for Normal Distributions",
         col = palette(rainbow(n_plots))[i] , xlab = "i index", ylab = "s*(i)",ylim = c(-0.25,0.25))
    if (i<n_plots){
      par(new=TRUE)
    } else if (i==n_plots){
      grid()
      abline(h = 0, col = "black")
    }
  }
  
  
  for (i in 1:n_plots){
    y <- rt(n=len, df=2)
    s_y <- std_res(y=y)
    plot(1:len,s_y, type = "l", pch=19, main="Empirical Process for t-Distributions",
         col = palette(rainbow(n_plots))[i] , xlab = "i index", ylab = "s*(i)",ylim = c(-0.25,0.25))
    if (i<n_plots){
      par(new=TRUE)
    } else if (i==n_plots){
      grid()
      abline(h = 0, col = "black")
    }
  }
  
  for (i in 1:n_plots){
    y <- rchisq(n=len, df=2)
    s_y <- std_res(y=y)
    plot(1:len,s_y, type = "l", pch=19, main="Empirical Process for Chi-Square Distributions",
         col = palette(rainbow(n_plots))[i] , xlab = "i index", ylab = "s*(i)",ylim = c(-0.25,0.25))
    if (i<n_plots){
      par(new=TRUE)
    } else if (i==n_plots){
      grid()
      abline(h = 0, col = "black")
    }
  }
  
  for (i in 1:n_plots){
    y <- rgamma(n=len, shape=2)
    s_y <- std_res(y=y)
    plot(1:len,s_y, type = "l", pch=19, main="Empirical Process for Gamma Distributions",
         col = palette(rainbow(n_plots))[i] , xlab = "i index", ylab = "s*(i)",ylim = c(-0.25,0.25))
    if (i<n_plots){
      par(new=TRUE)
    } else if (i==n_plots){
      grid()
      abline(h = 0, col = "black")
    }
  }
  
  ### Rotated kernel density plot describing empirical process
  source("ex_plot.R")
  plot_iter <- 100000
  
  d_norm <- ex_plot(len=len,type="Normal",iter=plot_iter,color="green")
  par(new=TRUE)
  d_t <- ex_plot(len=len,type="t",iter=plot_iter,df=2,color="red")
  par(new=TRUE)
  d_chi <- ex_plot(len=len,type="Chi Square",iter=plot_iter,df=2,color="blue")
  par(new=TRUE)
  d_uni <- ex_plot(len=len,type="Uniform",iter=plot_iter,color='purple')
  par(new=TRUE)
  d_exp <- ex_plot(len=len,type="Exponential",iter=plot_iter,color='black')
  par(new=TRUE)
  d_gamma <- ex_plot(len=len,type="Gamma",iter=plot_iter,df=2,color='orange')
  
  grid()
  legend(x = "topright",legend=c('Normal','t-distribution, df=2','Chi-Square, df=2','Uniform','Exponential','Gamma,df=2'),
         fill=c('green','red','blue','purple','black','orange'))
  
  
  
  ### Results plots
  
  
  plot(1:10, t_results[,1], type = "b", pch = 19, main="t-distribution Power Comparison",
       col = "green", xlab = "Degrees of Freedom", ylab = "Power",ylim = c(0,1))
  par(new=TRUE)
  plot(1:10, t_results[,2], type = "b", pch = 19,
       col = "red", xlab = "Degrees of Freedom", ylab = "Power",ylim = c(0,1))
  grid()
  legend(x = "topright",legend=c('PGLS','SW'), fill=c('green','red'))
  
  
  plot(1:10, chi_results[,1], type = "b", pch = 19, main="Chi-Square Distribution Power Comparison",
       col = "blue", xlab = "Degrees of Freedom", ylab = "Power",ylim = c(0,1))
  par(new=TRUE)
  plot(1:10, chi_results[,2], type = "b", pch = 19,
       col = "red", xlab = "Degrees of Freedom", ylab = "Power",ylim = c(0,1))
  grid()
  legend(x = "topright",legend=c('PGLS','SW'), fill=c('blue','red'))
  
  plot(1:10, g_results[,1], type = "b", pch = 19, main="Gamma-Distribution Power Comparison",
       col = "orange", xlab = "Degrees of Freedom", ylab = "Power",ylim = c(0,1))
  par(new=TRUE)
  plot(1:10, g_results[,2], type = "b", pch = 19,
       col = "red", xlab = "Degrees of Freedom", ylab = "Power",ylim = c(0,1))
  grid()
  legend(x = "topright",legend=c('PGLS','SW'), fill=c('orange','red'))
  
  
  plot(skew_tab[,1], gp_results_s[,1], type = "b", pch = 19, main="GPD L-Skew Power Comparison",
       col = "purple", xlab = "L-Skew Ratio", ylab = "Power",ylim = c(0,1))
  par(new=TRUE)
  plot(skew_tab[,1], gp_results_s[,2], type = "b", pch = 19,
       col = "red", xlab = "L-Skew Ratio", ylab = "Power",ylim = c(0,1))
  grid()
  legend(x = "bottomleft",legend=c('PGLS','SW'), fill=c('purple','red'))
  
  
  plot(kurt_tab[,1], gp_results_k[,1], type = "b", pch = 19, main="GPD L-Kurtosis Power Comparison",
       col = "black", xlab = "L-Kurtosis Ratio", ylab = "Power",ylim = c(0,1))
  par(new=TRUE)
  plot(kurt_tab[,1], gp_results_k[,2], type = "b", pch = 19,
       col = "red", xlab = "L-Kurtosis Ratio", ylab = "Power",ylim = c(0,1))
  grid()
  legend(x = "topleft",legend=c('PGLS','SW'), fill=c('black','red'))
  
}



