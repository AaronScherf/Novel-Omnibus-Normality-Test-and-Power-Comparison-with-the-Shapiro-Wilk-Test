#### Main Penalized General Least Squares Test ####

len <- 20
sim_iter <- 1000000
y_mean <- 0
y_sd <- 1
lambda <- 20

# application of standardization the probability integral transform

source("calc_sim_mat.R")

y <- rnorm(n=len,mean=y_mean,sd=y_sd)

# calculation of the empirical process via simulation

sim_mat <- calc_sim_mat(len=length(y),mean=mean(y),sd=sd(y),iter=sim_iter)

# covariance matrix of the standardized probability integral transform

cov_mat <- cov(sim_mat)
cov_inv <- solve(cov_mat)

# differences matrix with standard deviations per column vector of sim_mat, through n-1

D_sigma <- cbind(diag(19),numeric(19)) 
for (row in 1:dim(D_sigma)[1]){ 
  D_sigma[row,row+1] <- -1 
  var_av <- sqrt((cov_mat[row,row] + cov_mat[row+1,row+1])/2)
  D_sigma[row,] <- D_sigma[row,] * var_av 
}  

# penalty matrix of D_sig^T * D_sig; results in 20x20 square matrix

# calculation of eigenvectors and design matrix X

eigen_cov <- eigen(cov_mat)

X <- eigen_cov$vectors

# calculation of P matrix

P_mat <- t(X) %*% t(D_sigma) %*% D_sigma %*% X

interior <- t(X) %*% cov_inv %*% X 

interior_lambda <- interior + (lambda * P_mat)
int_inv <- solve(interior_lambda)

# Approximating the effective degrees of freedom via trace of H

H <- X %*% int_inv %*% t(X) %*% cov_inv
print(sum(diag(H)))

#### Normal Iterations Test for Expected Value ####

source("pgls.R")

iter <- 10000

beta_mat <- matrix(0,nrow=iter,ncol=len)

for (i in 1:iter) {
  y <- rnorm(n=len,mean=y_mean,sd=y_sd)
  beta_mat[i,] <- pgls(y=y,X=X,cov_inv=cov_inv,int_inv=int_inv)
} 

cov_beta <- cov(beta_mat) #C
cov_beta_inv <- solve(cov_beta) #C^-1

#### Calculating 95th percentile of distance metric by simulation for normal input vector ####

g_result <-  numeric(iter)

for (i in 1:iter) {
  y <- rnorm(n=len,mean=y_mean,sd=y_sd)
  beta_hat <- pgls(y=y,X=X,cov_inv=cov_inv,int_inv=int_inv)
  g_result[i] <- t(beta_hat) %*% cov_beta_inv %*% beta_hat
  
}

g_alpha <- quantile(g_result,probs = 0.95)

source("dist_iter.R")
source("rlambda.R")

#### Exponential Iterations Test for Expected Value ####

dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,g_alpha=g_alpha,type="exp",iter=iter,len=len,df=1)

#### Uniform Iterations Test for Expected Value ####

dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,g_alpha=g_alpha,type="uni",iter=iter,len=len,df=1)


#### t-distribution Iterations Test for Expected Value ####
df_index <- 1:10

t_results <- matrix(0,nrow=10,ncol=2)
for (i in df_index) {
  print(i)
  t_results[i,] <- dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,g_alpha=g_alpha,
                           type="t",iter=iter,len=len,df=i)
}

#### Chi-square Iterations Test for Expected Value ####
chi_results <- matrix(0,nrow=10,ncol=2)
for (i in 1:10) {
  chi_results[i,] <- dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,g_alpha=g_alpha,
                              type="chi",iter=iter,len=len,df=i)
}

#### Gamma Iterations Test for Expected Value ####
g_results <- matrix(0,nrow=10,ncol=2)
for (i in 1:10) {
  g_results[i,] <- dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,g_alpha=g_alpha,
                            type="gamma",iter=iter,len=len,df=i)
}

### Generalized Pareto Distribution ###
source("rpareto.R")
source("g_pareto_prep.R")

skew_tab <- g_pareto_prep(g_alpha = 0,g_beta = 1,g_delta=0.5,g_lambda = 4.256,type="skew")

kurt_tab <- g_pareto_prep(g_alpha = 0,g_beta = 1,g_delta=0.5,g_lambda = 4.256,type="kurt")

gp_results_s <- matrix(0,nrow=dim(skew_tab)[1],ncol=2)
for (i in 1:dim(skew_tab)[1]) {
  gp_results_s[i,] <- dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,g_alpha=g_alpha,
                                type="gpareto",iter=iter,len=len,df=skew_tab[i,1],
                                l1=skew_tab[i,2],l2=skew_tab[i,3],l3=skew_tab[i,4],l4=skew_tab[i,5])
}

gp_results_k <- matrix(0,nrow=dim(kurt_tab)[1],ncol=2)
for (i in 1:dim(kurt_tab)[1]) {
  gp_results_k[i,] <- dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,g_alpha=g_alpha,
                                type="gpareto",iter=iter,len=len,df=kurt_tab[i,1],
                                l1=kurt_tab[i,2],l2=kurt_tab[i,3],l3=kurt_tab[i,4],l4=kurt_tab[i,5])
}


# Plotting all results 

all_plots(len=len,t_results=t_results,chi_results=chi_results,g_results=g_results,
                      skew_tab=skew_tab,gp_results_s=gp_results_s,
                      kurt_tab=kurt_tab,gp_results_k=gp_results_k)
  