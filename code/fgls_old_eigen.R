#### Test Feasible General Least Squares ####

# Import the std_residual function
source("calc_sim_mat.R")


len <- 20
sim_iter <- 1000000
y_mean <- 0
y_sd <- 1

y <- rnorm(n=len,mean=y_mean,sd=y_sd)

sim_mat <- calc_sim_mat(len=length(y),mean=mean(y),sd=sd(y),iter=sim_iter)
cov_mat <- cov(sim_mat)
cov_inv <- solve(cov_mat)


lambda <- 6


# differences matrix with standard deviations per column vector of sim_mat, through n-1
D_sigma <- cbind(diag(19),numeric(19)) # start with a 19x19 identity, then append a column of 0's to make 19x20
for (row in 1:dim(D_sigma)[1]){ 
  D_sigma[row,row+1] <- -1 # for every row, add a -1 in the column to the right of the diagonal
  var_av <- sqrt((cov_mat[row,row] + cov_mat[row+1,row+1])/2)
  D_sigma[row,] <- D_sigma[row,] * var_av # multiply the row by the moving variance average from the covariance matrix
}  # results in 19x20 rectangular matrix

# penalty matrix of D_sig^T * D_sig; results in 20x20 square matrix
P_mat <- t(D_sigma) %*% D_sigma

# design matrix of function phi_i, determined by eigenvectors of penalty matrix P, sorted from smallest corresponding eigenvalue to largest
eigen_P <- eigen(P_mat)
X <- eigen_P$vectors
eigen_matrix <- diag(length(eigen_P$values)) * eigen_P$values

interior <- t(X) %*% cov_inv %*% X 
interior_lambda <- interior + (lambda * eigen_matrix)
int_inv <- solve(interior_lambda)

decomp_R <- chol(cov_mat)
decomp_R_inv <- solve(decomp_R)
decomp_R_t_inv <- solve(t(decomp_R))
P <- decomp_R_t_inv %*% X %*% int_inv %*% t(X) %*% decomp_R_inv
print(sum(diag(P)))

#### Normal Iterations Test for Expected Value ####

source("fgls.R")

iter <- 10000

beta_mat <- matrix(0,nrow=iter,ncol=len)

for (i in 1:iter) {
  y <- rnorm(n=len,mean=y_mean,sd=y_sd)
  beta_mat[i,] <- fgls(y=y,X=X,cov_inv=cov_inv,int_inv=int_inv)
} 

cov_beta <- cov(beta_mat)
cov_beta_inv <- solve(cov_beta)

#### Calculating 95th percentile of distance metric by simulation for normal input vector ####

D_result <-  numeric(iter)

for (i in 1:iter) {
  y <- rnorm(n=len,mean=y_mean,sd=y_sd)
  beta_hat <- fgls(y=y,X=X,cov_inv=cov_inv,int_inv=int_inv)
  D_result[i] <- t(beta_hat) %*% cov_beta_inv %*% beta_hat
  
}

D_alpha <- quantile(D_result,probs = 0.95)


source("dist_iter.R")

#### Exponential Iterations Test for Expected Value ####

dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,D_alpha=D_alpha,type="exp",iter=iter,len=len,df=1)

#### Uniform Iterations Test for Expected Value ####

dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,D_alpha=D_alpha,type="uni",iter=iter,len=len,df=1)

df_index <- 1:10
#### t-distribution Iterations Test for Expected Value ####
t_results <- {}
for (i in df_index) {
  t_results[i] <- dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,D_alpha=D_alpha,
                            type="t",iter=iter,len=len,df=i)
}


#### Chi-square Iterations Test for Expected Value ####
chi_results <- {}
for (i in 1:10) {
  chi_results[i] <- dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,D_alpha=D_alpha,
                              type="chi",iter=iter,len=len,df=i)
}

#### Gamma Iterations Test for Expected Value ####
g_results <- {}
for (i in 1:10) {
  g_results[i] <- dist_iter(X=X,cov_inv=cov_inv,int_inv=int_inv,D_alpha=D_alpha,
                            type="gamma",iter=iter,len=len,df=i)
}


plot(1:10, t_results, type = "b", pch = 19, main="Old Power Comparison by Disribution",
     col = "green", xlab = "degrees of freedom", ylab = "power",ylim = c(0,1))
par(new=TRUE)
plot(1:10,chi_results, type = "b", pch = 19,col = "red",ylab="",xlab="",axes=FALSE,ylim = c(0,1))
par(new=TRUE)
plot(1:10,g_results, type = "b", pch = 19,col = "blue",ylab="",xlab="",axes=FALSE,ylim = c(0,1))
grid()
legend(x = "topright",legend=c('t','chi', 'gamma'), fill=c('green','red', 'blue'))

