### Test Script for Residual Function ###

# Import the std_residual function
source("std_residual.R")

### Test for Error Handling ###

y <- 'bad_input'

norm_result <- std_res(y=y)

print(mean(norm_result))


### Single Test for Normal Distribution Sample ###

# Defining the sample length, mean, and standard deviation
# Theoretically, these shouldn't affect the outcome in a repeated test

sample_len <- 20
sample_mean <- 0
sample_sd <- 1

y <- rnorm(n=sample_len,mean=sample_mean,sd=sample_sd)

norm_result <- std_res(y=y)

print(mean(norm_result))
print(sd(norm_result))


### Repeated Test for Normal Distribution ###

iter <- 100000

sim_mat <- matrix(0,nrow=iter,ncol=sample_len)

for (i in 1:iter) {
  y <- rnorm(n=sample_len,mean=sample_mean,sd=sample_sd)
  sim_mat[i,] <- std_res(y=y)
} 

res_avg <- colMeans(sim_mat)
print(res_avg)
plot(res_avg)
title("Normal Distribution Residual Plot")
print(mean(res_avg))

### Single Test for Uniform Distribution Sample ###
sample_min <- -2
sample_max <- 2

y <- runif(n=sample_len,min=sample_min,max=sample_max)

unif_result <- std_res(y=y)

print(mean(unif_result))
print(sd(unif_result))


### Repeated Test for Uniform Distribution ###

sim_mat <- matrix(0,nrow=iter,ncol=sample_len)

for (i in 1:iter) {
  y <- runif(n=sample_len,min=sample_min,max=sample_max)
  sim_mat[i,] <- std_res(y=y)
} 

res_avg <- colMeans(sim_mat)
print(res_avg)
plot(res_avg)
title("Uniform Distribution Residual Plot")
print(mean(res_avg))

### Single Test for Exponential Distribution ###
rate <- 1
y <- rexp(n=sample_len, rate = rate) 

exp_result <- std_res(y=y)

print(mean(exp_result))
print(sd(exp_result))


### Repeated Test for Exponential Distribution ###

sim_mat <- matrix(0,nrow=iter,ncol=sample_len)

for (i in 1:iter) {
  y <- rexp(n=sample_len, rate = rate)
  sim_mat[i,] <- std_res(y=y)
} 

res_avg <- colMeans(sim_mat)
print(res_avg)
plot(res_avg)
title("Exponential Distribution Residual Plot")
print(mean(res_avg))