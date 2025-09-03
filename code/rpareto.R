# Generalized Pareto Distribution

rpareto <- function(alpha,beta,delta,lambda,n){
  tryCatch(
    {
      p<-runif(n=n)
      r_p <- alpha + beta*((1-delta)*((p^lambda-1)/lambda)-delta*(((1-p)^lambda-1)/lambda))
      return(r_p)
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

