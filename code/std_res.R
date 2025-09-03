### Standardized Residual Function for empirical process

std_res <- function(y){
  tryCatch(
    {
    n <- length(y)
    y_s <- (y - mean(y)) / sd(y)
    y_ss <- sort(y_s, decreasing = FALSE)
    phi_y <- pnorm(q=y_ss,mean=0,sd=1)
    i <- 1:20
    r <- phi_y - (i/(n+1))
    return(r)
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