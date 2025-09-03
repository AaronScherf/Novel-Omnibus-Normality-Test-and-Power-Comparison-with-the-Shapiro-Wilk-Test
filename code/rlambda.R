# Generalized Lambda Distribution

rlambda <- function(l1,l2,l3,l4,n){
  tryCatch(
    {
      p<-runif(n=n)
      r_p <- l1 + (p^l3 - (1-p)^l4)/l2
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

