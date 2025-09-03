# GPD L-moments

l_moments_gpd <- function(a,b,d,l){
  tryCatch(
    {
      l_temp <- matrix(0,nrow=1,ncol=4)
      l_temp[1,1] <- a - (b*(1-2*d)/(l+1))
      l_temp[1,2] <- b / ((l+1)*(l+2))
      l_temp[1,3] <- (((l-1)*(1-2*d))/(l+3))
      l_temp[1,4] <- (((l-1)*(l-2))/((l+3)*(l+4)))
      print(l_temp)
      return(l_temp)
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

