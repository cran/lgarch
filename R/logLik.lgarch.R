logLik.lgarch <-
function(object, arma=FALSE, ...)
{
  if(arma){
    result <- object$objective
  }else{
    sigma <- fitted.lgarch(object)
    result <- sum(object$aux$yzero*dnorm(object$aux$y,
      sd=sigma, log=TRUE))
  }
  attr(result, "nobs") <- length(object$aux$ynonzeron)
  attr(result, "df") <- length(object$par)
  class(result) <- "logLik"
  return(result)
}
