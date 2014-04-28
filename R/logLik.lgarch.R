logLik.lgarch <-
function(object, arma=FALSE, ...)
{
  if(arma==FALSE){
    sigma <- fitted.lgarch(object)
    result <- sum(object$aux$yzero*dnorm(object$aux$y,
      sd=sigma, log=TRUE))
  }
  if(arma==TRUE && object$aux$method=="ml"){
    result <- object$objective
  }
  if(arma==TRUE && object$aux$method=="ls"){
    uhat <- lgarchRecursion1(as.numeric(object$par.arma),
      object$aux)
    if(object$aux$yzeron > 0){
      uhat <- uhat[-object$aux$yzerowhere]
    }
    result <- sum(dnorm(uhat, sd=sd(uhat), log=TRUE))
  }
  attr(result, "nobs") <- length(object$aux$ynonzeron)
  attr(result, "df") <- length(object$par)
  class(result) <- "logLik"
  return(result)
}
