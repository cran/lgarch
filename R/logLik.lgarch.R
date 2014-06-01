logLik.lgarch <-
function(object, arma=FALSE, ...)
{
  if(arma==FALSE){
    sigma <- fitted.lgarch(object)
    result <- sum(object$aux$yzero*dnorm(object$aux$y,
      sd=sigma, log=TRUE))
    attr(result, "df") <- length(object$par)-1
  }
  if(arma==TRUE && object$aux$method=="ml"){
    result <- object$objective
    attr(result, "df") <- length(object$par)
  }
  if(arma==TRUE && object$aux$method=="ls"){
    uhat <- lgarchRecursion1(as.numeric(object$par.arma),
      object$aux)
    if(object$aux$yzeron > 0){
      uhat <- uhat[-object$aux$yzerowhere]
    }
    result <- sum(dnorm(uhat, sd=sd(uhat), log=TRUE))
    attr(result, "df") <- length(object$par)
  }
  attr(result, "nobs") <- length(object$aux$ynonzeron)
  class(result) <- "logLik"
  return(result)
}
