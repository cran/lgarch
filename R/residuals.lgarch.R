residuals.lgarch <-
function(object, arma=FALSE, ...)
{
  aux <- object$aux
  pars <- as.numeric(object$par.arma)
  if(arma){
    result <- lgarchRecursion1(pars, aux)
  }else{
    aux$verboseRecursion <- TRUE
    mUhat <- lgarchRecursion1(pars, aux)
    uhatnotzero <- as.numeric(mUhat[,1]!=0)
    lnz2adj <- mUhat[,1] + uhatnotzero*as.numeric(object$par["Elnz2"])
    logsigma2 <- mUhat[,2] - lnz2adj
    sigma <- exp(logsigma2/2)
    result <- aux$y/sigma
  }
  result <- zoo(result, order.by=aux$y.index)
  return(result)
}
