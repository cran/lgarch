fitted.lgarch <-
function(object, verbose=FALSE, ...)
{
  aux <- object$aux
  aux$verboseRecursion <- TRUE
  pars <- as.numeric(object$par.arma)
  mUhat <- lgarchRecursion1(pars, aux)
  uhatnotzero <- as.numeric(mUhat[,1]!=0)
  lnz2adj <- mUhat[,1] + uhatnotzero*as.numeric(object$par["Elnz2"])
  logsigma2 <- mUhat[,2] - lnz2adj
  sigma <- exp(logsigma2/2)
  if(verbose){
    zhat <- aux$y/sigma
    result <- cbind(sigma, logsigma2, zhat, mUhat[,1])
    colnames(result) <- c("sigma", "logsigma2", "zhat", "uhat")
  }else{
    result <- sigma
  }
  result <- zoo(result, order.by=aux$y.index)
  return(result)
}
