rss.lgarch <-
function(object, ...)
{
  if(object$aux$method=="ls"){
    result <- object$objective
  }else{
    uhat <- lgarchRecursion1(as.numeric(object$par.arma),
      object$aux)
    if(object$aux$yzeron > 0){
      uhat <- uhat[-object$aux$yzerowhere]
    }
    result <- sum(uhat^2)
  }
  return(result)
}
