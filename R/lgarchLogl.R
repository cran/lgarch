lgarchLogl <-
function(pars, aux)
{
  if( any(is.na(pars)) || any(pars<=aux$lower) || any(pars>=aux$upper) ){
    chk.conds <- FALSE
  }else{
    chk.conds <- TRUE
  }

  if(chk.conds){
    uadj <- lgarchRecursion1(pars, aux)
    if(aux$yzeron > 0){
      uadj <- uadj[-aux$yzerowhere]
    }

    sigma2u <- pars[aux$sigma2u.indx]
    logl <- -aux$ynonzeron*log(sigma2u)/2 - aux$ynonzeron*log(2*pi)/2 - sum(uadj^2/sigma2u)/2

    if(is.nan(logl) || is.na(logl) || abs(logl) == Inf){
      logl <- aux$logl.penalty
    }
  }else{
    logl <- aux$logl.penalty
  } #end if(chk.conds)
  return(logl)
}
