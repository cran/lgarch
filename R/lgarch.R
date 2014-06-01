lgarch <-
function(y, arch=1, garch=1, xreg=NULL,
  initial.values=NULL,
  backcast.values=list(lny2=NULL, lnz2=NULL, xreg=NULL),
  lower=NULL, upper=NULL, nlminb.control=list(), vcov=TRUE,
  method="ls", objective.penalty=NULL, solve.tol=.Machine$double.eps,
  c.code=TRUE)
{
  #arguments:
  if(arch < garch) stop("garch order cannot be greater than arch order, since estimation is via the arma representation")
  if(arch > 1) stop("Sorry, arch order cannot be greater than 1 in the current version of lgarch")

  #zoo:
  y <- as.zoo(y)
  y <- na.trim(y)
  y.index <- index(y)
  y <- coredata(y)

  #xts specifics:
  if(is.matrix(y)){
    if(NCOL(y) > 1) stop("Dependent variable not a 1-dimensional matrix")
     y <- y[,1]
  }
  y <- as.numeric(y)

  #xreg:
  if(!is.null(xreg)) xreg <- cbind(coredata(xreg))

  #begin aux list:
  aux <- list()
  aux$method <- method
  aux$y <- y
  aux$y.index <- y.index
  aux$n <- length(y)
  aux$yzero <- as.numeric(y != 0)
  aux$ynonzeron <- sum(aux$yzero)
  aux$yzeron <- aux$n - aux$ynonzeron
  if(aux$yzeron > 0) aux$yzerowhere <- which(y == 0)
  y2 <- y^2
  if(aux$yzeron > 0){
    miny2 <- min(y2[-aux$yzerowhere])
    y2[aux$yzerowhere] <- miny2
    #delete?:
    #lny2 <- log(y2)
    aux$lny2 <- log(y2)
    aux$Elny2 <- mean(aux$lny2[-aux$yzerowhere])
  }else{
    aux$lny2 <- log(y2)
    aux$Elny2 <- mean(aux$lny2)
  }

  #orders:
  aux$maxpq <- max(arch,garch)
  aux$nmaxpq <- aux$n + aux$maxpq
  aux$ar <- aux$maxpq
  aux$ma <- garch
  if(aux$ar>0){ aux$ar.indx <- 2:c(aux$ar+1) }else{ aux$ar.indx <- 0 }
  if(aux$ma>0){
    aux$ma.indx <- c(max(aux$ar.indx)+1):c(max(aux$ar.indx)+aux$ma)
  }else{
    aux$ma.indx <- 0
  }
  if(is.null(xreg)){
    aux$xreg.k <- 0
  }else{
    aux$xreg <- xreg
    aux$xreg.k <- ncol(aux$xreg)
    aux$xreg.indx <- c(max(1,aux$ar.indx,aux$ma.indx)+1):c(max(1,aux$ar.indx,aux$ma.indx)+aux$xreg.k)
  }
  if(method=="ml"){
    aux$sigma2u.indx <- 1 + aux$ar + aux$ma + aux$xreg.k + 1
  }

  #initial values:
  if(is.null(initial.values)){
    if(is.null(xreg)){
      xreg.initvals <- numeric(0)
      xregMean <- 0
    }else{
      xreg.initvals <- rep(0.01, aux$xreg.k)
      xregMean <- mean(aux$xreg %*% xreg.initvals)
    }
    ma.initvals <- rep(-0.8/aux$ma, aux$ma)
    ar.initvals <- rep(0.9/aux$ar, aux$ar)
    #to do: check ar.initvals for stability?
    constant <- (1-sum(ar.initvals))*aux$Elny2 - xregMean
    if(method=="ml"){ sigma2u <- 4.94 }else{ sigma2u <- NULL }
    initial.values <- c(constant, ar.initvals, ma.initvals,
      xreg.initvals, sigma2u)
  }else{
    max.indx <- max( c(1,aux$ar.indx,aux$ma.indx,aux$xreg.indx,aux$sigma2u.indx) )
    if( length(initial.values)!=max.indx ){
      stop("length(initial.values) not equal to no. of parameters to be estimated")
    } #end check length
  } #end if..else..
  aux$initial.values.arma <- initial.values

  #upper bounds:
  if(is.null(upper)){
    upper <- c(Inf, rep(1-.Machine$double.eps,aux$ar),
      rep(1-.Machine$double.eps,aux$ma), rep(Inf, aux$xreg.k))
    if(method=="ml"){ upper <- c(upper,Inf) }
  }else{
    if( length(upper)!=length(initial.values) )
      stop("length(upper) not equal to length(initial.values)")
  }
  aux$upper <- upper

  #lower bounds:
  if(is.null(lower)){
    lower <- c(-Inf, rep(.Machine$double.eps-1,aux$ar),
      rep(.Machine$double.eps-1,aux$ma), rep(-Inf, aux$xreg.k))
    if(method=="ml"){ lower <- c(lower,0) }
  }else{
    if( length(lower)!=length(initial.values) )
      stop("length(lower) not equal to length(initial.values)")
  }
  aux$lower <- lower

  #misc:
  aux$c.code <- c.code
  aux$solve.tol <- solve.tol
  aux$verboseRecursion <- FALSE
  aux$yzeroadj <- c(rep(1,max(1,aux$maxpq)), aux$yzero)
  aux$zerosaux <- rep(0,max(aux$nmaxpq, aux$n+1))
  if(is.null(objective.penalty)){
    aux$objective.penalty <- lgarchObjective(initial.values, aux)
  }else{
    aux$objective.penalty <- objective.penalty
  }

  #estimate:
  if(method=="ml"){
    objective.f <- function(pars, x=aux){ -lgarchObjective(pars,x) }
  }else{
    objective.f <- function(pars, x=aux){ lgarchObjective(pars,x) }
  }
  est <- nlminb(initial.values, objective.f, lower=lower,
    upper=upper, control=nlminb.control)
  if(method=="ml") est$objective <- -est$objective
  names(est)[2] <- "objective.arma"

  #parameters:
  uadj <- lgarchRecursion1(as.numeric(est$par), aux)
  if(aux$yzeron > 0){
    uadj <- uadj[-aux$yzerowhere]
  }
  if(method=="ml"){ Elnz2 <- -log(mean(exp(uadj - mean(uadj)))) }
  if(method=="ls"){ Elnz2 <- -log(mean(exp(uadj))) }
  par.lgarch <- Elnz2
  namesLgarch <- "Elnz2"
  par.arma <- est$par
  if(method=="ml"){
    namesArma <- "sigma2u"
  }else{ namesArma <- NULL }
  if(aux$xreg.k > 0){
    namesArma <- c(paste("xreg",1:aux$xreg.k,sep=""), namesArma)
    namesLgarch <- c(paste("xreg",1:aux$xreg.k,sep=""), namesLgarch)
    par.lgarch <- c(est$par[aux$xreg.indx], par.lgarch)
  }
  if(aux$ma > 0){
    namesArma <- c(paste("ma",1:aux$ma,sep=""), namesArma)
    namesLgarch <- c(paste("garch",1:aux$ma,sep=""), namesLgarch)
    par.lgarch <- c(-est$par[aux$ma.indx], par.lgarch)
  }
  if(aux$ar > 0){
    namesArma <- c(paste("ar",1:aux$ar,sep=""), namesArma)
    namesLgarch <- c(paste("arch",1:aux$ar,sep=""), namesLgarch)
    par.tmp <- c(est$par[aux$ma.indx], rep(0, aux$ar-aux$ma))
    par.lgarch <- c(est$par[aux$ar.indx]+par.tmp, par.lgarch)
  }
  namesArma <- c("intercept", namesArma)
  names(par.arma) <- namesArma
  namesLgarch <- c("intercept", namesLgarch)
  const.lgarch <- est$par[1] - (1+sum(est$par[aux$ma.indx]))*Elnz2
  par.lgarch <- c(const.lgarch,par.lgarch)
  names(par.lgarch) <- namesLgarch
  est$par <- par.lgarch
  est <- c(list(date=date(), par.arma=par.arma), est)

  #vcov matrix:
  if(vcov){
    #rewrite (look at mlgarch code):
    hessian.arma <- -optimHess(as.numeric(par.arma), objective.f)
    colnames(hessian.arma) <- namesArma
    rownames(hessian.arma) <- namesArma
    vcov.arma <- solve(-hessian.arma, tol=solve.tol)
    est <- c(list(aux=aux, hessian.arma=hessian.arma,
      vcov.arma=vcov.arma, vcov.lgarch=NULL), est)
    #est$vcov.lgarch <- vcov.lgarch(est, ...)
  }

  #out:
  if(is.null(est$aux)){
    est <- c(list(aux=aux),est)
  }
  class(est) <- "lgarch"
  return(est)
}
