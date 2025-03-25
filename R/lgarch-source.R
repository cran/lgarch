###########################################################
## This file contains the lgarch source code
##
## First build/version: 18 March 2014 (version 0.1)
##
## -- Part 1: G-functions
## NOT INCLUDED: rsep       #generate random draws from SEP density
## NOT INCLUDED: rsst       #generate random draws from ST density
## glag                     #Genaro's lag-function
## gdiff                    #Genaro's diff-function
##
## -- Part 2: univariate log-GARCH functions
## lgarchSim                #simulate from univariate log-garch(P,Q)
## lgarchRecursion1         #recursion
## lgarchObjective          #objective function
## lgarch                   #estimate; return object from lgarch class
##
## -- Part 3: S3 methods for lgarch class
## coef.lgarch              #S3 methods
## fitted.lgarch
## logLik.lgarch
## print.lgarch
## NOT INCLUDED YET: predict.lgarch
## residuals.lgarch
## rss
## summary.lgarch
## vcov.lgarch
##
##
###########################################################


###########################################################
## Part 1: G-FUNCTIONS
###########################################################

##Generate random draws from standardised Skew GED/EDP density
#rsep <- function(n, shape=2, skew=1, zero.prob=0)
#{
#if(shape==2){
#  z1 <- rnorm(n)
#}else{
#  lambda <- sqrt(2^(-2/shape) * gamma(1/shape)/gamma(3/shape))
#  r <- rgamma(n, 1/shape)
#  z1 <- lambda * (2 * r)^(1/shape) * sign(runif(n) - 1/2)
#}
#
#if(skew==1){
#  result <- z1
#}else{
#  weight <- skew/(skew + 1/skew)
#  z2 <- runif(n, -weight, 1 - weight)
#  signz2 <- sign(z2)
#  Xi <- skew^signz2
#  vRand <- -abs(z1)/Xi*signz2
#  g <- shape/(lambda * (2^(1 + 1/shape)) * gamma(1/shape))
#  m1 <- 2^(1/shape) * lambda * gamma(2/shape)/gamma(1/shape)
#  mu <- m1 * (skew - 1/skew)
#  sigma <- sqrt((1 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - 1)
#  vRand <- (vRand - mu)/sigma
#  result <- vRand
#}
#
#if(zero.prob!=0){
#  Ivar <- runif(n, -zero.prob, 1-zero.prob)
#  Ivar <- (Ivar > 0)
#  result <- result*Ivar/sqrt(1-zero.prob)
#}
#
#return(result)
#} #end rsep
#
##Generate random draws from standardised ST density
#rsst <- function(n, df=10, skew=1, zero.prob=0)
#{
#s2 <- sqrt(df/(df - 2))
#zstar <- rt(n=n,df=df)/s2
#
#if(skew==1){
#  result <- zstar
#}else{
#  if(!exists("beta")){
#    beta <- function(a,b){ exp( lgamma(a)+lgamma(b)-lgamma(a+b) ) }
#  }
#  weight <- skew/(skew + 1/skew)
#  z <- runif(n, -weight, 1 - weight)
#  Xi <- skew^sign(z)
#  Random <- -abs(zstar)/Xi * sign(z)
#  m1 <- 2 * sqrt(df - 2)/(df - 1)/beta(1/2, df/2)
#  mu <- m1 * (skew - 1/skew)
#  sigma <- sqrt((1 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - 1)
#  result <- (Random - mu)/sigma
#}
#
#if(zero.prob!=0){
#  Ivar <- runif(n, -zero.prob, 1-zero.prob)
#  Ivar <- (Ivar > 0)
#  result <- result*Ivar/sqrt(1-zero.prob)
#}
#
#return(result)
#
#} #end rsst

#################
glag <- function(x, k=1, pad=TRUE, pad.value=NA)
{
  #check arguments:
  if(k < 1) stop("Lag order k cannot be less than 1")

  #zoo-related:
  iszoo.chk <- is.zoo(x)
  x <- as.zoo(x)
  x.index <- index(x)
  x <- coredata(x)
  isvector <- is.vector(x)
  x <- cbind(x)
  x.n <- NROW(x)
  x.ncol <- NCOL(x)

  #do the lagging:
  x.nmink <- x.n - k
  xlagged <- matrix(x[1:x.nmink,], x.nmink, x.ncol)
  if(pad){
    xlagged <- rbind( matrix(pad.value,k,x.ncol) , xlagged)
  }

  #transform to vector?:
  if(x.ncol==1 && isvector==TRUE){
    xlagged <- as.vector(xlagged)
  }

  #if(is.zoo(x)):
  if(iszoo.chk){
    if(pad){
      xlagged <- zoo(xlagged, order.by=x.index)
    }else{
      xlagged <- zoo(xlagged, order.by=x.index[c(k+1):x.n])
    } #end if(pad)
  } #end if(iszoo.chk)

  #out:
  return(xlagged)
} #end glag

#glag <- cmpfun(glag)

#################
gdiff <- function(x, lag=1, pad=TRUE, pad.value=NA)
{
  #check arguments:
  if(lag < 1) stop("lag argument cannot be less than 1")

  #zoo-related:
  iszoo.chk <- is.zoo(x)
  x <- as.zoo(x)
  x.index <- index(x)
  x <- coredata(x)
  isvector <- is.vector(x)
  x <- cbind(x)
  x.n <- NROW(x)
  x.ncol <- NCOL(x)

  #do the differencing:
  xdiff <- x - glag(x, k=lag, pad=TRUE, pad.value=NA)

  #pad operations:
  if(!pad){
    xdiff <- na.trim(as.zoo(xdiff))
    xdiff <- coredata(xdiff)
  }else{
    whereisna <- is.na(xdiff)
    xdiff[whereisna] <- pad.value
  }

  #transform to vector?:
  if(x.ncol==1 && isvector==TRUE){
    xdiff <- as.vector(xdiff)
  }

  #if(is.zoo(x)):
  if(iszoo.chk){
    if(pad){
      xdiff <- zoo(xdiff, order.by=x.index)
    }else{
      xdiff <- zoo(xdiff, order.by=x.index[c(lag+1):x.n])
    } #end if(pad)
  } #end if(iszoo.chk)

  #out:
  return(xdiff)
} #end gdiff


###########################################################
## Part 2: log-garch functions
###########################################################

#################
lgarchSim <- function(n, constant=0, arch=0.05, garch=0.9, xreg=NULL,
  backcast.values=list(lnsigma2=NULL, lnz2=NULL, xreg=NULL),
  check.stability=TRUE, innovations=NULL, verbose=FALSE,
  c.code=TRUE)
{
  #check arguments:
  if(is.null(constant)) constant <- 0
  if(is.null(arch)) arch <- 0
  if(is.null(garch)) garch <- 0
  #if(is.null(asym)) asym <- 0

  #orders:
  arch.p <- length(arch)
  garch.q <- length(garch)
  maxpq <- max(arch.p, garch.q)
  nmaxpq <- n+maxpq
  #### NOTE: asym not implemented yet
  #asym.r <- length(asym)
  #maxpqr <- max(arch.p, garch.q, asym.r)
  #nmaxpqr <- n+maxpqr
  ####
  #### NOTE: c.code not implemented yet

  #make phi:
  phi <- c(arch,rep(0,maxpq-arch.p)) + c(garch,rep(0,maxpq-garch.q))

  #check stability:
  if(check.stability){
    roots <- polyroot(c(1,-phi))
    if(any(abs(roots) <= 1)){
      mssg <- paste("NOTE: The log-garch model may not be stable (one or more AR-roots is on or inside the unit circle)")
      print(mssg)
    }
  }

  #xreg, zoo-index:
  if(is.null(xreg)){
    result.index <- 1:n
    xreg <- rep(0,nmaxpq)
    xregMean <- 0
  }else{
    xreg <- as.zoo(xreg)
    result.index <- index(xreg)
    xreg <- coredata(xreg)
    if(is.null(backcast.values$xreg)){
      xregMean <- mean(xreg)
      xreg.init <- rep(xregMean,maxpq)
    }else{
      xreg.init <- backcast.values$xreg
    }
    xreg <- c(xreg.init, xreg)
  }

  #innovation:
  if(is.null(innovations)){
    z <- rnorm(n)
  }else{
    z <- innovations
  }
  z2 <- z^2
  lnz2 <- log(z2)
  if(is.null(backcast.values$lnz2)){
    Elnz2 <- mean(lnz2)
    lnz2 <- c(rep(Elnz2, maxpq), lnz2)
  }else{
    lnz2 <- c(backcast.values$lnz2, lnz2)
  }
  mLnz2 <- NULL
  for(k in 1:arch.p){
    mLnz2 <- cbind(mLnz2, glag(lnz2, k=k, pad=TRUE,
      pad.value=Elnz2))
  }

  ##innov series:
  if(arch.p > 0){
    lnz2.innov <- as.vector(mLnz2 %*% arch)
  }else{ lnz2.innov <- rep(0,nmaxpq) }
  innov <- constant + lnz2.innov + xreg

  ##lnsigma2:
  if(is.null(backcast.values$lnsigma2)){
    Elnsigma2 <- (constant+sum(arch)*Elnz2+xregMean)/(1-sum(phi))
    if(abs(Elnsigma2)==Inf){
      mssg <- paste("NOTE: The backcast value(s) of lnsigma2 is Inf")
      print(mssg)
    }
    lnsigma2 <- c(rep(Elnsigma2,maxpq),rep(0,n))
  }else{
    lnsigma2 <- c(backcast.values$lnsigma2,rep(0,n))
  }

  #recursion:
  maxpq1 <- maxpq+1
  phisum <- rep(0,nmaxpq)
  if( c.code ){
    tmp <- .C("LGARCHSIM", as.integer(maxpq), as.integer(nmaxpq),
      as.double(lnsigma2), as.double(phi), as.double(phisum), as.double(innov),
      PACKAGE = "lgarch")
    names(tmp) <- c("maxpq", "nmaxpq", "lnsigma2", "phi", "phisum", "innov")
    lnsigma2 <- tmp$lnsigma2
  }else{
    for(i in maxpq1:nmaxpq){
      phisum[i] <- sum(phi*lnsigma2[c(i-1):c(i-maxpq)])
      lnsigma2[i] <- innov[i] + phisum[i]
    } #end for loop
  } #end if(c.code)

  ##result:
  ##=======
  
  lnsigma2 <- lnsigma2[-c(1:maxpq)] #rm initial vals
  if(verbose){
    sigma <- exp(lnsigma2/2)
    y <- sigma*z
    lny2 <- log(y^2)
    lnz2 <- log(z^2)
    result <- cbind(y, lny2, sigma, lnsigma2, z, lnz2)
    colnames(result)[3:4] <- c("sd","lnsd2")
  }else{
    sigma <- exp(lnsigma2/2)
    result <- sigma*z
  }
  result <- zoo(result, order.by=result.index)
  return(result)

} #close lgarchSim()

#################
lgarchRecursion1 <- function(pars, aux)
{

  ##make vectors:
  ##=============
  if(aux$xreg.k == 0){
    innov <- pars[1]
    innovMean <- mean(innov)
    innov <- rep(innov, max(aux$n+1,aux$nmaxpq))
  }else{
    innov <- pars[1] + aux$xreg %*% pars[aux$xreg.indx]
    innovMean <- mean(innov)
    innov <- c(rep(innovMean, max(1,aux$maxpq)), innov)
  }
  if(aux$mean.correction){
    lny2adj <- c(rep(aux$Elny2,1:max(1,aux$maxpq)), aux$lny2mc)
  }else{
    lny2adj <- c(rep(aux$Elny2,1:max(1,aux$maxpq)), aux$lny2)
  }
  uadj <- aux$zerosaux

  ##recursion:
  ##==========
  if(aux$maxpq > 0){ phi1 <- pars[aux$ar.indx] }else{ phi1 <- 0 }
  if(aux$ma > 0){ theta1 <- pars[aux$ma.indx] }else{ theta1 <- 0 }
  if(aux$c.code){
    iStart <- max(2,aux$maxpq+1) - 1
    iEnd <- max(aux$nmaxpq, aux$n+1)
    tmp <- .C("ARMARECURSION1", as.integer(iStart), as.integer(iEnd),
      as.numeric(phi1), as.numeric(theta1), as.numeric(aux$yzeroadj),
      as.numeric(innov), as.numeric(lny2adj), as.numeric(uadj),
      PACKAGE="lgarch")
      names(tmp) <- c("iStart", "iEnd", "phi1", "theta1", "yzeroadj",
        "innov", "lny2adj", "uadj")
    uadj <- tmp$uadj
  }else{
    for(i in max(2,aux$maxpq+1):max(aux$nmaxpq, aux$n+1)){
      imin1 <- i - 1
      if(aux$yzeroadj[i]==0){
        lny2adj[i] <- innov[i] + phi1*lny2adj[imin1] + theta1*uadj[imin1]
        uadj[i] <- 0
      }else{
        uadj[i] <- lny2adj[i] - innov[i] - phi1*lny2adj[imin1] - theta1*uadj[imin1]
      } #end if(not-0)
    } #end for loop
  } #end if(c.code)

  ##result:
  ##=======
  if(aux$verboseRecursion){
    if(aux$c.code){
      result <- cbind(tmp$uadj, tmp$lny2adj)
    }else{
      result <- cbind(uadj, lny2adj)
    } #end if(c.code)
    result <- result[-c(1:max(1,aux$maxpq)),]
  }else{
    result <- uadj[-c(1:max(1,aux$maxpq))]
  } #end if(verboseRecursion)
  return(result)
  
} #close lgarchRecursion1()

#################
lgarchObjective <- function(pars, aux)
{
  #check parameters:
  if( any(is.na(pars)) || any(pars<=aux$lower) || any(pars>=aux$upper) ){
    chk.conds <- FALSE
  }else{
    chk.conds <- TRUE
  } #end if..else..

  if(chk.conds){
    #if mean-correction:
    if(aux$mean.correction){
      pars <- c(0,pars)
    }
    #recursion:
    uadj <- lgarchRecursion1(pars, aux)
    if(aux$yzeron > 0){
      uadj <- uadj[-aux$yzerowhere]
    }
    #compute objective value:
    if(aux$method=="ls"){
      objective.value <- sum(uadj^2)
    }
    if(aux$method=="ml"){
      sigma2u <- pars[aux$sigma2u.indx]
      objective.value <- -aux$ynonzeron*log(sigma2u)/2 - aux$ynonzeron*log(2*pi)/2 - sum(uadj^2/sigma2u)/2
    }
    if(aux$method=="cex2"){
      uadjElnz2 <- uadj + pars[aux$sigma2u.indx]
      objective.value <- -aux$ynonzeron*log(2*pi)/2 + sum(uadjElnz2 - exp(uadjElnz2))/2
    }
    #check objective value:
    if(is.na(objective.value) || abs(objective.value) == Inf){
      objective.value <- aux$objective.penalty
    }
  }else{
    objective.value <- aux$objective.penalty
  } #end if(chk.conds)
  return(objective.value)
} #end lgarchObjective

#################
lgarch <- function(y, arch=1, garch=1, xreg=NULL,
  initial.values=NULL, lower=NULL, upper=NULL,
  nlminb.control=list(), vcov=TRUE, method=c("ls","ml","cex2"),
  mean.correction=FALSE, objective.penalty=NULL,
  solve.tol=.Machine$double.eps, c.code=TRUE)
{
  #check/change arguments:
  if(arch < garch) stop("garch order cannot be greater than arch order, since estimation is via the arma representation")
  if(arch > 1) stop("Sorry, arch order cannot be greater than 1 in the current version of lgarch")
  method <- match.arg(method)
  if(method=="cex2"){ mean.correction <- TRUE }
  if(mean.correction==TRUE && method=="ls" && arch==0 && garch==0 && is.null(xreg) )
    stop("This combination is not possible. Try setting mean.correction=FALSE")

  #zoo and xts specific:
  y.name <- deparse(substitute(y))
  y <- as.zoo(cbind(y))
  if( is.null(y.name)){ y.name <- colnames(y)[1] }
  if( y.name[1] =="" ){ y.name <- "y" }
  y <- na.trim(y)
  y.n <- NROW(y)
  y.index <- index(y)
  t1 <- y.index[1]
  t2 <- y.index[y.n]
  y <- coredata(y)
  if(NCOL(y) > 1){
    stop("Dependent variable not 1-dimensional")
  }else{
    y <- y[,1]
  }

  ##xreg:
  if(!is.null(xreg)){
    xreg <- as.zoo(cbind(xreg))
    xreg.names <- colnames(xreg)
    if(is.null(xreg.names)){
      xreg.names <- paste("xreg", 1:NCOL(xreg), sep="")
    }
    if(any(xreg.names == "")){
      missing.colnames <- which(xreg.names == "")
      for(i in 1:length(missing.colnames)){
        xreg.names[i] <- paste("xreg", i, sep="")
      }
    }
    xreg <- window(xreg, start=t1, end=t2)
    xreg <- cbind(coredata(xreg))
  }

  ##begin aux list:
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
  if(mean.correction){ aux$lny2mc <- aux$lny2 - aux$Elny2 }

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
  if(method!="ls"){
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
    sigma2u <- NULL
    if(method=="ml"){ sigma2u <- 4.94 }
    if(method=="cex2"){ sigma2u <- -1.27 }
    initial.values <- c(constant, ar.initvals, ma.initvals,
      xreg.initvals, sigma2u)
  }else{
    max.indx <- max( c(1,aux$ar.indx,aux$ma.indx,aux$xreg.indx,aux$sigma2u.indx) )
    if( length(initial.values)!=max.indx ){
      stop("length(initial.values) not equal to no. of parameters to be estimated")
    } #end check length
  } #end if..else..

  #upper bounds:
  if(is.null(upper)){
    upper <- c(Inf, rep(1-.Machine$double.eps,aux$ar),
      rep(1-.Machine$double.eps,aux$ma), rep(Inf, aux$xreg.k))
    if(method!="ls"){ upper <- c(upper,Inf) }
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
    if(method=="cex2"){ lower <- c(lower,-Inf) }
  }else{
    if( length(lower)!=length(initial.values) )
      stop("length(lower) not equal to length(initial.values)")
  }
  aux$lower <- lower

  #misc:
  aux$c.code <- c.code
  aux$solve.tol <- solve.tol
  aux$mean.correction <- mean.correction
  aux$verboseRecursion <- FALSE
  aux$yzeroadj <- c(rep(1,max(1,aux$maxpq)), aux$yzero)
  aux$zerosaux <- rep(0,max(aux$nmaxpq, aux$n+1))
  if(mean.correction){
    initial.values <- initial.values[-1]
    aux$upper <- aux$upper[-1]
    aux$lower <- aux$lower[-1]
  }
  aux$initial.values.arma <- initial.values
  if(is.null(objective.penalty)){
    aux$objective.penalty <- lgarchObjective(initial.values, aux)
  }else{
    aux$objective.penalty <- objective.penalty
  }

  #estimate:
  if(method=="ls"){
    objective.f <- function(pars, x=aux){ lgarchObjective(pars,x) }
  }else{
    objective.f <- function(pars, x=aux){ -lgarchObjective(pars,x) }
  }
  est <- nlminb(initial.values, objective.f, lower=aux$lower,
    upper=aux$upper, control=nlminb.control)

  #post-estimation:
  if(method!="ls") est$objective <- -est$objective
  if(mean.correction){
    est$par <- c(0,est$par)
    if(aux$ar > 0){
      arsum <- sum(est$par[aux$ar.indx])
    }else{
      arsum <- 0
    }
    est$par[1] <- (1-arsum)*aux$Elny2
  }
  names(est)[2] <- "objective.arma"

  #Elnz2 and sigma2u parameters:
  if(mean.correction){
    uadj <- lgarchRecursion1(c(0,est$par[-1]), aux)
  }else{
    uadj <- lgarchRecursion1(est$par, aux)
  }
  if(aux$yzeron > 0){
    uadj <- uadj[-aux$yzerowhere]
  }
  if(method=="ls"){
    Elnz2 <- -log(mean(exp(uadj)))
    namesArma <- NULL
    sigma2u <- var(uadj)
  }
  if(method=="ml"){
    Elnz2 <- -log(mean(exp(uadj - mean(uadj))))
    namesArma <- "sigma2u"
    sigma2u <- est$par[aux$sigma2u.indx]
  }
  if(method=="cex2"){
    Elnz2<- est$par[aux$sigma2u.indx]
    namesArma <- "Elnz2"
    sigma2u <- var(uadj)
  }
  par.lgarch <- Elnz2
  namesLgarch <- "Elnz2"
  par.arma <- est$par

  #xreg, garch, arch and intercept parameters:
  if(aux$xreg.k > 0){
    namesArma <- c(xreg.names, namesArma)
    namesLgarch <- c(xreg.names, namesLgarch)
    #OLD:
    #namesArma <- c(paste("xreg",1:aux$xreg.k,sep=""), namesArma)
    #namesLgarch <- c(paste("xreg",1:aux$xreg.k,sep=""), namesLgarch)
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

  #intercept and completion:
  namesArma <- c("intercept", namesArma)
  #OLD: par.arma <- est$par
  names(par.arma) <- namesArma
  namesLgarch <- c("intercept", namesLgarch)
  const.lgarch <- est$par[1] - (1+sum(est$par[aux$ma.indx]))*Elnz2
  par.lgarch <- c(const.lgarch,par.lgarch)
  names(par.lgarch) <- namesLgarch
  est$par <- par.lgarch
  est <- c(list(date=date(), par.arma=par.arma), est)

  #vcov matrix:
  if(vcov){
    par.arma <- as.numeric(par.arma)
    if(mean.correction){
      par.arma <- par.arma[-1]
      namesArma <- namesArma[-1]
    }
    hessian.arma <- -optimHess(par.arma, objective.f)
    colnames(hessian.arma) <- namesArma
    rownames(hessian.arma) <- namesArma
    vcov.arma <- solve(-hessian.arma, tol=solve.tol)
    if(aux$method=="ls"){
      vcov.arma <- sigma2u*2*vcov.arma
    }
    est <- c(list(aux=aux, hessian.arma=hessian.arma,
      vcov.arma=vcov.arma, vcov.lgarch=NULL), est)
    est$vcov.lgarch <- vcov.lgarch(est, arma=FALSE)
  }else{
    est <- c(list(aux=aux, hessian.arma=NA,
      vcov.arma=NA, vcov.lgarch=NA), est)
  }

  #out:
  if(is.null(est$aux)){
    est <- c(list(aux=aux),est)
  }
  class(est) <- "lgarch"
  return(est)
} #end lgarch

###########################################################
## Part 3: S3 METHODS FOR LGARCH CLASS
###########################################################

##################
coef.lgarch <- function(object, arma=FALSE, ...)
{
  if(arma){
    result <- object$par.arma
  }else{
    result <- object$par
  }
  return(result)
} #end coef.lgarch

################### fitted values
fitted.lgarch <- function(object, verbose=FALSE, ...)
{
  object$aux$verboseRecursion <- TRUE
  pars <- as.numeric(object$par.arma)
  if(object$aux$mean.correction){ pars[1] <- 0 }
  mUhat <- lgarchRecursion1(pars, object$aux)
  lnz2adj <- mUhat[,1] + object$par["Elnz2"]
  if(object$aux$mean.correction){
    logsigma2 <- mUhat[,2] + object$aux$Elny2 - lnz2adj
  }else{
    logsigma2 <- mUhat[,2] - lnz2adj
  }
  sigma <- exp(logsigma2/2)
  if(verbose){
    zhat <- object$aux$y/sigma
    result <- cbind(sigma, logsigma2, zhat, mUhat[,1])
    colnames(result) <- c("sd", "lnsd2", "zhat", "uhat")
  }else{
    result <- sigma
  }
  result <- zoo(result, order.by=object$aux$y.index)
  return(result)
} #close fitted.lgarch()

#################### check if object is of class lgarch
#is.lgarch <- function(x)
#{
#if(class(x)=="lgarch"){ TRUE }else{ FALSE}
#}

################### return log-likelihood
logLik.lgarch <- function(object, arma=FALSE, ...)
{
  if(arma==FALSE){
    sigma <- fitted.lgarch(object)
    result <- sum(object$aux$yzero*dnorm(object$aux$y,
      sd=sigma, log=TRUE))
    attr(result, "df") <- length(object$par)-1
  }
  if(arma==TRUE && object$aux$method!="ls"){
    result <- object$objective.arma
    attr(result, "df") <- length(object$par.arma)-1
  }
  if(arma==TRUE && object$aux$method=="ls"){
    uhat <- lgarchRecursion1(as.numeric(object$par.arma),
      object$aux)
    if(object$aux$yzeron > 0){
      uhat <- uhat[-object$aux$yzerowhere]
    }
    result <- sum(dnorm(uhat, sd=sd(uhat), log=TRUE))
    attr(result, "df") <- length(object$par.arma)
  }
  attr(result, "nobs") <- length(object$aux$ynonzeron)
  class(result) <- "logLik"
  return(result)
} #end logLik.lgarch

################## forecast
#predict.lgarch <- function(object, n.ahead=1,
#  initial.values=NULL, n.sim=10000, verbose=FALSE, ...)
#{
#return(out)
#} #end predict.lgarch

################# print estimation result
print.lgarch <- function(x, arma=FALSE, ...)
{
  #out1:
  pars <- coef.lgarch(x, arma=arma)
  vcovmat <- vcov.lgarch(x, arma=arma)
  out1 <- rbind(pars, sqrt(diag(vcovmat)))
  rownames(out1) <- c("Estimate:", "Std. Error:")

  #out2:
  out2 <- as.data.frame(matrix(NA,4,1))
  out2[1,1] <- as.character(round(logLik.lgarch(x, arma=FALSE), digits=3))
  out2[2,1] <- as.character(round(logLik.lgarch(x, arma=TRUE), digits=3))
  out2[3,1] <- as.character(round(rss(x), digits=3))
  out2[4,1] <- as.character(round(x$aux$ynonzeron, digits=0))
  out2[5,1] <- as.character(round(x$aux$yzeron, digits=0))
  rownames(out2) <-   c("Log-likelihood (log-garch):",
    "Log-likelihood (arma):", "Sum Squared Resids. (arma):",
    "No. of non-zeros:", "No. of zeros:")
  colnames(out2) <- ""

  #header:
  cat("\n")
  cat("Date:", x$date, "\n")
  if(x$aux$method=="ls"){ methodmssg <- "LS" }
  if(x$aux$method=="ml"){ methodmssg <- "ML" }
  if(x$aux$method=="cex2"){ methodmssg <- "CEX2" }
  if(x$aux$mean.correction){
    methodmssg <- paste(methodmssg,
      " (w/mean-correction)", sep="")
  }
  cat("Method:", methodmssg, "\n")
  cat("Message (nlminb):", x$message, "\n")
  cat("No. of observations:", x$aux$n, "\n")
  cat("Sample:", as.character(x$aux$y.index[1]),
    "to", as.character(x$aux$y.index[x$aux$n]), "\n")
  cat("\n")
  cat("Coefficients:\n")
  print(out1)
  print(out2)
  cat("\n")
} #end print.lgarch

##################
residuals.lgarch <- function(object, arma=FALSE, ...)
{
  if(arma){
    pars <- as.numeric(object$par.arma)
    if(object$aux$mean.correction){ pars[1] <- 0 }
    result <- lgarchRecursion1(pars, object$aux)
  }else{
    sigma <- coredata(fitted.lgarch(object))
    result <- object$aux$y/sigma
  }
  result <- zoo(result, order.by=object$aux$y.index)
  return(result)
} #end residuals.lgarch

################# extract sum of squared residuals (arma)
rss <- function(object, ...)
{

  if( is(object)!="lgarch" && is(object)!="armax" )
    stop("object does not belong to the lgarch class")
#  OLD:
#  if(class(object)!="lgarch" && class(object)!="armax")
#    stop("object does not belong to the lgarch class")
  
  if( is(object)=="lgarch" ){
#  OLD:
#  if(class(object)=="lgarch"){
    pars <- as.numeric(object$par.arma)
    if(object$aux$mean.correction){ pars[1] <- 0 }
    uhat <- lgarchRecursion1(pars, object$aux)
    if(object$aux$yzeron > 0){
      uhat <- uhat[-object$aux$yzerowhere]
    }
  }
  
  return(sum(uhat^2))
} #end rss

################# summarise output
summary.lgarch <- function(object, ...)
{
  object.name <- deparse(substitute(object))
  xnames <- cbind(names(object))
  xrows <- nrow(xnames)
  colnames(xnames) <- ""
  rownames(xnames) <- rep(" ",xrows)
  cat("\n")
  cat("Items in list '", object.name, "':", sep="")
  cat("\n")
  print(xnames)
  cat("\n")
} #end summary.lgarch

################## return variance-covariance matrix
vcov.lgarch <- function(object, arma=FALSE, ...)
{
  #vcovarma:
  if(is.null(object$vcov.arma)){
    cat("vcov cannot be extracted (vcov=FALSE during estimation)\n")
    m.dim <- length(object$par.arma)
    vcovarma <- matrix(NA,m.dim,m.dim)
  }else{
    vcovarma <- object$vcov.arma
    if(object$aux$mean.correction){
      vcovarma <- cbind(rep(NA,NCOL(vcovarma)), vcovarma)
      vcovarma <- rbind(rep(NA,NCOL(vcovarma)), vcovarma)
    }
  }
  colnames(vcovarma) <- names(object$par.arma)
  rownames(vcovarma) <- names(object$par.arma)

  #vcovlgarch:
  if(arma){
    result <- vcovarma
  }else{

    if(!is.null(object$vcov.lgarch)){
      result <- object$vcov.lgarch
    }else{

      #make vcovlgarch matrix:
      vcovlgarch <- matrix(NA,length(object$par),
        length(object$par))
      colnames(vcovlgarch) <- names(object$par)
      rownames(vcovlgarch) <- names(object$par)

      #if garch:
      if(object$aux$ma > 0){
        vcovlgarch[object$aux$ma.indx, object$aux$ma.indx] <- vcovarma[object$aux$ma.indx, object$aux$ma.indx]
      } #end if(..ma > 0)

      #if arch:
      if(object$aux$ar > 0){
        ar.ma.diff <- object$aux$ar - object$aux$ma
        if(ar.ma.diff==0){
          varalpha <- vcovarma[2,2] + vcovarma[3,3] + 2*vcovarma[2,3]
          vcovlgarch[2,2] <- varalpha
          covalphabeta <- -vcovarma[2,3] - vcovarma[3,3]
          vcovlgarch[2,3] <- covalphabeta
          vcovlgarch[3,2] <- covalphabeta
        }else{
          vcovlgarch[2,2] <- vcovarma[2,2]
        } #end if(..diff==0)
      } #end if(..ar > 0)

      #if xreg:
      if(object$aux$xreg.k > 0){
        vcovlgarch[object$aux$xreg.indx, object$aux$xreg.indx] <- vcovarma[object$aux$xreg.indx, object$aux$xreg.indx]
        if(object$aux$ma > 0){
          vcovlgarch[object$aux$ma.indx, object$aux$xreg.indx] <- -vcovarma[object$aux$ma.indx,object$aux$xreg.indx]
          vcovlgarch[object$aux$xreg.indx, object$aux$ma.indx] <- -vcovarma[object$aux$xreg.indx, object$aux$ma.indx]
        }
        if(object$aux$ar > 0){
          vcovlgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcovarma[object$aux$ar.indx,object$aux$xreg.indx]
          if(object$aux$ma > 0){
            vcovlgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcovlgarch[object$aux$ar.indx, object$aux$xreg.indx] + vcovlgarch[object$aux$ma.indx, object$aux$xreg.indx]
          }
          vcovlgarch[object$aux$xreg.indx, object$aux$ar.indx] <- vcovlgarch[object$aux$ar.indx, object$aux$xreg.indx]
        } #end if(..ar > 0)
      } #end if(..k > 0)

      #Var(Elnz2^hat):
      if(object$aux$method=="cex2"){
        vcovlgarch[NROW(vcovlgarch), NCOL(vcovlgarch)] <- vcovarma[NROW(vcovarma), NCOL(vcovarma)]
      }else{

        zhat <- coredata(residuals.lgarch(object, verbose=TRUE))
        if(object$aux$yzeron > 0){
          zhat <- zhat[-object$aux$yzerowhere]
        }  #end if..
        zhat2 <- zhat^2
        avar <- var(zhat2 - log(zhat2))
        vcovlgarch[NROW(vcovlgarch), NCOL(vcovlgarch)] <- avar/length(zhat)

#        old:
#        uadj <- lgarchRecursion1(as.numeric(object$par.arma), object$aux)
#        if(object$aux$yzeron > 0){
#          uadj <- uadj[-object$aux$yzerowhere]
#        }  #end if..
#        expuadj <- exp(uadj)
#        uexpuadj <- uadj*exp(uadj)
#        avaruadj <- var(uadj) + var(expuadj)/mean(expuadj)^2 - 2*mean(uexpuadj)/mean(expuadj)
#        vcovlgarch[NROW(vcovlgarch), NCOL(vcovlgarch)] <- avaruadj/length(uadj)

      } #end if("cex2")else(..)

      #TO DO:
      #intercept:
      #if(object$aux$method=="cex2"){
      #} #end if("cex2")else(..)

      result <- vcovlgarch

    } #end if(!is.null(vcovlgarch))else(..)
  } #end if(arma)else(..)

  #out:
  return(result)
} #end vcov.lgarch
