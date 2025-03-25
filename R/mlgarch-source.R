###########################################################
## This file contains the source code of the mlgarch model
##
## First created: 1 June 2014 (version 0.3)
##
## CONTENTS:
##
## -- PART 1: MULTIVARIATE LOG-GARCH FUNCTIONS
##
##   rmnorm                   #fast generation from a multivariate normal
##   NOT INCLUDED YET: dmnorm #multivariate normal density
##   mlgarchSim               #simulate from multivariate log-garch(1,1)
##   mlgarchRecursion1        #recursion
##   mlgarchObjective         #objective function
##   mlgarch                  #estimate; return object from lgarch class
##
## -- Part 2: S3 methods for mlgarch class
##
##   coef.mlgarch             #S3 methods
##   fitted.mlgarch
##   logLik.mlgarch
##   print.mlgarch
##   NOT INCLUDED YET: predict.mlgarch
##   residuals.mlgarch
##   summary.mlgarch
##   vcov.mlgarch
##
###########################################################


###########################################################
## PART 1: MULTIVARIATE LOG-GARCH FUNCTIONS
###########################################################

#################
rmnorm <- function(n, mean=NULL, vcov=1)
{
  d <- NCOL(vcov)
  y <- matrix(rnorm(n * d), d, n)
  y <- crossprod(y,chol(vcov))
  if(!is.null(mean)){ y <- t(mean + t(y)) }
  return(y)
} #end rmnorm

#################
mlgarchSim <- function(n, constant=c(0,0), arch=diag(c(0.1,0.05)),
  garch=diag(c(0.7,0.8)), xreg=NULL,
  backcast.values=list(lnsigma2=NULL, lnz2=NULL, xreg=NULL),
  innovations=NULL, innovations.vcov=diag(rep(1,length(constant))),
  check.stability=TRUE, verbose=FALSE)
{
  #check/change arguments:
  if(is.null(constant)){ constant <- 0 }
  if(is.null(arch)){ arch <- 0 }
  if(is.null(garch)){ garch <- 0 }
  c.code <- FALSE #necessary until c.code implemented

  #prepare:
  npluss1 <- n+1
  constant <- as.vector(constant)
  iRows <- NROW(constant)

  #xreg:
  if(is.null(xreg)){
    constantx <- matrix(constant,iRows,npluss1)
    Econstantx <- constant
  }else{
    if(!is.matrix(xreg)) stop("xreg is not a matrix")
    if(ncol(xreg)!=iRows) stop("ncol(xreg) must equal length(constant)")
    if(nrow(xreg)!=n) stop("nrow(xreg) must equal n")
    constantx <- constant + t(xreg)
    Econstantx <- as.numeric(rowMeans(constantx))
    constantx <- cbind(Econstantx, constantx) #backcast values (1 period)
  }

  #arch:
  if(length(arch) == 1 && arch==0){
    arch <- matrix(0,iRows,iRows)
  }

  #garch:
  if(length(garch) == 1 && garch==0){
    garch <- matrix(0,iRows,iRows)
  }

  #make phi:
  phi <- arch + garch

  #check stability:
  if(check.stability){
    if(any(abs(eigen(phi)$values) >= 1)){
      mssg <- paste("The model may not be stable (one or more AR-roots is on or inside the unit circle)")
      print(mssg)
    }
  }

  #z series:
  if(is.null(innovations)){
    #for the future: check if diag(innovations.vcov)==1?
    innovations <- rmnorm(n, mean=rep(0,iRows),
      vcov=innovations.vcov)
  }
  mlnz2 <- log(innovations^2)
  Elnz2 <- colMeans(mlnz2)
  if(is.null(backcast.values$lnz2)){
    #to do: account for no lags
    mlnz2 <- rbind(as.numeric(Elnz2),mlnz2)
  }else{
    mlnz2 <- rbind(as.numeric(backcast.values$lnz2),mlnz2)
  }

  #make lnsigma2:
  lnsigma2 <- matrix(NA,npluss1,iRows)
  mI <- diag(rep(1,iRows))
  innovInit <- as.matrix(Econstantx) + arch %*% cbind(as.numeric(Elnz2))
  lnsigma2[1,] <- solve(mI-phi) %*% innovInit

  #recursion:
  if(c.code){
    #not available yet
  }else{
    for(i in 2:npluss1){
      #future code:
      #phiSum <- for(..) phi[[..]] etc.
      #archSum <- for(..) arch[[..]] etc.
      #lnsigma2[,] <- phiSum + archSum
      lnsigma2[i,] <- constantx[,i] + phi %*% lnsigma2[i-1,] + arch %*% mlnz2[i-1,]
    } #end for loop
  } #end if(c.code)

  #out:
  sigma <- exp(lnsigma2[-1,]/2)
  y <- sigma*innovations
  if(verbose){
    y <- cbind(y, sigma, innovations)
    colnames(y) <- c(paste("y",1:iRows,sep=""),
      paste("sd",1:iRows,sep=""),
      paste("innov",1:iRows,sep="") )
  }else{
    colnames(y) <- paste("y",1:iRows,sep="")
  }
  as.zoo(y)
} #end mlgarchSim


#################
mlgarchRecursion1 <- function(pars, aux)
{
  ##matrices:
  ##=========
  lny2adj <- aux$lny2adj
  uadj <- matrix(0,aux$nmaxpq,aux$m)
  mInnov <- t(matrix(rep(pars[aux$const.indx],aux$n),aux$m,aux$n))
  if(aux$xreg.k > 0){
    xpars <- matrix(pars[aux$xreg.indx],aux$m,aux$xreg.m)
    mInnov <- mInnov + aux$xreg %*% t(xpars)
    if(aux$maxpq > 0){
      innovMeans <- colMeans(mInnov)
      mInnov <- rbind(t(matrix(rep(innovMeans,aux$maxpq),aux$m,aux$maxpq)),
        mInnov)
    } #end if(..maxpq > 0)
  }else{
    if(aux$maxpq > 0){
      innovMeans <- pars[aux$const.indx]
      mInnov <- rbind(t(matrix(rep(innovMeans,aux$maxpq),aux$m,aux$maxpq)),
        mInnov)
    } #end if(..maxpq > 0)
  } #end if(..$xreg.k > 0)
  if(aux$ar > 0){
    phi1 <- matrix(pars[aux$ar.indx],aux$m,aux$m)
  }else{ phi1 <- matrix(0,aux$m,aux$m) }
  if(aux$ma > 0){
    theta1 <- matrix(pars[aux$ma.indx],aux$m,aux$m)
  }else{ theta1 <- matrix(0,aux$m,aux$m) }

  ##recursion:
  ##==========
  if(aux$c.code){
    tmp <- .C("VARMARECURSION1", as.integer(aux$maxpq),
      as.integer(aux$nmaxpq), as.integer(aux$m),
      as.double(uadj), as.double(lny2adj), as.double(mInnov),
      as.double(phi1), as.double(theta1),
      as.double(aux$yiszeroadj), PACKAGE="lgarch")
    names(tmp) <- c("iStart", "n", "m", "mU", "mY", "mInnov",
      "PHI", "THETA", "mYiszeroadj")  
    uadj <- matrix(tmp$mU,aux$nmaxpq,aux$m)
  }else{
    for(i in c(aux$maxpq+1):aux$nmaxpq){
      Elny2adj <- mInnov[i,] + phi1%*%lny2adj[c(i-1),] + theta1%*%uadj[c(i-1),]
      if(aux$yanyrowiszeroadj[i]==1){
        for(j in 1:aux$m){
          if(aux$yiszeroadj[i,j]==1){
            lny2adj[i,j] <- Elny2adj[j]
          } #end if(is zero)
        } #end for(j..)
      }
      uadj[i,] <- lny2adj[i,] - Elny2adj
    } #end for(i..) loop
  } #end if(aux$c.code)

  ##result:
  ##=======
  if(aux$maxpq > 0){
    uadj <- uadj[-c(1:aux$maxpq),]
  }
  if(aux$verboseRecursion){
    if(aux$c.code){
      lny2adj <- matrix(tmp$mY,aux$nmaxpq,aux$m)
    } #end if(aux$c.code)
    if(aux$maxpq > 0){
      lny2adj <- lny2adj[-c(1:aux$maxpq),]
    } #end if(aux$maxpq > 0)
    uadj <- cbind(uadj, lny2adj)
    colnames(uadj) <- c(paste("u",1:aux$m,sep=""),
      paste("lny2no",1:aux$m,sep=""))
  } #end if(aux$verboseRecursion)
  return(uadj)

} #close mlgarchRecursion1()

#################
mlgarchObjective <- function(pars, aux)
{
  #check parameters:
  if( any(is.na(pars)) || any(pars<=aux$lower) || any(pars>=aux$upper) ){
    chk.conds <- FALSE
  }else{
    chk.conds <- TRUE
  } #end if..else..

  if(chk.conds){
    #recursion:
    uadj <- mlgarchRecursion1(pars, aux)
    if(aux$yanyrowiszeron > 0){
      uadj <- uadj[-aux$yzerowhichrows,]
    }
    #compute objective value:
    mS <- matrix(NA,aux$m,aux$m) #m x m covariance matrix
    mS[lower.tri(mS)] <- pars[aux$cov.indx]
    mS[upper.tri(mS)] <- t(mS)[upper.tri(t(mS))]
    diag(mS) <- pars[aux$sigma2u.indx]
    mSlnDet <- as.numeric(determinant(mS, logarithm=TRUE)$modulus)
    mSinv <- solve(mS)
    tuadjmSuadj <- rowSums( (uadj%*%mSinv)*uadj )
    term1 <- -aux$ynonzerorowsn*aux$m*log(2*pi)/2
    term2 <- -aux$ynonzerorowsn*mSlnDet/2
    term3 <- -sum(tuadjmSuadj)/2
    objective.value <- term1 + term2 + term3
    #check objective value:
    if(is.na(objective.value) || abs(objective.value) == Inf){
      objective.value <- aux$objective.penalty
    }
  }else{
    objective.value <- aux$objective.penalty
  } #end if(chk.conds)
  return(objective.value)
} #end mlgarchObjective

#################
mlgarch <- function(y, arch=1, garch=1, xreg=NULL,
  initial.values=NULL, lower=NULL, upper=NULL,
  nlminb.control=list(), vcov=TRUE, objective.penalty=NULL,
  solve.tol=.Machine$double.eps, c.code=TRUE)
{
  #check/change arguments:
  if(is.null(arch)){ arch <- 0 }
  if(is.null(garch)){ garch <- 0 }
  if(arch < garch) stop("garch order cannot be greater than arch order, since estimation is via the varma representation")
  if(arch > 1) stop("Sorry, arch order cannot be greater than 1 in the current version of mlgarch")

  #zoo:
  y <- as.zoo(y)
  y <- na.trim(y)
  y.index <- index(y)
  y <- coredata(y)
  y.colnames <- colnames(y)
  colnames(y) <- NULL

  #xreg:
  if(!is.null(xreg)){
    if(NROW(xreg)!=NROW(y)) stop("NROW(xreg) must equal NROW(y)")
    xreg <- as.matrix(coredata(as.zoo(xreg)))
    xreg.colnames <- colnames(xreg)
    colnames(xreg) <- NULL
  }

  #rows and dimensions:
  aux <- list()
  aux$y <- y; #rm(y) in the future?
  aux$y.index <- y.index
  aux$n <- NROW(y)
  aux$m <- NCOL(y)
  if(!is.null(xreg)){
    ##OLD (w/bug detected by Rik Wienke) 
    #aux$xreg.m <- NCOL(aux$xreg)
    #aux$xreg <- xreg; #rm(xreg) in the future?
    ##NEW corrected (switch of order):
    aux$xreg <- xreg; #rm(xreg) in the future?
    aux$xreg.m <- NCOL(aux$xreg)
  }

  #orders:
  aux$maxpq <- max(arch,garch)
  aux$nmaxpq <- aux$n + aux$maxpq
  aux$ar <- aux$maxpq
  aux$ma <- garch

  #zeros and more:
  aux$yiszero <- matrix(as.numeric(y==0),aux$n,aux$m)
  aux$yiszeroadj <- rbind(matrix(0,aux$maxpq,aux$m),aux$yiszero)
  aux$yanyrowiszero <- as.numeric(rowSums(aux$yiszero)>0)
  aux$yanyrowiszeroadj <- as.numeric(rowSums(aux$yiszeroadj)>0)
  aux$yzerowhichrows <- which(aux$yanyrowiszero==1)
  aux$yanyrowiszeron <- sum(aux$yanyrowiszero)
  aux$ynonzerorowsn <- aux$n - aux$yanyrowiszeron
  aux$yzeron <- colSums(aux$yiszero)
  aux$ynonzeron <- aux$n - aux$yzeron
  aux$yzerowhere <- list()
  for(i in 1:aux$m){
    if(aux$yzeron[i] > 0){
      aux$yzerowhere[[i]] <- which(y[,i] == 0)
    }else{ aux$yzerowhere[[i]] <- numeric(0) }
  } #end aux$zerowhere list
  y2 <- y^2
  miny2 <- rep(NA,aux$m)
  for(i in 1:aux$m){
    if(aux$yzeron[i] > 0){
      miny2[i] <- min(y2[-aux$yzerowhere[[i]],i])
      y2[aux$yzerowhere[[i]],i] <- miny2[i]
    }
  }
  aux$lny2 <- log(y2)
  aux$Elny2 <- colMeans(aux$lny2)
  for(i in 1:aux$m){
    if(aux$yzeron[i] > 0){
      aux$Elny2[i] <- mean(aux$lny2[-aux$yzerowhere[[i]],i])
    }
  }
  aux$lny2adj <- aux$lny2
  if(aux$maxpq > 0){
    aux$lny2adj <- rbind(t(matrix(rep(aux$Elny2,aux$maxpq),aux$m,aux$maxpq)),
      aux$lny2adj)
#    for(i in 1:aux$maxpq){
#      aux$lny2adj <- rbind(aux$Elny2,aux$lny2adj)
#    }
  }
  rm(y2) #remove y2 for memory-efficiency reasons

  #indices:
  aux$const.indx <- 1:aux$m
  if(aux$ar>0){
    aux$ar.indx <- c(aux$m+1):c(aux$m+aux$ar*aux$m^2)
  }else{
    aux$ar.indx <- 0
  }
  if(aux$ma>0){
    aux$ma.indx <- c(max(aux$ar.indx)+1):c(max(aux$ar.indx)+aux$ma*aux$m^2)
  }else{
    aux$ma.indx <- 0
  }
  if(is.null(aux$xreg)){
    aux$xreg.k <- 0
    aux$xreg.indx <- 0
  }else{
    aux$xreg.k <- aux$xreg.m * aux$m
    aux$xreg.indx <- c(max(aux$m,aux$ar.indx,aux$ma.indx)+1):c(max(aux$m,aux$ar.indx,aux$ma.indx)+aux$xreg.k)
  }
  aux$sigma2u.indx <- c(max(aux$m,aux$ar.indx,aux$ma.indx,aux$xreg.indx)+1):c(max(aux$m,aux$ar.indx,aux$ma.indx,aux$xreg.indx)+aux$m)
  if(aux$m>1){
    aux$cov.indx <- c(max(aux$sigma2u.indx)+1):c(max(aux$sigma2u.indx)+(aux$m^2-aux$m)/2)
  }else{
    aux$cov.k <- 0
    aux$cov.indx <- 0
  }

  #initial values:
  if(is.null(initial.values)){
    if(aux$ma > 0){
      ma.initvals <- as.vector( diag(rep(-0.8/aux$ma, aux$m)) )
    }else{ ma.initvals <- NULL }
    if(aux$ar > 0){
      ar.initvals <- as.vector( diag(rep(0.9/aux$ar, aux$m)) )
    }else{ ar.initvals <- NULL }
    #future: check ar.initvals for stability?
    if(is.null(aux$xreg)){
      xreg.initvals <- numeric(0)
      aux$xregMeans <- 0
    }else{
      xreg.initvals <- rep(0.01, aux$xreg.k)
      aux$xregMeans <- as.numeric(colMeans(aux$xreg))
    }
    IminPhi1 <- diag(rep(1,aux$m))
    if(aux$ar > 0){
      IminPhi1 <- IminPhi1 - matrix(ar.initvals,aux$m,aux$m)
    }
    const.initvals <- IminPhi1%*%cbind(aux$Elny2)
    if(aux$xreg.k > 0){
      const.initvals <- const.initvals - matrix(xreg.initvals,aux$m,aux$xreg.m) %*% aux$xregMeans
    }
    const.initvals <- as.vector(const.initvals)
    sigma2u.initvals <- rep(4.94,aux$m)
    if(aux$m>1){
      cov.initvals <- var(aux$lny2)[lower.tri(var(aux$lny2))]
    }
    initial.values <- c(const.initvals, ar.initvals, ma.initvals,
      xreg.initvals, sigma2u.initvals, cov.initvals)
  }else{
    if( length(initial.values)!=max( aux$cov.indx ) ){
      stop("length(initial.values) not equal to no. of parameters to be estimated")
    } #end check length(initial.values)
  } #end if..else..
  aux$initial.values <- initial.values

  #upper bounds:
  if(is.null(upper)){
    upper <- rep(Inf,aux$m) #constant
    upper <- c(upper, rep(1-.Machine$double.eps,aux$ar*aux$m^2)) #ar parameters
    upper <- c(upper, rep(1-.Machine$double.eps,aux$ma*aux$m^2)) #ma parameters
    upper <- c(upper, rep(Inf,aux$xreg.k))
    upper <- c(upper, rep(Inf,aux$m)) #variances
    upper <- c(upper, rep(Inf,(aux$m^2-aux$m)/2)) #covariances
  }else{
    if( length(upper)!=length(initial.values) )
      stop("length(upper) not equal to length(initial.values)")
  }
  aux$upper <- upper

  #lower bounds:
  if(is.null(lower)){
    lower <- rep(-Inf,aux$m) #constant
    lower <- c(lower, rep(-1+.Machine$double.eps,aux$ar*aux$m^2)) #ar parameters
    lower <- c(lower, rep(-1+.Machine$double.eps,aux$ma*aux$m^2)) #ma parameters
    lower <- c(lower, rep(-Inf,aux$xreg.k))
    lower <- c(lower,rep(0,aux$m)) #variances
    lower <- c(lower, rep(-Inf,(aux$m^2-aux$m)/2)) #covariances
  }else{
    if( length(lower)!=length(initial.values) )
      stop("length(lower) not equal to length(initial.values)")
  }
  aux$lower <- lower

  #misc:
  aux$c.code <- c.code
  aux$solve.tol <- solve.tol
  aux$verboseRecursion <- FALSE
#  aux$yzeroadj <- c(rep(1,max(1,aux$maxpq)), aux$yzero)
#  aux$zerosaux <- rep(0,max(aux$nmaxpq, aux$n+1))
  if(is.null(objective.penalty)){
    aux$objective.penalty <- mlgarchObjective(initial.values, aux)
  }else{
    aux$objective.penalty <- objective.penalty
  }

  #estimate:
  objective.f <- function(pars, x=aux){ -mlgarchObjective(pars,x) }
  est <- nlminb(initial.values, objective.f, lower=lower,
    upper=upper, control=nlminb.control)
  est$objective <- -est$objective
  names(est)[2] <- "objective.varma"

  #parameters:
  uadj <- mlgarchRecursion1(as.numeric(est$par), aux)
  Elnz2 <- rep(NA,aux$m)
  for(i in 1:aux$m){
    if(aux$yzeron[i] > 0){
      uadjtmp <- uadj[-aux$yzerowhere[[i]],i]
    }else{
      uadjtmp <- uadj[,i]
    }
    Elnz2[i] <- -log(mean(exp(uadjtmp - mean(uadjtmp))))
  }
  rm(uadjtmp) #remove object
  parMlgarch <- Elnz2
  namesMlgarch <- paste("Elnz2no",1:aux$m,sep="")
  parVarma <- est$par
  matIndices <- matrix(NA,aux$m,aux$m)
  for(j in 1:aux$m){
    for(i in 1:aux$m){
      matIndices[i,j] <- paste(i,j,sep="")
    }
  }
  namesVarma <- matIndices[lower.tri(matIndices)]
  namesVarma <- c(diag(matIndices),namesVarma)
  namesVarma <- paste("cov",namesVarma,sep="")
#  parMlgarch <- c(parVarma[aux$cov.indx],parMlgarch)
#  parMlgarch <- c(parVarma[aux$sigma2u.indx],parMlgarch)
#  namesMlgarch <- c(namesVarma,namesMlgarch)
  if(aux$xreg.k > 0){
    xregIndices <- matrix(NA,aux$m,aux$xreg.m)
    for(j in 1:aux$xreg.m){
      for(i in 1:aux$m){
#        xregIndices[i,j] <- paste(i,j,sep="")
        xregIndices[i,j] <- paste(i,"no",j,sep="")
      }
    }
    xregNames <- paste("xreg",as.character(xregIndices),sep="")
    namesVarma <- c(xregNames, namesVarma)
    namesMlgarch <- c(xregNames, namesMlgarch)
    parMlgarch <- c(est$par[aux$xreg.indx], parMlgarch)
  }
  if(aux$ma > 0){
    namesVarma <- c(paste("ma",as.character(matIndices),
      ".1",sep=""), namesVarma)
    namesMlgarch <- c(paste("garch",as.character(matIndices),
      ".1",sep=""), namesMlgarch)
    parMlgarch <- c(-parVarma[aux$ma.indx], parMlgarch)
  }
  if(aux$ar > 0){
    namesVarma <- c(paste("ar",as.character(matIndices),
      ".1",sep=""), namesVarma)
    namesMlgarch <- c(paste("arch",as.character(matIndices),
      ".1",sep=""), namesMlgarch)
    if(aux$ma > 0){
      parTmp <- c(parVarma[aux$ar.indx]+parVarma[aux$ma.indx])
    }else{
      parTmp <- parVarma[aux$ar.indx]
    }
    parMlgarch <- c(parTmp, parMlgarch)
  }
  namesVarma <- c(paste("intercept",1:aux$m,sep=""),namesVarma)
  namesMlgarch <- c(paste("intercept",1:aux$m,sep=""),namesMlgarch)
  if(aux$ma > 0){
    tmpMlgarch <- parVarma[aux$const.indx] - (diag(rep(1,aux$m))+matrix(parVarma[aux$ma.indx],aux$m,aux$m))%*%Elnz2
  }else{
    tmpMlgarch <- parVarma[aux$const.indx] - diag(rep(1,aux$m))%*%Elnz2
  }
  parMlgarch <- c(tmpMlgarch, parMlgarch)
  names(parVarma) <- namesVarma
  names(parMlgarch) <- namesMlgarch
  est$par <- parMlgarch
  est <- c(list(date=date(), par.varma=parVarma), est)

  #vcov matrix:
  if(vcov){
    hessian.varma <- -optimHess(as.numeric(parVarma), objective.f)
    colnames(hessian.varma) <- namesVarma
    rownames(hessian.varma) <- namesVarma
    vcov.varma <- solve(-hessian.varma, tol=aux$solve.tol)
    est <- c(list(aux=aux, hessian.varma=hessian.varma,
      vcov.varma=vcov.varma, vcov.mlgarch=NULL), est)
    est$vcov.mlgarch <- vcov.mlgarch(est)
  }

  #out:
  if(is.null(est$aux)){
    est <- c(list(aux=aux),est)
  }
  class(est) <- "mlgarch"
  return(est)
} #end mlgarch
    
    
###########################################################
## PART 2: S3 METHODS FOR MLGARCH CLASS
###########################################################

##################
coef.mlgarch <- function(object, varma=FALSE, ...)
{
  if(varma){
    result <- object$par.varma
  }else{
    result <- object$par
  }
  return(result)
} #end coef.mlgarch

################### fitted values
fitted.mlgarch <- function(object, varma=FALSE, verbose=FALSE, ...)
{
  aux <- object$aux
  aux$verboseRecursion <- TRUE
  pars <- as.numeric(object$par.varma)
  mUhat <- mlgarchRecursion1(pars, aux)
  if(varma){
    result <- mUhat[,c(aux$m+1):c(2*aux$m)]
  }else{
    logsigma2 <- matrix(NA,aux$n,aux$m)
    sigma <- matrix(NA,aux$n,aux$m)
    i2 <- length(object$par)
    i1 <- i2 - aux$m + 1
    Elnz2 <- as.numeric(object$par[i1:i2])
    for(i in 1:aux$m){
      uhatnotzero <- as.numeric(mUhat[,i]!=0)
      lnz2adj <- mUhat[,i] + uhatnotzero*Elnz2[i]
      logsigma2[,i] <- mUhat[,c(i+aux$m)] - lnz2adj
      sigma[,i] <- exp(logsigma2[,i]/2)
    }
#    colnames(sigma) <- paste("sigmaFit",seq(1:aux$m),sep="")
    colnames(sigma) <- paste("sd",seq(1:aux$m),sep="")
    result <- sigma
    if(verbose){
#      colnames(logsigma2) <- paste("logsigma2Fit",seq(1,aux$m),sep="")
      colnames(logsigma2) <- paste("lnsd2no",seq(1,aux$m),sep="")
      result <- cbind(result,logsigma2)
    } #end if(verbose)
  } #end if(varma)else(..)
  result <- zoo(result, order.by=aux$y.index)
  return(result)
} #end fitted.mlgarch

###################### check if object is of class lgarch
##is.mlgarch <- function(x)
##{
##  if(class(x)=="mlgarch"){ TRUE }else{ FALSE}
##}

#################### return log-likelihood
logLik.mlgarch <- function(object, varma=FALSE, ...)
{
  if(varma==TRUE){
    result <- object$objective.varma
    attr(result, "df") <- length(object$par.varma)
  }else{
    mZhat <- residuals.mlgarch(object)
    mSigmaFit <- fitted.mlgarch(object)
    if(object$aux$yanyrowiszeron > 0){
      mZhat <- mZhat[-object$aux$yzerowhichrows,]
      mSigmaFit <- mSigmaFit[-object$aux$yzerowhichrows,]
    }
    mR <- cor(mZhat)
    mRlnDet <- as.numeric(determinant(mR, logarithm=TRUE)$modulus)
    mRinv <- solve(mR)
    tmZhatmRmZhat <- rowSums( (mZhat%*%mRinv)*mZhat )
    term1 <- -object$aux$ynonzerorowsn*(object$aux$m*log(2*pi) + mRlnDet)/2
    term2 <- -sum(rowSums(log(mSigmaFit)))
    term3 <- -sum(tmZhatmRmZhat)/2
    result <- term1 + term2 + term3
    attr(result, "df") <- length(object$par) - object$aux$m
  }
#  attr(result, "nobs") <- length(object$aux$n)
  attr(result, "nobs") <- length(object$aux$ynonzerorowsn)
  class(result) <- "logLik"
  return(result)
} #end logLik.mlgarch

################### forecast
##predict.mlgarch <- function(object, n.ahead=1,
##  initial.values=NULL, n.sim=10000, verbose=FALSE, ...)
##{
##return(out)
##} #end predict.lgarch

################## summarise output
print.mlgarch <- function(x, varma=FALSE, ...)
{
  pars <- coef.mlgarch(x, varma=varma)
  vcovmat <- vcov.mlgarch(x, varma=varma)
  out1 <- rbind(pars, sqrt(diag(vcovmat)))
  rownames(out1) <- c("Estimate:", "Std. Error:")
  out2 <- as.data.frame(matrix(NA,4,1))
  out2[1,1] <- as.character(round(logLik.mlgarch(x, varma=FALSE), digits=3))
  out2[2,1] <- as.character(round(logLik.mlgarch(x, varma=TRUE), digits=3))
  out2[3,1] <- as.character(round(x$aux$ynonzerorowsn, digits=0))
  out2[4,1] <- as.character(round(x$aux$yanyrowiszeron, digits=0))
  rownames(out2) <-   c("Log-likelihood (log-mgarch):",
    "Log-likelihood (varma):", "No. of obs. without zeros:",
    "No. of obs. with zeros:")
  colnames(out2) <- ""
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Method: Multivariate ML \n")
  cat("Message (nlminb):", x$message, "\n")
  cat("No. of observations:", x$aux$n, "\n")
  cat("Sample:", as.character(x$aux$y.index[1]),
    "to", as.character(x$aux$y.index[x$aux$n]), "\n")
  cat("\n")
  cat("Coefficients:\n")
  cat("\n")
  print(out1)
  print(out2)
  cat("\n")
} #end print.mlgarch

###################
residuals.mlgarch <- function(object, varma=FALSE, ...)
{
  aux <- object$aux
  aux$verboseRecursion <- TRUE
  pars <- as.numeric(object$par.varma)
  mUhat <- mlgarchRecursion1(pars, aux)
  if(varma){
    result <- mUhat[,1:aux$m]
  }else{
    result <- matrix(NA,aux$n,aux$m)
    i2 <- length(object$par)
    i1 <- i2 - aux$m + 1
    Elnz2 <- as.numeric(object$par[i1:i2])
    for(i in 1:aux$m){
      uhatnotzero <- as.numeric(mUhat[,i]!=0)
      lnz2adj <- mUhat[,i] + uhatnotzero*Elnz2[i]
      logsigma2 <- mUhat[,c(i+aux$m)] - lnz2adj
      sigma <- exp(logsigma2/2)
      result[,i] <- aux$y[,i]/sigma
    }
    colnames(result) <- paste("z",seq(1,aux$m),sep="")
  }
  result <- zoo(result, order.by=aux$y.index)
  return(result)
} #end residuals.mlgarch

################## summarise output
summary.mlgarch <- function(object, ...)
{
  xnames <- cbind(names(object))
  xrows <- nrow(xnames)
  colnames(xnames) <- ""
  rownames(xnames) <- rep(" ",xrows)
  cat("\n")
  cat("Items in list:")
  cat("\n")
  print(xnames)
  cat("\n")
} #end summary.mlgarch

################### return variance-covariance matrix
vcov.mlgarch <- function(object, varma=FALSE, ...)
{
  #if varma:
  if(varma){
    if(is.null(object$vcov.varma)){
      cat("vcov cannot be extracted (vcov=FALSE during estimation)\n")
      m.dim <- length(object$par.varma)
      result <- matrix(NA,m.dim,m.dim)
      colnames(result) <- names(object$par.varma)
      rownames(result) <- names(object$par.varma)
    }else{ result <- object$vcov.varma }
  }

  #if mlgarch:
  if(varma==FALSE){
    if(!is.null(object$vcov.mlgarch)){
      result <- object$vcov.mlgarch
    }else{

      #extract vcov.varma:
      if(is.null(object$vcov.varma)){
        cat("vcov cannot be extracted (vcov=FALSE during estimation)\n")
        m.dim <- length(object$par)
        result <- matrix(NA,m.dim,m.dim)
        colnames(result) <- names(object$par)
        rownames(result) <- names(object$par)
      }else{
        vcov.varma <- object$vcov.varma

        #make vcov.mlgarch:
        vcov.mlgarch <- matrix(NA,length(object$par),length(object$par))
        colnames(vcov.mlgarch) <- names(object$par)
        rownames(vcov.mlgarch) <- names(object$par)

        #make garch part:
        if(object$aux$ma > 0){
        vcov.mlgarch[object$aux$ma.indx, object$aux$ma.indx] <- vcov.varma[object$aux$ma.indx, object$aux$ma.indx]
        } #end if(..ma > 0)

        #make arch part:
        if(object$aux$ar > 0){
          ar.ma.diff <- object$aux$ar - object$aux$ma
          if(ar.ma.diff==0){
            vcov.ar <- diag(vcov.varma[object$aux$ar.indx,object$aux$ar.indx])
            vcov.ma <- diag(vcov.varma[object$aux$ma.indx,object$aux$ma.indx])
            vcov.ar.ma <- diag(vcov.varma[object$aux$ar.indx,object$aux$ma.indx])
            vcov.tmp <- matrix(NA,length(vcov.ar),length(vcov.ar))
            diag(vcov.tmp) <- vcov.ar + vcov.ma + 2*vcov.ar.ma
            vcov.mlgarch[object$aux$ar.indx,object$aux$ar.indx] <- vcov.tmp
          }else{
            vcov.mlgarch[object$aux$ar.indx,object$aux$ar.indx] <- vcov.varma[object$aux$ar.indx,object$aux$ar.indx]
          } #end if(..diff==0)
        } #end if(..ar > 0)

        #make xreg part:
        if(object$aux$xreg.k > 0){
          vcov.mlgarch[object$aux$xreg.indx, object$aux$xreg.indx] <- vcov.varma[object$aux$xreg.indx, object$aux$xreg.indx]
#           #univariate code:
#         if(object$aux$ma > 0){
#           vcov.lgarch[object$aux$ma.indx, object$aux$xreg.indx] <- -vcov.arma[object$aux$ma.indx,object$aux$xreg.indx]
#           vcov.lgarch[object$aux$xreg.indx, object$aux$ma.indx] <- -vcov.arma[object$aux$xreg.indx, object$aux$ma.indx]
#         }
#         if(object$aux$ar > 0){
#           vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcov.arma[object$aux$ar.indx,object$aux$xreg.indx]
#           if(object$aux$ma > 0){
#             vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] + vcov.lgarch[object$aux$ma.indx, object$aux$xreg.indx]
#           }
#           vcov.lgarch[object$aux$xreg.indx, object$aux$ar.indx] <- vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx]
#         } #end if(..ar > 0)
        } #end if(..k > 0)

        #Var(Elnz2^hat):
        uadj <- mlgarchRecursion1(as.numeric(object$par.varma), object$aux)
        uadj.m <- NCOL(uadj)
        vcov.m <- NROW(vcov.mlgarch)
        for(i in 1:uadj.m){
          uadjtmp <- uadj[,i]
          if(object$aux$yzeron[i] > 0){
            uadjtmp <- uadjtmp[ -object$aux$yzerowhere[[i]] ]
          }
          expuadj <- exp(uadjtmp)
          uexpuadj <- uadjtmp*exp(uadjtmp)
          avaruadj <- var(expuadj)/mean(expuadj)^2 + var(uadjtmp) - 2*mean(uexpuadj)/mean(expuadj)
          dim1 <- vcov.m - uadj.m + i
          dim2 <- dim1
          vcov.mlgarch[dim1, dim2] <- avaruadj/length(uadjtmp)
        } #end for(i in..)

      result <- vcov.mlgarch
      } #end if(is.null(vcov.varma))else(..)
    } #end if(!is.null(vcov.mlgarch))else(..)
  } #end if(varma==FALSE)else(..)

  #out:
  return(result)
} #end vcov.mlgarch
