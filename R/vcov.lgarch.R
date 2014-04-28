vcov.lgarch <-
function(object, arma.version=FALSE,
  full.matrix=FALSE, ...)
{
  #compute vcov:
  if(is.null(object$vcov.arma)){
    result <- NULL
  }else{
    vcov.arma <- object$vcov.arma
    if(arma.version){
      result <- vcov.arma
    }else{

      #make vcov.lgarch matrix:
      vcov.lgarch <- matrix(NA,length(object$par),
        length(object$par))
      colnames(vcov.lgarch) <- names(object$par)
      rownames(vcov.lgarch) <- names(object$par)

      if(is.null(object$aux)){
        cat(" Auxiliary list 'aux' not stored, vcov cannot be extracted")
      }else{

        if(object$aux$ma > 0){
          vcov.lgarch[object$aux$ma.indx, object$aux$ma.indx] <- vcov.arma[object$aux$ma.indx, object$aux$ma.indx]
        } #end if(..ma > 0)

        if(object$aux$ar > 0){
          ar.ma.diff <- object$aux$ar - object$aux$ma
          if(ar.ma.diff==0){
            varalpha <- vcov.arma[2,2] + vcov.arma[3,3] + 2*vcov.arma[2,3]
            vcov.lgarch[2,2] <- varalpha
            covalphabeta <- -vcov.arma[2,3] - vcov.arma[3,3]
            vcov.lgarch[2,3] <- covalphabeta
            vcov.lgarch[3,2] <- covalphabeta
          }else{
            vcov.lgarch[2,2] <- vcov.arma[2,2]
          } #end if(..diff==0)
        } #end if(..ar > 0)

        if(object$aux$xreg.k > 0){
          vcov.lgarch[object$aux$xreg.indx, object$aux$xreg.indx] <- vcov.arma[object$aux$xreg.indx, object$aux$xreg.indx]
          if(object$aux$ma > 0){
            vcov.lgarch[object$aux$ma.indx, object$aux$xreg.indx] <- -vcov.arma[object$aux$ma.indx,object$aux$xreg.indx]
            vcov.lgarch[object$aux$xreg.indx, object$aux$ma.indx] <- -vcov.arma[object$aux$xreg.indx, object$aux$ma.indx]
          }
          if(object$aux$ar > 0){
            vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcov.arma[object$aux$ar.indx,object$aux$xreg.indx]
            if(object$aux$ma > 0){
              vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] <- vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx] + vcov.lgarch[object$aux$ma.indx, object$aux$xreg.indx]
            }
            vcov.lgarch[object$aux$xreg.indx, object$aux$ar.indx] <- vcov.lgarch[object$aux$ar.indx, object$aux$xreg.indx]
          } #end if(..ar > 0)

        } #end if(..k > 0)

        #Var(Elnz2^hat):
        uadj <- lgarchRecursion1(as.numeric(object$par.arma), object$aux)
        if(object$aux$yzeron > 0){
          uadj <- uadj[-object$aux$yzerowhere]
        }  #end if..
        expuadj <- exp(uadj)
        uexpuadj <- uadj*exp(uadj)
        avaruadj <- var(expuadj)/mean(expuadj)^2 + var(uadj) - 2*mean(uexpuadj)/mean(expuadj)
        vcov.lgarch[NROW(vcov.lgarch), NCOL(vcov.lgarch)] <- avaruadj/length(uadj)

      } #end if(is.null(object$aux))else(..)

      #out:
      if(full.matrix){
        result <- vcov.lgarch
      }else{
        vcov.lgarch <- vcov.lgarch[-NROW(vcov.lgarch),-NCOL(vcov.lgarch)]
        vcov.lgarch <- vcov.lgarch[-1,-1]
        result <- vcov.lgarch
      }
    } #end if(arma.version)
  } #end if(is.null(vcov.arma)))
  return(result)
}
