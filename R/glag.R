glag <-
function(x, k=1, pad=FALSE, pad.value=NA){
  #check arguments:
  if(k < 1) stop("Lag order k cannot be less than 1")

  #zoo-related operations:
  iszoo.chk <- is.zoo(x)
  x <- as.zoo(cbind(x))
  x.n <- NROW(x)
  x.ncol <- NCOL(x)
  x.index <- index(x)
  x <- coredata(x)

  #do the lagging:
  x.nmink <- x.n - k
  xlagged <- matrix(x[1:x.nmink,], x.nmink, x.ncol)
  if(pad){
    xlagged <- rbind( matrix(pad.value,k,x.ncol) , xlagged)
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
  if(x.ncol==1){
    xlagged <- as.vector(xlagged)
  }
  return(xlagged)
}
