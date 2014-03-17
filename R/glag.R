glag <-
function(x, k=1, pad=FALSE, pad.value=NA){
  x.n <- length(x)
  x.nmink <- x.n - k
  xlagged <- x[1:x.nmink]
  if(pad){
    xlagged <- c(rep(pad.value,k), xlagged)
  }
  return(xlagged)
}
