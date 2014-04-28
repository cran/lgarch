print.lgarch <-
function(x, arma.version=FALSE, ...)
{
  pars <- coef.lgarch(x, arma=arma.version)
  vcovmat <- vcov.lgarch(x, arma.version=arma.version,
    full.matrix=TRUE)
  if(is.null(vcovmat)){
    vcovmat <- matrix(NA, length(pars), length(pars))
    colnames(vcovmat) <- names(pars)
    rownames(vcovmat) <- names(pars)
  }
  out1 <- rbind(pars, sqrt(diag(vcovmat)))
  rownames(out1) <- c("Estimate:", "Std. Error:")
  out2 <- as.data.frame(matrix(NA,4,1))
  out2[1,1] <- as.character(round(logLik.lgarch(x, arma=FALSE), digits=3))
  out2[2,1] <- as.character(round(logLik.lgarch(x, arma=TRUE), digits=3))
  out2[3,1] <- as.character(round(rss.lgarch(x), digits=3))
  out2[4,1] <- as.character(round(x$aux$ynonzeron, digits=0))
  out2[5,1] <- as.character(round(x$aux$yzeron, digits=0))
  rownames(out2) <-   c("Log-likelihood (log-garch):",
    "Log-likelihood (arma):", "Sum Squared Resids. (arma):",
    "No. of non-zeros:", "No. of zeros:")
  colnames(out2) <- ""
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Message (nlminb):", x$message, "\n")
  cat("No. of observations:", x$aux$n, "\n")
#  cat("Sample:", object$aux$y.index[1], "to",
#    object$aux$y.index[object$aux$n], "\n")
  cat("\n")
  cat("Coefficients:\n")
  print(out1)
  print(out2)
  cat("\n")
}
