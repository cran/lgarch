.onAttach <- function(libname, pkgname)
{
txt <- c("\n",
  paste(sQuote("lgarch"), "version 0.7\n"),
  "\n",
  paste0("Simulation and Estimation of Log-GARCH Models"),
  "\n",
  paste("CRAN website: https://CRAN.R-project.org/package=lgarch"),
  "\n")
  
  ##print message:
  if(interactive() || getOption("verbose")){
    packageStartupMessage(paste(strwrap(txt, indent = 2,
      exdent = 4), collapse = "\n"))
  }
} #close .onAttach
