\name{glag}
\alias{glag}

\title{ Lag vector or matrix}

\description{
Similar to the \code{\link{lag}} function from the stats package, but glag enables padding (e.g. NAs or 0s) of the lost entries.
}

\usage{
glag(x, k = 1, pad = FALSE, pad.value = NA)
}

\arguments{
  \item{x}{a numeric vector}
  \item{k}{integer equal to the lag (the default is 1)}
  \item{pad}{logical. If TRUE, then the lost entries are padded with pad.value. If FALSE (default), then no padding is undertaken}
  \item{pad.value}{numeric. The padding-value}
}

\value{A vector with the lagged values}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\note{Empty}

\seealso{
\code{\link{lag}}
}

\examples{
##lag series with NA for the missing entries:
x <- rnorm(5)
xlag1 <- glag(x, pad=TRUE)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
