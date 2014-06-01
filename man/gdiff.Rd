\name{gdiff}
\alias{gdiff}

\title{ Difference a vector or a matrix, with special treatment of zoo objects}

\description{
Similar to the \code{\link{diff}} function from the base package, but gdiff enables padding (e.g. NAs or 0s) of the lost entries. Contrary to the \code{\link{diff}} function in the base package, however, the default in gdiff is to pad (with NAs). The gdiff function is particularly suited for zoo objects, since their indexing is retained
}

\usage{
gdiff(x, lag = 1, pad = TRUE, pad.value = NA)
}

\arguments{
  \item{x}{a numeric vector or matrix}
  \item{lag}{integer equal to the difference-length (the default is 1)}
  \item{pad}{logical. If TRUE (default), then the lost entries are padded with pad.value. If FALSE, then no padding is undertaken}
  \item{pad.value}{numeric. The pad-value}
}

\value{A vector or matrix with the differenced values}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\note{Empty}

\seealso{
\code{\link{diff}}, \code{\link{glag}}, \code{\link{lag}}
}

\examples{
##1st difference of a series:
x <- rnorm(5)
gdiff(x)

##1st difference with no padding:
gdiff(x, pad=FALSE)

##1st difference retaining the original zoo-index ordering:
gdiff(as.zoo(x))

##1st difference of a matrix:
y <- matrix(rnorm(8),4,2)
gdiff(y)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
