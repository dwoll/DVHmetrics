\name{getBED}
\alias{getBED}
\title{Calculate biologically effective dose (BED)}
\description{
Calculate biologically effective dose (BED) according to the linear-quadratic model.
}
\usage{
getBED(D=NULL, fd=NULL, fn=NULL, ab=NULL)
}
\arguments{
  \item{D}{Total dose. If \code{NULL}, \code{nf} must be given.}
  \item{fd}{Fractional dose.}
  \item{fn}{Number of fractions. If \code{NULL}, \code{D} must be given.}
  \item{ab}{alpha/beta ratio for the relevant tissue.}
}
\value{
The (vector of) BED value(s)
}
\references{
Fowler, J. F. (2010). 21 years of Biologically Effective Dose. British Journal of Radiology, 83, 554-568.
}
\seealso{
\code{\link[DVHmetrics]{getEQD2}},
\code{\link[DVHmetrics]{getIsoEffD}}
}
\examples{
getBED(D=50, fd=2.5, ab=c(2, 3, 4))
}
