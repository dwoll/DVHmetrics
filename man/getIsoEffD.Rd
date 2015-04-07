\name{getIsoEffD}
\alias{getIsoEffD}
\title{2Gy fractions biologically equivalent dose (EQD2)}
\description{
Convert given (fractional) dose into an isoeffective (fractional) dose for a different total dose / fractionation schedule according to the linear-quadratic model.
}
\usage{
getIsoEffD(D1=NULL,  D2=NULL,
           fd1=NULL, fd2=NULL,
           fn1=NULL, fn2=NULL, ab=NULL)
}
\arguments{
  \item{D1}{\code{numeric} vector. Total dose 1.}
  \item{D2}{\code{numeric} vector. Total dose 2.}
  \item{fd1}{\code{numeric} vector. Fractional dose 1.}
  \item{fd2}{\code{numeric} vector. Fractional dose 2.}
  \item{fn1}{\code{numeric} vector. Number of fractions 1.}
  \item{fn2}{\code{numeric} vector. Number of fractions 2.}
  \item{ab}{\code{numeric} vector. alpha/beta ratio for the relevant tissue in the linear-quadatic model.}
}
\details{
Either calculate \code{D2} based on \code{D1}, \code{fd1}, \code{fd2}, and \code{ab}. Or calculate \code{fd2} based on \code{D1}, \code{D2}, \code{fd1}, and \code{ab}.

Instead of providing \code{D1} (or \code{D2}) directly, \code{fd1} AND \code{fn1} (or \code{fd2} AND \code{fn2}) can be given.
}
\value{
The (vector of) isoeffective (fractional) dose value(s)
}
\references{
Barendsen, G. W. (1982) Dose fractionation, dose-rate and iso-effect relationships for normal-tissue response. International Journal of Radiation Oncology Biology Physics 8(11), 1981-1997.
}
\seealso{
\code{\link[DVHmetrics]{getBED}},
\code{\link[DVHmetrics]{getEQD2}}
}
\examples{
getIsoEffD(D1=50, fd1=2.5, fd2=3, ab=c(2, 3, 4))
}
