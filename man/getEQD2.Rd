\name{getEQD2}
\alias{getEQD2}
\title{2Gy fractions biologically equivalent dose (EQD2)}
\description{
Calculate dose in 2Gy fractions biologically equivalent dose according to the linear-quadratic model.
}
\usage{
getEQD2(D=NULL, fd=NULL, fn=NULL, ab=NULL)
}
\arguments{
  \item{D}{Total dose. If \code{NULL}, \code{nf} must be given.}
  \item{fd}{Fractional dose.}
  \item{fn}{Number of fractions. If \code{NULL}, \code{D} must be given.}
  \item{ab}{alpha/beta ratio for the relevant tissue.}
}
\details{
EQD2 is a special case of isoeffective dose calculation with fractional dose \code{d2=2}, see \code{\link[DVHmetrics]{getIsoEffD}}.
}
\value{
The (vector of) EQD2 value(s)
}
\references{
IAEA, & ICRU. (2008). Relative biological effectiveness in ion-beam therapy (Tech. Rep. No. IAEA-TR 461). Vienna, Austria: IAEA (International Atomic Energy Agency) and ICRU (International
Commission on Radiation Units and Measurements).
}
\seealso{
\code{\link[DVHmetrics]{getBED}},
\code{\link[DVHmetrics]{getIsoEffD}}
}
\examples{
getEQD2(D=50, fd=2.5, ab=c(2, 3, 4))
}
