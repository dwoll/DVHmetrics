\name{getIsoEffD}
\alias{getIsoEffD}
\alias{getIsoEffD.default}
\alias{getIsoEffD.DVHs}
\alias{getIsoEffD.DVHLst}
\alias{getIsoEffD.DVHLstLst}
\title{Isoeffective dose calculation}
\description{
Convert given (fractional) dose into a corresponding (fractional) dose for a different total dose / fractionation schedule according to the linear-quadratic model.
}
\usage{
getIsoEffD(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL)

\method{getIsoEffD}{default}(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL)

\method{getIsoEffD}{DVHs}(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL)

\method{getIsoEffD}{DVHLst}(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL)

\method{getIsoEffD}{DVHLstLst}(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL)
}
\arguments{
  \item{D1}{Default: \code{numeric} vector. Total dose 1. Alternative: One cumulative DVH (object of class \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{D2}{\code{numeric} vector. Total dose 2. Ignored if \code{D} is some kind of \code{DVH} object.}
  \item{fd1}{\code{numeric} vector. Fractional dose 1. If \code{D} is some kind of \code{DVH} object, only the first element will be used.}
  \item{fd2}{\code{numeric} vector. Fractional dose 2. If \code{D} is some kind of \code{DVH} object, only the first element will be used.}
  \item{ab}{\code{numeric} vector. alpha/beta ratio for the relevant tissue in the linear-quadatic model. If \code{D} is some kind of \code{DVH} object, only the first element will be used.}
}
\details{
\code{DVH} methods: Calculate \code{D2} based on \code{D1}, \code{fd1}, \code{fd2}, and \code{ab}. The default method can also calculate \code{fd2} based on \code{D1}, \code{D2}, \code{fd1}, and \code{ab}.
}
\value{
The (vector of) corresponding (fractional) dose value(s). If \code{D} is some kind of \code{DVH} object, the same kind of object is returned with the individual dose values converted to \code{D2}.
}
\references{
IAEA, & ICRU. (2008). Relative biological effectiveness in ion-beam therapy (Tech. Rep. No. IAEA-TR 461). Vienna, Austria: IAEA (International Atomic Energy Agency) and ICRU (International
Commission on Radiation Units and Measurements).
}
\seealso{
\code{\link[DVHmetrics]{getBED}},
\code{\link[DVHmetrics]{getEQD2}}
}
\examples{
# reference: 70Gy in 2Gy fractions
# new fractionation: 3Gy fractions
# calculate corresponding dose
(D2 <- getIsoEffD(D1=70, fd1=2, fd2=3, ab=c(3.5, 10)))

getIsoEffD(D1=dataMZ[[c(1, 1)]], fd1=1.8, fd2=2, ab=3.5)
}
