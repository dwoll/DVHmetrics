\name{getEQD2}
\alias{getEQD2}
\alias{getEQD2.default}
\alias{getEQD2.DVHs}
\alias{getEQD2.DVHLst}
\alias{getEQD2.DVHLstLst}
\title{2Gy fractions biologically equivalent dose (EQD2)}
\description{
Calculate dose in 2Gy fractions biologically equivalent dose according to the linear-quadratic model.
}
\usage{
getEQD2(D=NULL, fd=NULL, fn=NULL, ab=NULL)

\method{getEQD2}{default}(D=NULL, fd=NULL, fn=NULL, ab=NULL)

\method{getEQD2}{DVHs}(D=NULL, fd=NULL, fn=NULL, ab=NULL)

\method{getEQD2}{DVHLst}(D=NULL, fd=NULL, fn=NULL, ab=NULL)

\method{getEQD2}{DVHLstLst}(D=NULL, fd=NULL, fn=NULL, ab=NULL)
}
\arguments{
  \item{D}{Default: Total dose. If \code{NULL}, \code{fn} must be given. Alternative: One cumulative DVH (object of class \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{fd}{Fractional dose. If \code{D} is some kind of \code{DVH} object, only the first element will be used.}
  \item{fn}{Number of fractions. If \code{NULL}, \code{D} must be given. Ignored if \code{D} is some kind of \code{DVH} object.}
  \item{ab}{alpha/beta ratio for the relevant tissue.  If \code{D} is some kind of \code{DVH} object, only the first element will be used.}
}
\details{
EQD2 is a special case of isoeffective dose calculation with fractional dose \code{d2=2}, see \code{\link[DVHmetrics]{getIsoEffD}}.
}
\value{
Default method: A data frame with variables \code{EQD2}, \code{fractDose}, \code{ab}.

If \code{D} is some kind of \code{DVH} object, the same kind of object is returned with the individual dose values converted to EQD2.
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
getEQD2(dataMZ[[c(1, 1)]], fd=1.8, ab=3)
}
