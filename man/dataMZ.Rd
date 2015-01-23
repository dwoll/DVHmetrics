\name{dataMZ}
\alias{dataMZ}
\docType{data}
\title{DVH data from 3 patients}
\description{
Data from 3 patients with radiotherapy after breast cancer surgery. DVHs for 7 heart structures.}
\usage{data(dataMZ)}
\format{
  Object of class \code{DVHLstLst} with 3 components corresponding to 3 patients.
  \describe{
    \item{\code{P123}}{Object of class \code{DVHLst}. 7 objects of class \code{DVHs} for structures AMYOCL (left anterior heart wall), AMYOCR (right anterior heart wall), AOVALVE (aortic valve), AVNODE (AV node), HEART (complete heart), PULMVALVE (pulmonary valve), MYOCARD (heart wall)}
    \item{\code{P234}}{Object of class \code{DVHLst}. 7 objects of class \code{DVHs} for the same structures as patient \code{P123}.}
    \item{\code{P345}}{Object of class \code{DVHLst}. 7 objects of class \code{DVHs} for the same structures as patient \code{P123}.}
  }
}
\details{
Data courtesy of Heiko Karle, Department of Radiation Oncology (Prof. Dr. Schmidberger), University Medical Center Mainz, Germany.

See \code{\link[DVHmetrics]{readDVH}} for classes \code{DVHLstLst}, \code{DVHLst}, and \code{DVHs}.
}
\seealso{
\code{\link[DVHmetrics]{readDVH}},
\code{\link[DVHmetrics]{print.DVHs}}
}
\examples{
print(dataMZ, verbose=TRUE)
}
\keyword{datasets}
