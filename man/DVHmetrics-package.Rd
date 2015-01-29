\name{DVHmetrics-package}
\alias{DVHmetrics-package}
\alias{DVHmetrics}
\docType{package}
\title{Analyze Dose-Volume Histograms and Check Constraints}
\description{
Functionality for analyzing dose-volume histograms
(DVH) in radiation oncology: Read DVH text files,
calculate DVH metrics, show DVH diagrams,
check and visualize quality assurance constraints
for the DVH. Includes web-based graphical user interface.
}
\details{
\tabular{ll}{
Package: \tab DVHmetrics\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-01-29\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Daniel Wollschlaeger and Heiko Karle

Maintainer: Daniel Wollschlaeger <wollschlaeger@uni-mainz.de>
}
\references{
For a solution that also reads files in DICOM-RT format and has a broader scope in supporting radiation oncology, see the \code{RadOnc} package: \url{http://CRAN.R-project.org/package=RadOnc}.
}
\keyword{package}
\examples{
data(dataMZ)
showDVH(dataMZ[[1]])
checkConstraint(dataMZ, "D1CC < 10Gy")
}
