\name{DVHmetrics-package}
\alias{DVHmetrics-package}
\alias{DVHmetrics}
\docType{package}
\title{Analyze Dose-Volume Histograms and Check Constraints}
\description{
Functionality for analyzing dose-volume histograms
(DVH) in radiation oncology: Read DVH text files,
calculate DVH metrics, gEUD, BED, EQD2, NTCP, TCP, show DVH
diagrams, check and visualize quality assurance constraints
for the DVH. Includes web-based graphical user interface.
}
\details{
\tabular{ll}{
Package: \tab DVHmetrics\cr
Type: \tab Package\cr
Version: \tab 0.3.7.0.999\cr
Date: \tab 2018-10-08\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Daniel Wollschlaeger and Heiko Karle

Maintainer: Daniel Wollschlaeger <wollschlaeger@uni-mainz.de>
}
\references{
For a solution that also reads files in DICOM-RT format, see the \code{RadOnc} package: \url{https://CRAN.R-project.org/package=RadOnc}.
}
\keyword{package}
\examples{
showDVH(dataMZ[[1]])
checkConstraint(dataMZ, "D1CC < 10Gy")
}
