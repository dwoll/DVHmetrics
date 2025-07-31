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
Version: \tab 0.4.3\cr
Date: \tab 2025-07-31\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Daniel Wollschlaeger and Heiko Karle

Maintainer: Daniel Wollschlaeger <wollschlaeger@uni-mainz.de>
}
\references{
For solutions that also read files in DICOM format, see packages \code{espadon} (\url{https://espadon.cnrs.fr/}, \url{https://CRAN.R-project.org/package=espadon}) and \code{RadOnc} (\url{https://CRAN.R-project.org/package=RadOnc}).
}
\keyword{package}
\examples{
showDVH(dataMZ[[1]])
checkConstraint(dataMZ, "D1CC < 10Gy")
}
