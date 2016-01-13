\name{getDMEAN}
\alias{getDMEAN}
\alias{getDMEAN.DVHs}
\alias{getDMEAN.DVHLst}
\alias{getDMEAN.DVHLstLst}
\title{DMEAN and other dose metrics}
\description{
Calculate DMEAN and other dose metrics from the (interpolated) differential DVH without relying on the values exported by the TPS.
}
\usage{
getDMEAN(x, interp=c("linear", "spline", "ksmooth", "smoothSpl"),
         nodes=5001L)

\method{getDMEAN}{DVHs}(x, interp=c("linear", "spline", "ksmooth", "smoothSpl"),
         nodes=5001L)

\method{getDMEAN}{DVHLst}(x, interp=c("linear", "spline", "ksmooth", "smoothSpl"),
         nodes=5001L)

\method{getDMEAN}{DVHLstLst}(x, interp=c("linear", "spline", "ksmooth", "smoothSpl"),
         nodes=5001L)
}
\arguments{
  \item{x}{One DVH (object of class \code{DVHs}, multiple DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{interp}{\code{character}. Method of interpolation between DVH points: Linear interpolation applies to the cumulative DVH (recommended). Spline interpolation with \code{\link[stats]{splinefun}}, local polynomial regression with \code{\link[KernSmooth]{locpoly}}, and smoothing splines with \code{\link[stats]{smooth.spline}} apply to the differential DVH (not recommended).}
  \item{nodes}{\code{numeric}. Minimum number of nodes to use in interpolation. Number of available nodes is kept as is for \code{NULL} or if larger than \code{nodes}.}
}
\value{
A data frame with the following value(s).
\describe{
\item{patID}{Patient ID.}
\item{structure}{Structure name.}
\item{doseMin}{Minimum dose.}
\item{doseMax}{Maximum dose.}
\item{doseAvg}{Mean dose.}
\item{doseMed}{Median dose.}
\item{doseSD}{Dose standard deviation.}
\item{doseMode}{Dose mode.}
\item{doseAvgTPS}{Mean dose as exported from the TPS (if available).}
\item{doseMedTPS}{Median dose as exported from the TPS (if available).}
\item{doseMinTPS}{Minimum dose as exported from the TPS (if available).}
\item{doseMaxTPS}{Maximum dose as exported from the TPS (if available).}
}
}
\seealso{
\code{\link[DVHmetrics]{getMetric}},
\code{\link[DVHmetrics]{convertDVHsmooth}},
\code{\link[stats]{approxfun}},
\code{\link[stats]{splinefun}},
\code{\link[stats]{smooth.spline}},
\code{\link[KernSmooth]{dpill}},
\code{\link[KernSmooth]{locpoly}}
}
\examples{
getDMEAN(dataMZ[[1]], interp="linear")
}
