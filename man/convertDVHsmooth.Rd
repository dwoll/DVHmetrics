\name{convertDVHsmooth}
\alias{convertDVHsmooth}
\alias{convertDVHsmooth.matrix}
\alias{convertDVHsmooth.DVHs}
\alias{convertDVHsmooth.DVHLst}
\alias{convertDVHsmooth.DVHLstLst}
\title{Convert between differential and cumulative DVH}
\description{Convert between differential and cumulative DVH as well as between dose units, using smoothng of the differential DVH.}
\usage{
convertDVHsmooth(x,
                 toType=c("asis", "cumulative", "differential"),
                 toDoseUnit=c("asis", "GY", "CGY"),
                 interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
                 nodes=NULL, rangeD=NULL, perDose=TRUE)

\method{convertDVHsmooth}{matrix}(x,
                 toType=c("asis", "cumulative", "differential"),
                 toDoseUnit=c("asis", "GY", "CGY"),
                 interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
                 nodes=NULL, rangeD=NULL, perDose=TRUE)

\method{convertDVHsmooth}{DVHs}(x,
                 toType=c("asis", "cumulative", "differential"),
                 toDoseUnit=c("asis", "GY", "CGY"),
                 interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
                 nodes=NULL, rangeD=NULL, perDose=TRUE)

\method{convertDVHsmooth}{DVHLst}(x,
                 toType=c("asis", "cumulative", "differential"),
                 toDoseUnit=c("asis", "GY", "CGY"),
                 interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
                 nodes=NULL, rangeD=NULL, perDose=TRUE)

\method{convertDVHsmooth}{DVHLstLst}(x,
                 toType=c("asis", "cumulative", "differential"),
                 toDoseUnit=c("asis", "GY", "CGY"),
                 interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
                 nodes=NULL, rangeD=NULL, perDose=TRUE)
}
\arguments{
  \item{x}{One DVH (object of class \code{matrix} or \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{toType}{\code{character}. Convert the DVH to this type. \code{"asis"} keeps the current DVH type.}
  \item{toDoseUnit}{\code{character}. Convert the DVH to this dose unit. \code{"asis"} keeps the current dose unit.}
  \item{interp}{\code{character}. Interpolation method for the differential DVH. \code{"asis"} and \code{"linear"} for no interpolation. \code{"spline"} for spline interpolation using \code{\link[stats]{splinefun}} (\code{"fmm"} for differential, \code{"monoH.FC"} for cumulative DVHs), \code{"ksmooth"} for local polynomial regression using \code{\link[KernSmooth]{locpoly}} with kernel bandwidth chosen by the direct plug-in method using \code{\link[KernSmooth]{dpill}}, \code{"smoothSpl"} for a smoothing spline using \code{\link[stats]{smooth.spline}}, with the smoothing parameter chosen by generalized crossvalidation.}
  \item{nodes}{\code{numeric}. Minimum number of nodes to use in interpolation for method \code{"ksmooth"}. Number of available nodes is kept as is for \code{NULL} or if larger than \code{nodes}.}
  \item{rangeD}{\code{numeric}. Dose range for interpolation methods \code{"linear"}, \code{"spline"}, \code{"smoothSpl"}. If \code{NULL} it is determined individally for each DVH.}
  \item{perDose}{\code{logical}. Are the differential DVH volume values per unit dose?}
}
\value{Depending on the input, an object of class \code{matrix}, \code{DVHs}, \code{DVHLst}, or \code{DVHLstLst}.}
\seealso{
\code{\link[DVHmetrics]{convertDVH}},
\code{\link[DVHmetrics]{readDVH}},
\code{\link[DVHmetrics]{showDVH}}
}
\examples{
res <- convertDVHsmooth(dataMZ[[c(1, 1)]],
                  toType="cumulative",
                  toDoseUnit="CGY")
}
