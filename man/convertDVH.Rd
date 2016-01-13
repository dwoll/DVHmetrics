\name{convertDVH}
\alias{convertDVH}
\alias{convertDVH.matrix}
\alias{convertDVH.DVHs}
\alias{convertDVH.DVHLst}
\alias{convertDVH.DVHLstLst}
\title{Convert between differential and cumulative DVH}
\description{Convert between differential and cumulative DVH as well as between dose units.}
\usage{
convertDVH(x, toType=c("asis", "cumulative", "differential"),
           toDoseUnit=c("asis", "GY", "CGY"),
           interp=c("asis", "linear"),
           nodes=NULL, rangeD=NULL, perDose=TRUE)

\method{convertDVH}{matrix}(x, toType=c("asis", "cumulative", "differential"),
           toDoseUnit=c("asis", "GY", "CGY"),
           interp=c("asis", "linear"),
           nodes=NULL, rangeD=NULL, perDose=TRUE)

\method{convertDVH}{DVHs}(x, toType=c("asis", "cumulative", "differential"),
           toDoseUnit=c("asis", "GY", "CGY"),
           interp=c("asis", "linear"),
           nodes=NULL, rangeD=NULL, perDose=TRUE)

\method{convertDVH}{DVHLst}(x, toType=c("asis", "cumulative", "differential"),
           toDoseUnit=c("asis", "GY", "CGY"),
           interp=c("asis", "linear"),
           nodes=NULL, rangeD=NULL, perDose=TRUE)

\method{convertDVH}{DVHLstLst}(x, toType=c("asis", "cumulative", "differential"),
           toDoseUnit=c("asis", "GY", "CGY"),
           interp=c("asis", "linear"),
           nodes=NULL, rangeD=NULL, perDose=TRUE)
}
\arguments{
  \item{x}{One DVH (object of class \code{matrix} or \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{toType}{\code{character}. Convert the DVH to this type. \code{"asis"} keeps the current DVH type.}
  \item{toDoseUnit}{\code{character}. Convert the DVH to this dose unit. \code{"asis"} keeps the current dose unit.}
  \item{interp}{\code{character}. Interpolation method for the cumulative DVH. \code{"asis"} for no interpolation and \code{"linear"} for linear interpolation.}
  \item{nodes}{\code{numeric}. Minimum number of nodes to use in linear interpolation. Number of available nodes is kept as is for \code{NULL} or if larger than \code{nodes}.}
  \item{rangeD}{\code{numeric}. Dose range for linear interpolation method. If \code{NULL} it is determined individally for each DVH.}
  \item{perDose}{\code{logical}. Are the differential DVH volume values per unit dose?}
}
\value{Depending on the input, an object of class \code{matrix}, \code{DVHs}, \code{DVHLst}, or \code{DVHLstLst}.}
\seealso{
\code{\link[DVHmetrics]{convertDVHsmooth}},
\code{\link[DVHmetrics]{readDVH}},
\code{\link[DVHmetrics]{showDVH}}
}
\examples{
res <- convertDVH(dataMZ[[c(1, 1)]],
                  toType="cumulative",
                  toDoseUnit="CGY")
}
