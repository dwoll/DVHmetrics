\name{getMetric}
\alias{getMetric}
\alias{getMetric.DVHs}
\alias{getMetric.DVHLst}
\alias{getMetric.DVHLstLst}
\title{Calculate dose-volume-histogram metrics}
\description{
Simultaneously calculates multiple metrics for multiple cumulative DVHs.
}
\usage{
getMetric(x, metric, patID, structure,
           sortBy=c("none", "observed", "patID", "structure", "metric"),
          splitBy=c("none", "patID", "structure", "metric"),
           interp=c("linear", "spline"), ...)

\method{getMetric}{DVHs}(x, metric, patID, structure,
           sortBy=c("none", "observed", "patID", "structure", "metric"),
          splitBy=c("none", "patID", "structure", "metric"),
           interp=c("linear", "spline"), ...)

\method{getMetric}{DVHLst}(x, metric, patID, structure,
           sortBy=c("none", "observed", "patID", "structure", "metric"),
          splitBy=c("none", "patID", "structure", "metric"),
           interp=c("linear", "spline"), ...)

\method{getMetric}{DVHLstLst}(x, metric, patID, structure,
           sortBy=c("none", "observed", "patID", "structure", "metric"),
          splitBy=c("none", "patID", "structure", "metric"),
           interp=c("linear", "spline"), ...)
}
\arguments{
  \item{x}{One cumulative DVH (object of class \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{metric}{\code{character} vector defining one or more DVH metrics. See Details for their definition. For metrics involving the relative dose, the DVH must contain the prescription dose.}
  \item{patID}{\code{character} vector. Calculate given DVH metrics for these patients only. If missing, DVH metrics are calculated for all patients. Can be a regular expression, see \code{\link[base]{regex}}.}
  \item{structure}{\code{character} vector. Calculate given DVH metrics for these structures only. If missing, DVH metrics are calculated for all structures. Can be a regular expression, see \code{\link[base]{regex}}.}
  \item{sortBy}{\code{character} vector giving the sorting criteria for the output data frame.}
  \item{splitBy}{\code{character} vector. Split results into a list of data frames where list components are defined by groups from combining these variables.}
  \item{interp}{\code{character}. Method of interpolation between DVH points.}
  \item{...}{Further arguments passed to code{\link[base]{grep}}. Use \code{fixed=TRUE} for exact matching of \code{patID} and \code{structure}.}
}
\details{
A \emph{pre-specified} DVH metric is one of the following character strings:
\itemize{
\item \code{"DMEAN"}: The volume-weighted mean dose of the structure.
\item \code{"DMIN"}: The minimum dose of the non-zero-dose voxels in the structure.
\item \code{"DMAX"}: The maximum dose of the non-zero-dose voxels in the structure.
\item \code{"DSD"}: The standard deviation of the dose in the structure.
\item \code{"DRX"}: The prescription dose.
}
A \emph{free} DVH metric is a character string which has three mandatory elements and one optional element in the following order:
\itemize{
\item 1st letter \code{"D"} or \code{"V"}: \code{"D"} If the requested value is a dose, \code{"V"} if it is a volume.
\item 2nd element \code{<number>}: If the first letter is \code{"D"}, this gives the volume for which the dose value of the cumulative DVH should be reported. If the first letter is \code{"V"}, this gives the dose for which the volume value of the cumulative DVH should be reported.
\item 3rd element \code{<measurement unit>}: The measurement unit for the 2nd element of the metric. Absolute volumes are indicated by \code{"CC"} for cubic centimeter, relative volumes by \code{"\%"}. Absolute doses are indicated by \code{"Gy"} for Gray or \code{"cGy"} for Centigray, relative doses by \code{"\%"}.
\item Optional 4th element \code{_<measurement unit>}: The measurement unit of the output value. Possible units are as for the 3rd element. If missing, dose is reported as absolute dose in the measurement unit used in the DVH. Volume is reported as relative volume in \%.
}
Examples:
\itemize{
\item \code{"D1\%"}: Minimal absolute dose for the "hottest" 1\% of the structure, i.e., the maximally irradiated 1\% of the structure was exposed to at least this absolute dose.
\item \code{"D1CC_\%"}: Minimal relative dose (\% of prescription dose) for the maximally irradiated cm^3 of the structure.
\item \code{"V500cGy"}: Relative structure volume in \% that was exposed to at least 500cGy.
\item \code{"V10\%_CC"}: Absolute structure volume in cm^3 that was exposed to at least 10\% of prescription dose.
}

If volume or dose values outside the range of possible values for a structure are requested, metrics cannot be calculated, and the result will be \code{NA} with a warning.
}
\value{
A data frame or a list with details on the calculated metrics.
\item{patID}{Patient ID}
\item{structure}{Structure}
\item{metric}{The calculated DVH metric}
\item{observed}{The observed value for the DVH metric}
}
\seealso{
\code{\link[DVHmetrics]{saveMetric}}
}
\examples{
data(dataMZ)
getMetric(dataMZ, c("D1CC", "V10\%_CC"),
          sortBy=c("metric", "structure", "observed"))

# matching patients are P123 and P234
# matching structures are AMYOCL and AMYOCR
getMetric(dataMZ, c("D1CC", "V10\%_CC"),
          patID="23",
          structure=c("AMYOC", "VALVE"),
          splitBy="patID")
}
