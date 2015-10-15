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
           interp=c("linear", "spline", "ksmooth"), fixed=TRUE, ...)

\method{getMetric}{DVHs}(x, metric, patID, structure,
           sortBy=c("none", "observed", "patID", "structure", "metric"),
          splitBy=c("none", "patID", "structure", "metric"),
           interp=c("linear", "spline", "ksmooth"), fixed=TRUE, ...)

\method{getMetric}{DVHLst}(x, metric, patID, structure,
           sortBy=c("none", "observed", "patID", "structure", "metric"),
          splitBy=c("none", "patID", "structure", "metric"),
           interp=c("linear", "spline", "ksmooth"), fixed=TRUE, ...)

\method{getMetric}{DVHLstLst}(x, metric, patID, structure,
           sortBy=c("none", "observed", "patID", "structure", "metric"),
          splitBy=c("none", "patID", "structure", "metric"),
           interp=c("linear", "spline", "ksmooth"), fixed=TRUE, ...)
}
\arguments{
  \item{x}{One cumulative DVH (object of class \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{metric}{\code{character} vector defining one or more DVH metrics. See Details for their definition. For metrics involving the relative dose, the DVH must contain the prescription dose.}
  \item{patID}{\code{character} vector. Calculate given DVH metrics for these patients only. If missing, DVH metrics are calculated for all patients. Can be a regular expression if additional argument \code{fixed=FALSE} is supplied as well, see \code{\link[base]{regex}}.}
  \item{structure}{\code{character} vector. Calculate given DVH metrics for these structures only. If missing, DVH metrics are calculated for all structures. Can be a regular expression if additional argument \code{fixed=FALSE} is supplied as well, see \code{\link[base]{regex}}.}
  \item{sortBy}{\code{character} vector giving the sorting criteria for the output data frame.}
  \item{splitBy}{\code{character} vector. Split results into a list of data frames where list components are defined by groups from combining these variables.}
  \item{interp}{\code{character}. Method of interpolation between DVH points: Linear interpolation using \code{\link[stats]{approx}}, monotone Hermite spline interpolation using \code{\link[stats]{splinefun}}, or local polynomial regression using \code{\link[KernSmooth]{locpoly}} with kernel bandwidth chosen by the direct plug-in method using \code{\link[KernSmooth]{dpill}}.}
  \item{fixed}{\code{logical}. Use \code{fixed=FALSE} for regular expression matching of \code{patID} and \code{structure}.}
  \item{...}{Further arguments passed to \code{\link[DVHmetrics]{getEUD}} (for \code{metric="DEUD"}), \code{\link[DVHmetrics]{getTCP}} (for \code{metric="DTCP"}), or \code{\link[DVHmetrics]{getNTCP}} (for \code{metric="DNTCP"}).}
}
\details{
A \emph{pre-specified} DVH metric is one of the following character strings:
\itemize{
\item \code{"DMEAN"}: The volume-weighted mean dose of the structure.
\item \code{"DMEDIAN"}: Median dose, equal to D50\%
\item \code{"DMIN"}: The minimum dose of the non-zero-dose voxels in the structure.
\item \code{"DMAX"}: The maximum dose of the non-zero-dose voxels in the structure.
\item \code{"DSD"}: The standard deviation of the dose in the structure.
\item \code{"DRX"}: The prescription dose.
\item \code{"DEUD"}: The generalized equivalent uniform dose (gEUD). See \code{\link[DVHmetrics]{getEUD}} for mandatory and optional parameters.
\item \code{"DNTCP"}: The normal tissue complication probability (NTCP). See \code{\link[DVHmetrics]{getNTCP}} for mandatory and optional parameters.
\item \code{"DTCP"}: The tumor control probability (TCP). See \code{\link[DVHmetrics]{getNTCP}} for mandatory and optional parameters.
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

DMEAN, DMEDIAN, DMIN, DMAX, DSD are taken from the exported DVH if present. Otherwise, the differential DVH is generated and used for calculating these metrics.
}
\value{
A data frame or a list with details on the calculated metrics.
\item{patID}{Patient ID}
\item{structure}{Structure}
\item{metric}{The calculated DVH metric}
\item{observed}{The observed value for the DVH metric}
}
\references{
Rancati et al. (2004). Fitting late rectal bleeding data using different NTCP models: results from an Italian multi-centric study (AIROPROS0101). Radiotherapy Oncology, 73, 21-32.

Wu et al. (2002). Optimization of intensity modulated radiotherapy plans based on the equivalent uniform dose. International Journal of Radiation Oncology Biology Physics, 52, 224-235.
}
\seealso{
\code{\link[DVHmetrics]{saveMetric}},
\code{\link[DVHmetrics]{getEUD}},
\code{\link[DVHmetrics]{getNTCP}},
\code{\link[DVHmetrics]{getTCP}},
\code{\link[DVHmetrics]{getEQD2}},
\code{\link[stats]{approxfun}},
\code{\link[stats]{splinefun}},
\code{\link[KernSmooth]{dpill}},
\code{\link[KernSmooth]{locpoly}}
}
\examples{
getMetric(dataMZ, c("D1CC", "V10\%_CC"),
          sortBy=c("metric", "structure", "observed"))

# matching patients are P123 and P234
# matching structures are AMYOCL and AMYOCR
getMetric(dataMZ, c("D1CC", "V10\%_CC"),
          patID="23",
          structure=c("AMYOC", "VALVE"),
          splitBy="patID",
          fixed=FALSE)

# gEUD with a=2
getMetric(dataMZ[[c(1, 1)]], "DEUD", EUDa=2)

# gEUD based on EQD2 with a=2, 20 fractions
getMetric(dataMZ[[c(1, 1)]], "DEUD", EUDa=2, EUDfd=1.8)

# NTCP Lyman probit model with TD50=20, m=4, n=0.5
getMetric(dataMZ[[c(1, 1)]], "DNTCP",
          NTCPtd50=20, NTCPm=4, NTCPn=0.5, NTCPtype="probit")
}
