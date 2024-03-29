\name{checkConstraint}
\alias{checkConstraint}
\alias{checkConstraint.DVHs}
\alias{checkConstraint.DVHLst}
\alias{checkConstraint.DVHLstLst}
\title{Check constraints on dose-volume histograms (DVH)}
\description{
Simultaneously checks one or more quality assurance constraints on one or more DVHs. Reports compliance with each constraint as well as observed difference between linearly interpolated DVH and the given constraints in terms of (relative) dose, (relative) volume, and (relative) minimal Euclidean distance.
}
\usage{
checkConstraint(x, constr, byPat=TRUE, semSign=FALSE,
                sortBy=c("none", "observed", "compliance", "structure",
                         "constraint", "patID", "deltaV", "deltaD",
                         "dstMin", "dstMinRel"),
                interp=c("linear", "spline", "smooth"), ...)

\method{checkConstraint}{DVHs}(x, constr, byPat=TRUE, semSign=FALSE,
                sortBy=c("none", "observed", "compliance", "structure",
                         "constraint", "patID", "deltaV", "deltaD",
                         "dstMin", "dstMinRel"),
                interp=c("linear", "spline", "smooth"), ...)

\method{checkConstraint}{DVHLst}(x, constr, byPat=TRUE, semSign=FALSE,
                sortBy=c("none", "observed", "compliance", "structure",
                         "constraint", "patID", "deltaV", "deltaD",
                         "dstMin", "dstMinRel"),
                interp=c("linear", "spline", "smooth"), ...)

\method{checkConstraint}{DVHLstLst}(x, constr, byPat=TRUE, semSign=FALSE,
                sortBy=c("none", "observed", "compliance", "structure",
                         "constraint", "patID", "deltaV", "deltaD",
                         "dstMin", "dstMinRel"),
                interp=c("linear", "spline", "smooth"), ...)
}
\arguments{
  \item{x}{A single DVH (object of class \code{DVHs}), multiple DVHs from one patient/structure (object of class \code{DVHLst}), or multiple DVHs from many patients/structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{constr}{One or more constraints - given as a \code{character} vector or as a \code{data.frame}. See Details.}
  \item{byPat}{\code{logical}. Relevant if multiple DVHs are given. If \code{x} has class \code{DVHLst}: \code{byPat=TRUE} means that the DVHs are for one patient with multiple structures. \code{byPat=FALSE} means that the DVHs are for one structure from multiple patients. If \code{x} has class \code{DVHLstLst}: \code{byPat=TRUE} means that the DVHs are for multiple patients (list components of \code{x}) with multiple structures. \code{byPat=FALSE} means that the DVHs are for multiple structures (list components of \code{x}) from multiple patients.}
  \item{semSign}{\code{logical}. Meaning of the sign of the observed dose/volume differences between DVHs and constraints. \code{semSign=TRUE} means that negative differences indicate constraint compliance, positive differences indicate constraint violations. With \code{semSign=FALSE}, the algebraic differences are returned as is.}
  \item{sortBy}{\code{character} vector. Sorting criteria for the output data frame.}
  \item{interp}{\code{character}. Method of interpolation between DVH points: Linear interpolation using \code{\link[stats:approxfun]{approx}}, monotone Hermite spline interpolation using \code{\link[stats:splinefun]{spline}}, or local polynomial regression using \code{\link[KernSmooth]{locpoly}} with kernel bandwidth chosen by the direct plug-in method using \code{\link[KernSmooth]{dpill}}.}
  \item{...}{Additional parameters passed to \code{\link[DVHmetrics]{getMetric}}. Use for constraints on EUD (see \code{\link[DVHmetrics]{getEUD}} for parameter names), TCP (see \code{\link[DVHmetrics]{getTCP}}), and NTCP (see \code{\link[DVHmetrics]{getNTCP}}).}
}
\details{
A DVH constraint is a character string that consists of three parts: The DVH metric, the comparison operator (\code{<}, \code{>}, \code{<=}, \code{>=}), and the reference value together with the measurement unit. See \code{\link[DVHmetrics]{getMetric}} for defining a DVH metric, as well as for possible measurement units for dose and volume. For constraints involving the relative dose, the DVH must contain the prescription dose.

Some example constraints are \code{"V10Gy > 80\%"} (more than 80\% of the structure should have received 10Gy), \code{"V20\% < 10CC"} (less than 10cm^3 of the structure should have received 20\% of the prescription dose), or \code{"D10CC > 500cGy"} (The "hottest" 10cm^3 of the structure should have received more than 500cGy).

For constraints on \code{DEUD}, \code{DNTCP} and \code{DTCP} (see \code{\link[DVHmetrics]{getMetric}}), the reference measurement unit currently must be \code{Gy}, \code{cGy}, even though NTCP and TCP are probabilities. Example: \code{"DNTCP < 0.5Gy"}. The fix for this limitation is an open TODO.

A DVH constraint can apply to a specific patient or to all patients, and to a specific structure or to all structures.
\itemize{
\item If constraints apply to all patients/structures in \code{x}, \code{constr} can be a \code{character} vector with elements such as the examples above.
\item If constraints apply only to some patients/structures, \code{constr} must be a data frame with variables \code{constraint}, \code{patID} and \code{structure}. Each row then defines one constraint and its scope: \code{constraint} must be a character string with one constraint definition as in the examples above. \code{patID} must be either a character string with a valid patient ID or \code{"*"} if the the constraint applies to all patients. \code{structure} must be either a character string with a valid structure or \code{"*"} if the the constraint applies to all structures. If variable \code{patID} is missing from the data frame, the constraints apply to all available patients. If variable \code{structure} is missing from the data frame, the constraints apply to all available structures. See \code{\link[DVHmetrics]{readConstraint}} for reading appropriate constraint \code{data.frames} from external text files.
}

For calculating the minimal Euclidean distance between the constraint point and the DVH, the constraint point is orthogonally projected onto each DVH segment between (interpolated) DVH nodes. The relative Euclidean distance is the minimum of these distances divided by the distance of the constraint point to the closer one of both axes (dose and volume).

If volume or dose values outside the range of possible values for a structure are requested, metrics cannot be calculated, and the result will be \code{NA} with a warning.
}
\value{
A data frame with details on constraint compliance / violation.
\item{patID}{Patient ID}
\item{structure}{Structure}
\item{constraint}{The checked constraint}
\item{observed}{The observed value for the metric given in the constraint}
\item{compliance}{Does the DVH satisfy the constraint?}
\item{deltaV}{Volume difference between constraint and observed DVH (for the constraint dose) in measurement unit specified by constraint}
\item{deltaVpc}{Percent volume difference between constraint and observed DVH (for the constraint dose) relative to constraint volume}
\item{deltaD}{Dose difference between constraint and observed DVH (for the constraint volume) in measurement unit specified by constraint}
\item{deltaDpc}{Percent dose difference between constraint and observed DVH (for the constraint volume) relative to constraint dose}
\item{dstMin}{Minimal Euclidean distance between constraint and the cumulative DVH, using linear interpolation}
\item{ptMinD}{Dose coordinate of closest point on cumulative DVH to constraint}
\item{ptMinV}{Volume coordinate of closest point on cumulative DVH to constraint}
}
\seealso{
\code{\link[DVHmetrics]{getMetric}},
\code{\link[DVHmetrics]{getEUD}},
\code{\link[DVHmetrics]{getNTCP}},
\code{\link[DVHmetrics]{getTCP}},
\code{\link[DVHmetrics]{readConstraint}},
\code{\link[DVHmetrics]{saveConstraint}},
\code{\link[DVHmetrics]{showConstraint}}
}
\examples{
res <- checkConstraint(dataMZ, c("D10CC < 10Gy", "V20Gy < 20\%"))
head(res)

# define constraints
constr <- data.frame(
    patID=c("P123", "P234"),
    structure=c("HEART", "*"),
    constraint=c("D1CC < 20Gy", "V10\% > 8CC"),
    stringsAsFactors=FALSE)       # this is important
checkConstraint(dataMZ, constr=constr)
}
