\name{getMeanDVH}
\alias{getMeanDVH}
\alias{getMeanDVH.DVHs}
\alias{getMeanDVH.DVHLst}
\alias{getMeanDVH.DVHLstLst}
\title{Point-wise mean DVH with point-wise SDs}
\description{
Returns the point-wise mean DVH as well as point-wise 1-standard deviation and 2-standard deviations around this mean for a given list of input DVHs.
}
\usage{
getMeanDVH(x, cumul=TRUE, purge=3, byPat=TRUE, patID=NULL, structure=NULL,
           interp=c("linear", "spline", "smoothSpl"), fixed=TRUE)

\method{getMeanDVH}{DVHs}(x, cumul=TRUE, purge=3, byPat=TRUE, patID=NULL, structure=NULL,
           interp=c("linear", "spline", "smoothSpl"), fixed=TRUE)

\method{getMeanDVH}{DVHLst}(x, cumul=TRUE, purge=3, byPat=TRUE, patID=NULL, structure=NULL,
           interp=c("linear", "spline", "smoothSpl"), fixed=TRUE)

\method{getMeanDVH}{DVHLstLst}(x, cumul=TRUE, purge=3, byPat=TRUE, patID=NULL, structure=NULL,
           interp=c("linear", "spline", "smoothSpl"), fixed=TRUE)
}
\arguments{
  \item{x}{A single DVH (object of class \code{DVHs}), multiple DVHs from one patient/structure (object of class \code{DVHLst}), or multiple DVHs from many patients/structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{cumul}{\code{logical}. Get point-wise mean and SD for cumulative or differential (per unit dose) DVH?}
  \item{purge}{\code{numeric}. The number of DVH nodes (dose values) is reduced by 1/\code{purge} of the maximum number of nodes in \code{x} before interpolating and averaging.}
  \item{byPat}{\code{logical}. Relevant if multiple DVHs are given. \code{byPat=TRUE} means that for each patient, DVHs for multiple structures are averaged point wise. \code{byPat=FALSE} means that for each structure, DVHs for multiple patients averaged point wise.}
  \item{patID}{\code{character} vector. Include DVHs for these patients only when calculating mean/SD. If missing, all patients are used. Can be a regular expression with \code{fixed=FALSE}, see \code{\link[base]{regex}}.}
  \item{structure}{\code{character} vector. Include DVHs for these structures only when calculating mean/SD. If missing, all structures are used. Can be a regular expression with \code{fixed=FALSE}, see \code{\link[base]{regex}}.}
  \item{interp}{\code{character}. Interpolation method for the differential DVH. \code{"asis"} and \code{"linear"} for no interpolation. \code{"spline"} for spline interpolation using \code{\link[stats]{splinefun}} (\code{"fmm"} for differential, \code{"monoH.FC"} for cumulative DVHs), \code{"smoothSpl"} for a smoothing spline using \code{\link[stats]{smooth.spline}}, with the smoothing parameter chosen by generalized crossvalidation.}
  \item{fixed}{\code{logical}. Use \code{fixed=FALSE} for regular expression matching of \code{patID} and \code{structure}.}
}
\details{
Before calculating the point-wise mean andSD, DVHs in \code{x} are first linearly interpolated with \code{\link[DVHmetrics]{convertDVH}} using the same set of nodes.
}
\value{
Returns a data frame with point-wise mean DVH averaged over structures (\code{byPat=TRUE}) or over patients (\code{byPat=FALSE}).
}
\seealso{
\code{\link[DVHmetrics]{showDVH}},
\code{\link[DVHmetrics]{convertDVH}}
}
\examples{
res1 <- getMeanDVH(dataMZ, byPat=TRUE, structure=c("HEART", "AMYOCL"))
head(res1)

# average differential DVHs, use kernel smoothing
# matches patients P123 and P234
res2 <- getMeanDVH(dataMZ, cumul=FALSE, byPat=FALSE,
                   patID="23", interp="smoothSpl", fixed=FALSE)
head(res2)
}