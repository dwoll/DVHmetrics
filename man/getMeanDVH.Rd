\name{getMeanDVH}
\alias{getMeanDVH}
\alias{getMeanDVH.DVHs}
\alias{getMeanDVH.DVHLst}
\alias{getMeanDVH.DVHLstLst}
\title{Point-wise mean DVH with point-wise SDs}
\description{
Returns the point-wise mean and median DVH with the point-wise standard deviation for a given list of input DVHs. Other point-wise measures may be calculated as well.
}
\usage{
getMeanDVH(x, fun=list(mean=mean, median=median, sd=sd),
           cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
           fixed=TRUE, returnDVHObj=FALSE)

\method{getMeanDVH}{DVHs}(x, fun=list(mean=mean, median=median, sd=sd),
           cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
           fixed=TRUE, returnDVHObj=FALSE)

\method{getMeanDVH}{DVHLst}(x, fun=list(mean=mean, median=median, sd=sd),
           cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
           fixed=TRUE, returnDVHObj=FALSE)

\method{getMeanDVH}{DVHLstLst}(x, fun=list(mean=mean, median=median, sd=sd),
           cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
           fixed=TRUE, returnDVHObj=FALSE)
}
\arguments{
  \item{x}{A single DVH (object of class \code{DVHs}), multiple DVHs from one patient/structure (object of class \code{DVHLst}), or multiple DVHs from many patients/structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{fun}{Named \code{list} of functions that should be applied to yield 1 point-wise DVH measure. Functions must have exactly 1 return value.}
  \item{cumul}{\code{logical}. Get point-wise mean and SD for cumulative or differential (per unit dose) DVH?}
  \item{thin}{\code{numeric}. The number of DVH nodes (dose values) is reduced by 1/\code{thin} of the maximum number of nodes in \code{x} before interpolating and averaging.}
  \item{byPat}{\code{logical}. Relevant if multiple DVHs are given. \code{byPat=TRUE} means that for each patient, DVHs for multiple structures are averaged point wise. \code{byPat=FALSE} means that for each structure, DVHs for multiple patients averaged point wise.}
  \item{patID}{\code{character} vector. Include DVHs for these patients only when calculating mean/SD. If missing, all patients are used. Can be a regular expression with \code{fixed=FALSE}, see \code{\link[base]{regex}}.}
  \item{structure}{\code{character} vector. Include DVHs for these structures only when calculating mean/SD. If missing, all structures are used. Can be a regular expression with \code{fixed=FALSE}, see \code{\link[base]{regex}}.}
  \item{fixed}{\code{logical}. Use \code{fixed=FALSE} for regular expression matching of \code{patID} and \code{structure}.}
  \item{returnDVHObj}{\code{logical}. With \code{returnDVHObj=TRUE}, a regular DVH object is returned. In that case, only the first component of \code{fun} is used which should be \code{mean} or \code{median} (not checked).}
}
\details{
Before calculating the point-wise mean and SD, DVHs in \code{x} are first linearly interpolated with \code{\link[DVHmetrics]{convertDVH}} using the same set of nodes.
}
\value{
By default (\code{returnDVHObj=FALSE}) returns a data frame with point-wise mean DVH averaged over structures (\code{byPat=TRUE}) or over patients (\code{byPat=FALSE}) including the point-wise standard deviation or other measures as controlled by \code{fun}. With \code{returnDVHObj=TRUE}, a DVH object is returned that is equivalent to a DVH as imported from a file. In particular, functions like \code{\link[DVHmetrics]{showDVH}} or \code{\link[DVHmetrics]{getMetric}} can be used on such an object.
}
\seealso{
\code{\link[DVHmetrics]{showDVH}},
\code{\link[DVHmetrics]{convertDVH}}
}
\examples{
res1 <- getMeanDVH(dataMZ, byPat=TRUE, structure=c("HEART", "AMYOCL"))
head(res1)

# average differential DVHs
# matches patients P123 and P234
res2 <- getMeanDVH(dataMZ, fun=list(min=min, max=max),
                   cumul=FALSE, byPat=FALSE,
                   patID="23", fixed=FALSE)
head(res2)
}
