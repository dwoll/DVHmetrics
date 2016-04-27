\name{showMeanDVH}
\alias{showMeanDVH}
\title{Show average dose volume histograms}
\description{
Displays average dose volume histograms.
}
\usage{
showDVH(x, byPat=TRUE, patID=NULL, structure=NULL,
        rel=TRUE, guessX=TRUE, thresh=1, show=TRUE, fixed=TRUE,
        showSD=TRUE, color=TRUE)
}
\arguments{
  \item{x}{A data frame as returned by \code{\link[DVHmetrics]{getMeanDVH}} or a list of such data frames. See Details.}
  \item{byPat}{\code{logical}. Relevant if multiple DVHs are given. If \code{x} has class \code{DVHLstLst}: \code{byPat=TRUE} means that one diagram shows DVHs from one patient with multiple structures. \code{byPat=FALSE} means that one diagram shows DVHs for one structure from multiple patients.}
  \item{patID}{\code{character} vector. Show diagram for these patients only. If missing, all patients are shown. Can be a regular expression with \code{fixed=FALSE}, see \code{\link[base]{regex}}.}
  \item{structure}{\code{character} vector. Show diagram for these structures only. If missing, all structures are shown. Can be a regular expression with \code{fixed=FALSE}, see \code{\link[base]{regex}}.}
  \item{rel}{\code{logical}. Show relative volume?}
  \item{guessX}{\code{logical}. Try to clip the x-axis for better visibility of main DVH range?}
  \item{thresh}{\code{numeric} value. Relative volume threshold used with \code{guessX=TRUE}. Clip x-axis (+10\%) such that the "highest" DVH is cut off at this relative volume.}
  \item{show}{\code{logical}. If \code{TRUE}, diagrams are shown, if \code{FALSE} diagrams are not shown - only \code{\link[ggplot2]{ggplot}} diagram objects are silently returned.}
  \item{fixed}{\code{logical}. Use \code{fixed=FALSE} for regular expression matching of \code{patID} and \code{structure}.}
  \item{showSD}{\code{logical}. If \code{TRUE}, diagram shows shaded areas for point-wise 1-standard deviation and 2-standard deviations around this mean. See details.}
  \item{color}{\code{logical}. If \code{TRUE}, diagram uses color to distinguish groups. If \code{FALSE}, colors are greyscale, and line types are used to distinguish groups.}
}
\details{
  TODO
}
\value{
Silently returns a \code{\link[ggplot2]{ggplot}} diagram object, or - when multiple diagrams are constructed - a list of \code{\link[ggplot2]{ggplot}} diagram objects.
}
\seealso{
\code{\link[ggplot2]{ggplot}},
\code{\link[DVHmetrics]{showDVH}},
\code{\link[DVHmetrics]{getMeanDVH}}
}
\examples{
# TODO
}
