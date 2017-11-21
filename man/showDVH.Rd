\name{showDVH}
\alias{showDVH}
\alias{showDVH.DVHs}
\alias{showDVH.DVHLst}
\alias{showDVH.DVHLstLst}
\title{Display dose volume histograms}
\description{
Displays dose volume histograms: Either one diagram per patient - including multiple structures. Or one diagram per structure - including multiple patients.
}
\usage{
showDVH(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
        rel=TRUE, guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, fixed=TRUE)

\method{showDVH}{DVHs}(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
        rel=TRUE, guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, fixed=TRUE)

\method{showDVH}{DVHLst}(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
        rel=TRUE, guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, fixed=TRUE)

\method{showDVH}{DVHLstLst}(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
        rel=TRUE, guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, fixed=TRUE)
}
\arguments{
  \item{x}{A single DVH (object of class \code{DVHs}), multiple DVHs from one patient/structure (object of class \code{DVHLst}), or multiple DVHs from many patients/structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}. See Details.}
  \item{cumul}{\code{logical}. Show cumulative or differential (per unit dose) DVH?}
  \item{byPat}{\code{logical}. Relevant if multiple DVHs are given. If \code{x} has class \code{DVHLstLst}: \code{byPat=TRUE} means that one diagram shows DVHs from one patient with multiple structures. \code{byPat=FALSE} means that one diagram shows DVHs for one structure from multiple patients.}
  \item{patID}{\code{character} vector. Show diagram for these patients only. If missing, all patients are shown. Can be a regular expression with \code{fixed=FALSE}, see \code{\link[base]{regex}}.}
  \item{structure}{\code{character} vector. Show diagram for these structures only. If missing, all structures are shown. Can be a regular expression with \code{fixed=FALSE}, see \code{\link[base]{regex}}.}
  \item{rel}{\code{logical}. Show relative volume?}
  \item{guessX}{\code{logical}. Try to clip the x-axis for better visibility of main DVH range?}
  \item{thresh}{\code{numeric} value. Relative volume threshold used with \code{guessX=TRUE}. Clip x-axis (+10\%) such that the "highest" DVH is cut off at this relative volume.}
  \item{addMSD}{\code{logical}. If \code{TRUE}, diagram shows the point-wise mean DVH as well as shaded areas for point-wise 1-standard deviation and 2-standard deviations around this mean. See details.}
  \item{show}{\code{logical}. If \code{TRUE}, diagrams are shown, if \code{FALSE} diagrams are not shown - only \code{\link[ggplot2]{ggplot}} diagram objects are silently returned.}
  \item{fixed}{\code{logical}. Use \code{fixed=FALSE} for regular expression matching of \code{patID} and \code{structure}.}
}
\details{
If multiple diagrams are produced, they are shown in the same graphics device. If interactive inspection is required, make sure you use an R development environment that saves previous diagrams and allows navigating between them - e.g., \href{https://www.rstudio.com/products/rstudio/download/}{RStudio} or \href{https://www.getarchitect.io/}{OpenAnalytics Architect}.

For \code{addMSD=TRUE}, the number of DVH nodes (dose values) is reduced by 1/3 of the maximum number of nodes in \code{x}. Before calculating the point-wise mean and SD, DVHs in \code{x} are first linearly interpolated using the same set of nodes.
}
\value{
Silently returns a \code{\link[ggplot2]{ggplot}} diagram object, or - when multiple diagrams are constructed - a list of \code{\link[ggplot2]{ggplot}} diagram objects.
}
\seealso{
\code{\link[ggplot2]{ggplot}},
\code{\link[DVHmetrics]{readDVH}},
\code{\link[DVHmetrics]{saveDVH}},
\code{\link[DVHmetrics]{getMeanDVH}}
}
\examples{
showDVH(dataMZ, byPat=TRUE, structure=c("HEART", "AMYOCL"))

# matches patients P123 and P234
showDVH(dataMZ, byPat=FALSE, patID="23", fixed=FALSE)
}
