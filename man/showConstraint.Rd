\name{showConstraint}
\alias{showConstraint}
\alias{showConstraint.DVHs}
\alias{showConstraint.DVHLst}
\alias{showConstraint.DVHLstLst}
\title{Display constraints for cumulative dose-volume histograms}
\description{
Displays quality assurance constraints for cumulative dose-volume histograms: Either one diagram per patient - including multiple structures. Or one diagram per structure - including multiple patients.
}
\usage{
showConstraint(x, constr, byPat=TRUE, rel=TRUE, guessX=TRUE, guessY=TRUE,
               thresh=1, show=TRUE, visible=FALSE)

\method{showConstraint}{DVHs}(x, constr, byPat=TRUE, rel=TRUE, guessX=TRUE, guessY=TRUE,
               thresh=1, show=TRUE, visible=FALSE)

\method{showConstraint}{DVHLst}(x, constr, byPat=TRUE, rel=TRUE, guessX=TRUE, guessY=TRUE,
               thresh=1, show=TRUE, visible=FALSE)

\method{showConstraint}{DVHLstLst}(x, constr, byPat=TRUE, rel=TRUE, guessX=TRUE, guessY=TRUE,
               thresh=1, show=TRUE, visible=FALSE)
}
\arguments{
  \item{x}{A single DVH (object of class \code{DVHs}), multiple DVHs from one patient/structure (object of class \code{DVHLst}), or multiple DVHs from many patients/structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}. See Details.}
  \item{constr}{One or more constraints - given as a \code{character} vector or as a \code{data.frame}. See \code{\link[DVHmetrics]{checkConstraint}} for their definition.}
  \item{byPat}{\code{logical}. Relevant if multiple DVHs are given. If \code{x} has class \code{DVHLstLst}: \code{byPat=TRUE} means that one diagram shows DVHs from one patient with multiple structures. \code{byPat=FALSE} means that one diagram shows DVHs for one structure from multiple patients.}
  \item{rel}{\code{logical}. Show relative volume?}
  \item{guessX}{\code{logical}. Try to guess the best x-axis limits for better visibility of main DVH range? If \code{FALSE}, x-axis runs from 0 to maximum dose. If \code{TRUE}, x-axis runs from 0 to dose value where volume approaches 0. If a single number is given, it is interpreted as the maximum value. If a vector of two numbers is given, it is interpreted as the range of the axis.}
  \item{guessY}{\code{logical}. Try to guess the best y-axis limits? If a single number is given, it is interpreted as the maximum value. If a vector of two numbers is given, it is interpreted as the range of the axis.}
  \item{thresh}{\code{numeric} value. Relative volume threshold used with \code{guessX=TRUE}. Clip x-axis (+10\%) such that the "highest" DVH is cut off at this relative volume.}
  \item{show}{\code{logical}. If \code{TRUE}, diagrams are shown, if \code{FALSE} diagrams are not shown - only \code{\link[ggplot2]{ggplot}} diagram objects are silently returned.}
  \item{visible}{\code{logical}. Return \code{\link[ggplot2]{ggplot}} diagram object visibly or invisibly. \code{show=FALSE} with \code{visible=TRUE} is useful for zooming in shiny apps.}
}
\details{
Constraints are shown as points in the cumulative DVH with an additional arrow indicating where the cumulative DVH curve should lie relative to the constraint. On each DVH curve, the point with the minimal Euclidean distance to the constraint is indicated. Note that, visually, this point only has the minimal apparent distance if the aspect ratio of the diagram is 1.

If multiple diagrams are produced, they are shown in the same graphics device. If interactive inspection is required, make sure you use an R development environment that saves previous diagrams and allows navigating between them - e.g., \href{https://www.rstudio.com/products/rstudio/download/}{RStudio} or \href{https://www.getarchitect.io/}{OpenAnalytics Architect}.
}
\value{
Silently returns a \code{\link[ggplot2]{ggplot}} diagram object, or - when multiple diagrams are constructed - a list of \code{\link[ggplot2]{ggplot}} diagram objects.
}
\seealso{
\code{\link[DVHmetrics]{checkConstraint}},
\code{\link[DVHmetrics]{saveDVH}}
}
\examples{
data(dataMZ)

# define constraints
constr <- data.frame(
    patID=c("P123", "P234"),
    structure=c("HEART", "*"),
    constraint=c("D1CC < 20Gy", "V10\% > 8CC"),
    stringsAsFactors=FALSE)       # this is important
showConstraint(dataMZ, constr=constr, byPat=FALSE)
}
