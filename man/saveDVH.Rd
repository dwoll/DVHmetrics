\name{saveDVH}
\alias{saveDVH}
\title{Save DVH diagram to file}
\description{Saves one or multiple DVH diagrams to file.}
\usage{
saveDVH(x, file="", ...)
}
\arguments{
  \item{x}{A single \code{\link[ggplot2]{ggplot}} object or a list of multiple \code{\link[ggplot2]{ggplot}} objects as returned by \code{\link[DVHmetrics]{showDVH}} or \code{\link[DVHmetrics]{showConstraint}}.}
  \item{file}{\code{character}. Path to file. The file-ending determines what kind of file is written, e.g., \code{"filename.pdf"} will write a pdf document, \code{"filename.jpg"} a JPEG image.}
  \item{...}{Further arguments passed to \code{\link[ggplot2]{ggsave}}, e.g., \code{width} and \code{height} to determine the figure size.}
}
\details{This is a wrapper for \code{\link[ggplot2]{ggsave}}.}
\value{If \code{x} is a \code{list} of \code{\link[ggplot2]{ggplot}} objects, one file is written for each list component. If \code{x} is a single \code{\link[ggplot2]{ggplot}} object, one file is written.}
\seealso{
\code{\link[ggplot2]{ggsave}},
\code{\link[DVHmetrics]{showDVH}},
\code{\link[DVHmetrics]{showConstraint}}
}
\examples{
res <- showDVH(dataMZ, byPat=TRUE, structure=c("HEART", "AMYOCL"))
\dontrun{
saveDVH(res, "out.pdf")
}
}
