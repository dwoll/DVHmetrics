\name{saveConstraint}
\alias{saveConstraint}
\title{Save constraint result to file}
\description{Saves results from \code{\link[DVHmetrics]{checkConstraint}} to a text file.}
\usage{saveConstraint(x, ...)}
\arguments{
  \item{x}{\code{data.frame} - the result from \code{\link[DVHmetrics]{checkConstraint}}.}
  \item{...}{Further arguments passed to \code{\link[utils]{write.table}} - e.g., \code{file="<filename>"} for the output filename, \code{dec="."} to define the decimal separator as point or \code{sep="\t"} to define the column sepator as tab.}
}
\details{This is a wrapper for \code{\link[utils]{write.table}}.}
\seealso{
\code{\link[utils]{write.table}},
\code{\link[DVHmetrics]{checkConstraint}}
}
\examples{
data(dataMZ)
res <- checkConstraint(dataMZ, c("D10CC < 10Gy", "V20Gy < 20\%"))
\dontrun{
# not run
saveConstraint(res, file="constrResults.txt", sep="\t")
}
}
