\name{saveMetric}
\alias{saveMetric}
\alias{saveMetric.data.frame}
\alias{saveMetric.list}
\title{Save DVH metrics to file}
\description{Saves results from \code{\link[DVHmetrics]{getMetric}} to a text file.}
\usage{
saveMetric(x, file = "", ...)

\method{saveMetric}{data.frame}(x, file = "", ...)

\method{saveMetric}{list}(x, file = "", ...)
}
\arguments{
  \item{x}{\code{data.frame} or \code{list} - the result from \code{\link[DVHmetrics]{getMetric}}.}
  \item{file}{\code{character} string. Path to file.}
  \item{...}{Further arguments passed to \code{\link[utils]{write.table}} - e.g., \code{dec="."} to define the decimal separator as point or \code{sep="\t"} to define the column sepator as tab.}
}
\details{This is a wrapper for \code{\link[utils]{write.table}}.}
\value{If \code{x} is a \code{list}, one text file is written for each list component. If \code{x} is a \code{data.frame}, one file is written.}
\seealso{
\code{\link[utils]{write.table}},
\code{\link[DVHmetrics]{getMetric}}
}
\examples{
data(dataMZ)
res <- getMetric(dataMZ, c("D1CC", "V10\%_CC"),
                 sortBy=c("metric", "structure"),
                 splitBy="patID")
\dontrun{
# not run
saveMetric(res, file="metricsResults.txt", sep="\t")
}
}
