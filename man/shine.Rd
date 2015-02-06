\name{shine}
\alias{shine}
\title{Open DVHshiny app in browser}
\description{
Opens the DVHshiny web application in an external browser.
}
\usage{
shine(...)
}
\arguments{
  \item{...}{Arguments passed to \code{\link[shiny]{runApp}}. Supply \code{port=80} if a web browser refuses to connect to the randomly chosen port for security reasons.}
}
\details{
This function calls \code{\link[shiny]{runApp}} to run the included DVHshiny application. See \code{vignette("DVHshiny")} for documentation.
}
\seealso{
\code{\link[shiny]{runApp}}
}
\examples{
\dontrun{
shine()
}
}
