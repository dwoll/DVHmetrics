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
  \item{...}{Arguments passed to \code{\link[shiny]{runApp}}.}
}
\details{
This function calls \code{\link[shiny]{runApp}} to run the included DVHshiny application.
}
\seealso{
\code{\link[shiny]{runApp}}
}
\examples{
\dontrun{
shine()
}
}
