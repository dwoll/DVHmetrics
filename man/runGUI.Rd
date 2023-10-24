\name{runGUI}
\alias{runGUI}
\title{Open web-based GUI in browser}
\description{
Opens the web-based GUI in an external browser.
}
\usage{
runGUI(...)
}
\arguments{
  \item{...}{Arguments passed to \code{\link[shiny]{runApp}}. Supply \code{port=80} if a web browser refuses to connect to the randomly chosen port for security reasons.}
}
\details{
This function is a wrapper for \code{\link[shiny]{runApp}} which runs the included DVHshiny application. See \code{vignette("DVHshiny")} for documentation.
}
\seealso{
\code{\link[shiny]{runApp}}
}
\examples{
\dontrun{
runGUI()
}
}
