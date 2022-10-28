\name{run_gui}
\alias{run_gui}
\title{Open web-based GUI in browser}
\description{
Opens the web-based GUI in an external browser.
}
\usage{
run_gui(...)
}
\arguments{
  \item{...}{Arguments passed to \code{\link[shiny]{runApp}}. Supply \code{port=80} if a web browser refuses to connect to the randomly chosen port for security reasons.}
}
\details{
This function calls \code{\link[shiny]{runApp}} to run the included DVHshiny application.
}
\seealso{
\code{\link[shiny]{runApp}}
}
\examples{
\dontrun{
run_gui()
}
}
