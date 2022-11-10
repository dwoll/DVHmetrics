\name{run_gui}
\alias{run_gui}
\title{Open web-based GUI in browser}
\description{
Opens a web-based GUI in an external browser.
}
\usage{
run_gui(...)
}
\arguments{
  \item{...}{Arguments passed to \code{\link[shiny]{runApp}}. Supply \code{port=80} if a web browser refuses to connect to the randomly chosen port for security reasons.}
}
\details{
This function calls \code{\link[shiny]{runApp}} to run the included \code{MeshAgreement} application. Be aware that packages \code{shiny}, \code{DT}, \code{bs4Dash}, \code{sortable} need to be installed. Packages \code{PolygonSoup} and \code{Boov} have to be installed as well if union / intersection volume should be computed via \code{\link[Boov]{MeshesUnion}} / \code{\link[Boov]{MeshesIntersection}}.
}
\seealso{
\code{\link[shiny]{runApp}}
}
\examples{
\dontrun{
run_gui()
}
}
