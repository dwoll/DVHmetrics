\name{readConstraint}
\alias{readConstraint}
\title{Read constraint definitions from text file}
\description{Reads the definition of quality assurance constraints from a text file.}
\usage{readConstraint(x, ...)}
\arguments{
  \item{x}{\code{character} string giving the path to a single text file with the constraint definition. May contain globbing symbols understood by \code{\link[base]{Sys.glob}}. If missing and in interactive mode, \code{readDVH} opens a file selector widget. See Details.}
  \item{...}{Further arguments passed to \code{\link[utils]{read.table}}, e.g., \code{sep="\t"} to define the column sepator as tab.}
}
\details{
This is a wrapper for \code{\link[utils]{read.table}}.

The text file should contain three columns with the column names \code{patID}, \code{structure}, \code{constraint} in the first line. Each further line then defines one constraint and the scope it applies to in terms of patients and structures. See \code{\link[DVHmetrics]{checkConstraint}} for the definition of a constraint and for the definition of a scope. Example content:

\code{"patID" "structure" "constraint"}\cr
\code{"*"     "HEART"     "D1CC < 20Gy"}\cr
\code{"234"   "*"         "V10\% > 8CC"}
}
\value{
A \code{data.frame} with columns \code{patID}, \code{structure}, \code{constraint} that can be used in functions \code{\link[DVHmetrics]{checkConstraint}} and \code{\link[DVHmetrics]{showConstraint}}.
}
\seealso{
\code{\link[utils]{read.table}},
\code{\link[DVHmetrics]{checkConstraint}},
\code{\link[DVHmetrics]{saveConstraint}},
\code{\link[DVHmetrics]{showConstraint}}
}
\examples{
\dontrun{
readConstraint("constraint.txt")
readConstraint()
}
}
