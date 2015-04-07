\name{dataConstr}
\alias{dataConstr}
\docType{data}
\title{Constraint data frame}
\description{
Data frame with quality assurance constraints to use with built-in DVH object \code{\link[DVHmetrics]{dataMZ}}.}
\usage{data(dataConstr)}
\format{
  A data frame with 6 entries for the following 3 variables.
  \describe{
    \item{\code{constraint}}{The constraint \code{character} string.}
    \item{\code{patID}}{The patient ID \code{character} string or \code{*} wildcard.}
    \item{\code{structure}}{The structure \code{character} string or \code{*} wildcard.}
  }
}
\details{
See \code{\link[DVHmetrics]{checkConstraint}} for the definition of a constraint.
}
\seealso{
\code{\link[DVHmetrics]{readConstraint}}
\code{\link[DVHmetrics]{checkConstraint}}
\code{\link[DVHmetrics]{showConstraint}}
}
\examples{
checkConstraint(dataMZ, constr=dataConstr)
}
\keyword{datasets}
