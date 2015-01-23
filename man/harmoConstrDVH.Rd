\name{harmoConstrDVH}
\alias{harmoConstrDVH}
\alias{harmoConstrDVH.DVHLst}
\alias{harmoConstrDVH.DVHLstLst}
\title{Harmonize DVH list and constraint list}
\description{
Harmonize DVH list and constraint list such that they apply to the same set of patients and structures.
}
\usage{
harmoConstrDVH(x, constr, byPat=TRUE)

\method{harmoConstrDVH}{DVHLst}(x, constr, byPat=TRUE)

\method{harmoConstrDVH}{DVHLstLst}(x, constr, byPat=TRUE)
}
\arguments{
  \item{x}{Multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{constr}{One or more constraints - given as a \code{character} vector or as a \code{data.frame}. See Details.}
  \item{byPat}{Relevant if multiple DVHs are given. If \code{x} has class \code{DVHLst}: \code{byPat=TRUE} means that the DVHs are for one patient with multiple structures. \code{byPat=FALSE} means that the DVHs are for one structure from multiple patients. If \code{x} has class \code{DVHLstLst}: \code{byPat=TRUE} means that the DVHs are for multiple patients (list components of \code{x}) with multiple structures. \code{byPat=FALSE} means that the DVHs are for multiple structures (list components of \code{x}) from multiple patients.}
}
\value{
A list with those DVHs and constraints that apply to the same patients/structures.
\item{x}{DVH list - object of class \code{DVHLst} or \code{DVHLstLst}}
\item{constr}{Parsed and expanded constraint list}
}
\seealso{
\code{\link[DVHmetrics]{checkConstraint}},
\code{\link[DVHmetrics]{showConstraint}}
}
\examples{
res <- harmoConstrDVH(dataMZ, dataConstr, byPat=TRUE)
res$x
names(res$constr)
}
