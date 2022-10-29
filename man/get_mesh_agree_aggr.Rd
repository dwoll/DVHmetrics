\name{get_mesh_agree_aggr}
\alias{get_mesh_agree_aggr}
\title{Aggregate agreement measures for 3D mesh pairs}
\description{Aggregates the agreement measures calculated for pairs of 3D meshes.}
\usage{
get_mesh_agree_aggr(x, na.rm=FALSE)
}
\arguments{
  \item{x}{\code{data.frame} with agreement measures as returned by \code{\link[MeshAgreement]{get_mesh_agree}}.}
  \item{na.rm}{\code{logical}. Remove missings before calculating statistical measures?}
}
\value{\code{data.frame} with the following variables:
\item{metric}{The calculated agreement measure.}
\item{Mean}{The arithmetic mean over mesh pairs.}
\item{Median}{The median over mesh pairs.}
\item{SD}{The standard deviation over mesh pairs.}
\item{CV}{The coefficient of variation calculated as SD/M.}
\item{CV_ln}{The coefficient of variation calculated assuming log-normality.}
}
\seealso{
\code{\link[MeshAgreement]{get_mesh_agree}}
}
\examples{
d_agree_pair <- get_mesh_agree(data_heart_meshL)
get_mesh_agree_aggr(d_agree_pair)
}
