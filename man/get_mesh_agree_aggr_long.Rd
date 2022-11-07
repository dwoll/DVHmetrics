\name{get_mesh_agree_aggr_long}
\alias{get_mesh_agree_aggr_long}
\title{Shape aggregated agreement measures to long format}
\description{Utility function to shape data frame with aggregated mesh agreement measures from wide format to long format.}
\usage{
get_mesh_agree_aggr_long(x)
}
\arguments{
  \item{x}{\code{data.frame} with aggregated agreement measures in wide format as returned by \code{\link[MeshAgreement]{get_mesh_agree_aggr}}.}
}
\details{Front-end for \code{\link[stats]{reshape}}.}
\value{\code{data.frame} with the following variables:
\item{group}{The grouping variable indicating which pairs belong to the same structure.}
\item{metric}{The calculated agreement measure.}
\item{statistic}{The aggregate statistic over mesh pairs. This includes the mean, median, standard deviation, and coefficient of variation.}
\item{observed}{The value observed for the statistic.}
}
\seealso{
\code{\link[MeshAgreement]{get_mesh_agree_aggr}},
\code{\link[stats]{reshape}}
}
\examples{
heartL      <- mesh3dL_to_cgalMeshL(data_heart_obsL)
agreeW      <- get_mesh_agree(heartL, do_ui=FALSE)
agree_aggrW <- get_mesh_agree_aggr(agreeW)
get_mesh_agree_aggr_long(agree_aggrW)
}
