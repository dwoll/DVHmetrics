\name{data_heart_meshL}
\alias{data_heart_meshL}
\docType{data}
\title{Heart structure delineated by 3 observers}
\description{3D meshes for the complete heart from one patient, delineated by 3 observers.}
\usage{data(data_heart_meshL)}
\format{Object of class \code{list} with 3 components of class \code{\link[rgl]{mesh3d}}.
  \describe{
    \item{\code{HEART_OBS1}}{Object of class \code{\link[rgl]{mesh3d}} with a mesh of the complete heart, delineated by observer 1.}
    \item{\code{HEART_OBS2}}{Object of class \code{\link[rgl]{mesh3d}} with a mesh of the complete heart, delineated by observer 2.}
    \item{\code{HEART_OBS2}}{Object of class \code{\link[rgl]{mesh3d}} with a mesh of the complete heart, delineated by observer 3.}
  }
}
\details{Data courtesy of Department of Radiation Oncology, University Medical Center Mainz, Germany.}
\seealso{
\code{\link[rgl]{mesh3d}},
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[MeshAgreement]{print_mesh}}
}
\examples{
print_mesh(data_heart_meshL)
}
\keyword{datasets}
