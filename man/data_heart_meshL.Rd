\name{data_heart_meshL}
\alias{data_heart_meshL}
\docType{data}
\title{Heart structure delineated by 3 observers}
\description{3D meshes for the complete heart from one patient, as delineated by 3 observers. No gold standard.}
\usage{data(data_heart_meshL)}
\format{Object of class \code{list} with 3 components, each representing the complete heart as delineated by one observer.
  \describe{
    \item{\code{HEART_OBS1}}{Mesh object - a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[rgl]{mesh3d}} object.}
    \item{\code{HEART_OBS2}}{Mesh object with the same structure as \code{HEART_OBS1}.}
    \item{\code{HEART_OBS2}}{Mesh object with the same structure as \code{HEART_OBS1}.}
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
