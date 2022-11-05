\name{data_heart_obsL}
\alias{data_heart_obsL}
\docType{data}
\title{Heart structures delineated by 3 observers}
\description{3D meshes for the complete heart and for a surrogate volume of the aortic valve from one patient, as delineated by 3 observers. No gold standard.}
\usage{data(data_heart_obsL)}
\format{Object of class \code{list} with 3 components (\code{Obs01}, \code{Obs02}, \code{Obs03}), each representing one observer. Each of these lists contains two meshes - the complete heart, and a surrogate volume for the aortic valve. For \code{Obs01}, the list has this structure:
  \describe{
    \item{\code{Obs01_HEART}}{Mesh object - a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[rgl]{mesh3d}} object.}
    \item{\code{Obs01_AOKL}}{Mesh object with the same structure as \code{Obs01_HEART}.}
  }
}
\details{Data courtesy of Department of Radiation Oncology, University Medical Center Mainz, Germany.}
\seealso{
\code{\link[rgl]{mesh3d}},
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[MeshAgreement]{print_mesh}}
}
\examples{
print_mesh(data_heart_obsL)
}
\keyword{datasets}
