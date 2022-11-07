\name{read_mesh_one}
\alias{read_mesh_one}
\title{Read one mesh file}
\description{Read one mesh object from a file.}
\usage{read_mesh_one(x, name,
              reconstruct=c("No", "AFS", "Poisson"),
              spacing=1)}
\arguments{
  \item{x}{\code{character}. File name.}
  \item{name}{\code{character}. Name that should be assigned to the mesh object. If missing, the name is taken from the file name.}
  \item{reconstruct}{\code{character}. Use \code{"No"} for no reconstruction. \code{"AFS"} uses the \code{AFSreconstruction} method of the \code{\link[cgalMeshes]{cgalMesh}} class to re-build the triangles from the vertices. Useful for meshes with defects such as not being closed or having self-intersections. \code{"Poisson"} durrently not implemented and alls back to \code{AFS}.}
  \item{spacing}{Controls vertex spacing for \code{reconstruct="Poisson"}. Currently not implemented.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF.}
\value{Mesh object - a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[cgalMeshes]{cgalMesh}} object. Use \code{\link[MeshAgreement]{read_mesh}} to read a list of several mesh objects each for several observers. Use \code{\link[MeshAgreement]{read_mesh_obs}} to read a list of several mesh objects each for 1 observer.}
\seealso{
\code{\link[cgalMeshes]{cgalMesh}},
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[MeshAgreement]{read_mesh_obs}},
}
\examples{
\dontrun{
ff <- list.files("c:/meshes/", pattern="PLY$", full.names=TRUE)
meshL <- read_mesh_one(ff[1], reconstruct="AFS")
}
}
