\name{read_mesh_one}
\alias{read_mesh_one}
\title{Read one mesh file}
\description{Read one mesh object from a file name.}
\usage{read_mesh_one(x, name, reconstruct=TRUE)}
\arguments{
  \item{x}{\code{character}. File name.}
  \item{name}{\code{character}. Name that should be assigned to the mesh object. If missing, the name is taken from the file name.}
  \item{reconstruct}{\code{logical}. Use \code{\link[SurfaceReconstruction]{AFSreconstruction}} to re-build the triangles from the vertices. Useful for meshes with defects.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF.}
\value{Mesh object - a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[rgl]{mesh3d}} object. Use \code{\link[MeshAgreement]{read_mesh}} to read a list of mesh objects.}
\seealso{
\code{\link[rgl]{mesh3d}},
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[SurfaceReconstruction]{AFSreconstruction}}
}
\examples{
\dontrun{
f <- list.files("c:/meshes/", pattern="PLY$", full.names=TRUE)
meshL <- read_mesh_one(f[1])
}
}
