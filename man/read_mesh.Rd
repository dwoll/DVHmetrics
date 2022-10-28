\name{read_mesh}
\alias{read_mesh}
\title{Read mesh files}
\description{Read a list of mesh objects from file names.}
\usage{read_mesh(x), name, reconstruct=TRUE}
\arguments{
  \item{x}{\code{character} vector of file names.}
  \item{name}{\code{character} vector of names that should be assigned to the mesh objects. If missing, names are taken from the file names.}
  \item{reconstruct}{Use \code{\link[SurfaceReconstruction]{AFSreconstruction}} to re-build the triangles from the vertices. Useful for meshes with defects.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF.}
\value{List of \code{\link[rgl]{mesh3d}} objects. Use \code{\link[MeshAgreement]{read_mesh_one}} to read a single mesh object.}
\seealso{
\code{\link[rgl]{mesh3d}},
\code{\link[MeshAgreement]{read_mesh_one}},
\code{\link[SurfaceReconstruction]{AFSreconstruction}}
}
\examples{
1
}
