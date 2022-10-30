\name{read_mesh}
\alias{read_mesh}
\title{Read mesh files}
\description{Read a list of mesh objects from file names.}
\usage{read_mesh(x, name,
          reconstruct=c("none", "AFS", "Poisson"),
          spacing)}
\arguments{
  \item{x}{\code{character} vector of file names.}
  \item{name}{\code{character} vector of names that should be assigned to the mesh objects. If missing, names are taken from the file names.}
  \item{reconstruct}{\code{character}. Use \code{"none"} for no reconstruction. \code{"AFS"} uses \code{\link[SurfaceReconstruction]{AFSreconstruction}}, and\code{"Poisson"} uses \code{\link[SurfaceReconstruction]{PoissonReconstruction}} to re-build the triangles from the vertices. Useful for meshes with defects such as not being closed or having self-intersections.}
  \item{spacing}{Controls vertex spacing for \code{reconstruct="Poisson"}, see \code{\link[SurfaceReconstruction]{PoissonReconstruction}}.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF.}
\value{List of mesh objects. Each object is a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[rgl]{mesh3d}} object. Use \code{\link[MeshAgreement]{read_mesh_one}} to read a single mesh file.}
\seealso{
\code{\link[rgl]{mesh3d}},
\code{\link[MeshAgreement]{read_mesh_one}},
\code{\link[SurfaceReconstruction]{AFSreconstruction}},
\code{\link[SurfaceReconstruction]{PoissonReconstruction}}
}
\examples{
\dontrun{
f <- list.files("c:/meshes/", pattern="PLY$", full.names=TRUE)
meshL <- read_mesh(f, reconstruct="AFS")
}
}
