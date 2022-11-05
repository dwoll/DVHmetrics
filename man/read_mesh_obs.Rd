\name{read_mesh_obs}
\alias{read_mesh_obs}
\title{Read mesh files}
\description{Read a list of mesh objects from one observer from files.}
\usage{read_mesh_obs(x, name,
              reconstruct=c("No", "AFS", "Poisson"),
              spacing=1)}
\arguments{
  \item{x}{\code{character} vector of file names.}
  \item{name}{\code{character} vector of names that should be assigned to the mesh objects. If missing, names are taken from the file names.}
  \item{reconstruct}{\code{character}. Use \code{"No"} for no reconstruction. \code{"AFS"} uses \code{\link[SurfaceReconstruction]{AFSreconstruction}}, and\code{"Poisson"} uses \code{\link[SurfaceReconstruction]{PoissonReconstruction}} to re-build the triangles from the vertices. Useful for meshes with defects such as not being closed or having self-intersections.}
  \item{spacing}{Controls vertex spacing for \code{reconstruct="Poisson"}, see \code{\link[SurfaceReconstruction]{PoissonReconstruction}}.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF.}
\value{List of mesh objects. Each object is a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[rgl]{mesh3d}} object. Use \code{\link[MeshAgreement]{read_mesh}} to read a list of several mesh objects each for several observers. Use \code{\link[MeshAgreement]{read_mesh_one}} to read a single mesh file.}
\seealso{
\code{\link[rgl]{mesh3d}},
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[MeshAgreement]{read_mesh_one}},
\code{\link[SurfaceReconstruction]{AFSreconstruction}},
\code{\link[SurfaceReconstruction]{PoissonReconstruction}}
}
\examples{
\dontrun{
ff <- list.files("c:/meshes/", pattern="PLY$", full.names=TRUE)
meshL <- read_mesh_obs(ff, reconstruct="AFS")
}
}
