\name{read_mesh}
\alias{read_mesh}
\title{Read mesh files}
\description{Read several mesh objects each for several observers from files.}
\usage{read_mesh(x, name,
          reconstruct=c("No", "AFS", "Poisson"),
          spacing=1)}
\arguments{
  \item{x}{\code{list} of \code{character} vector of file names. One vector for the files from one observer.}
  \item{name}{\code{character} vector of names that should be assigned to the mesh objects. If missing, names are taken from the file names.}
  \item{reconstruct}{\code{character}. Use \code{"No"} for no reconstruction. \code{"AFS"} uses \code{\link[SurfaceReconstruction]{AFSreconstruction}}, and\code{"Poisson"} uses \code{\link[SurfaceReconstruction]{PoissonReconstruction}} to re-build the triangles from the vertices. Useful for meshes with defects such as not being closed or having self-intersections.}
  \item{spacing}{Controls vertex spacing for \code{reconstruct="Poisson"}, see \code{\link[SurfaceReconstruction]{PoissonReconstruction}}.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF. For subsequent comparison, the first file from observer 1 will be paired with the first file from observer 2, and so on.}
\value{\code{list} (each component representing one observer) of \code{list}s of mesh objects. Each object is a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[rgl]{mesh3d}} object. Sample object given in \code{\link[MeshAgreement]{data_heart_obsL}}. Use \code{\link[MeshAgreement]{read_mesh_obs}} to read several mesh files from one observer. Use \code{\link[MeshAgreement]{read_mesh_one}} to read a single mesh file.}
\seealso{
\code{\link[rgl]{mesh3d}},
\code{\link[MeshAgreement]{data_heart_obsL}},
\code{\link[MeshAgreement]{read_mesh_obs}},
\code{\link[MeshAgreement]{read_mesh_one}},
\code{\link[SurfaceReconstruction]{AFSreconstruction}},
\code{\link[SurfaceReconstruction]{PoissonReconstruction}}
}
\examples{
\dontrun{
ff1 <- list.files("c:/meshes/obs1", pattern="PLY$", full.names=TRUE)
ff2 <- list.files("c:/meshes/obs1", pattern="PLY$", full.names=TRUE)
ff3 <- list.files("c:/meshes/obs1", pattern="PLY$", full.names=TRUE)
obsL <- read_mesh(list(Obs01=ff1, Obs02=ff2, Obs03=ff3),
                  reconstruct="AFS")
}
}
