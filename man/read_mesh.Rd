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
  \item{reconstruct}{\code{character}. Use \code{"No"} for no reconstruction. \code{"AFS"} uses the \code{AFSreconstruction} method of the \code{\link[cgalMeshes]{cgalMesh}} class to re-build the triangles from the vertices. Useful for meshes with defects such as not being closed or having self-intersections. \code{"Poisson"} durrently not implemented and alls back to \code{AFS}.}
  \item{spacing}{Controls vertex spacing for \code{reconstruct="Poisson"}. Currently not implemented.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF. For subsequent comparison, the first file from observer 1 will be paired with the first file from observer 2, and so on.}
\value{\code{list} (each component representing one observer) of \code{list}s of mesh objects. Each object is a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[cgalMeshes]{cgalMesh}} object. Sample object given in \code{\link[MeshAgreement]{data_heart_obsL}}. Use \code{\link[MeshAgreement]{read_mesh_obs}} to read several mesh files from one observer. Use \code{\link[MeshAgreement]{read_mesh_one}} to read a single mesh file.}
\seealso{
\code{\link[cgalMeshes]{cgalMesh}},
\code{\link[MeshAgreement]{data_heart_obsL}},
\code{\link[MeshAgreement]{read_mesh_obs}},
\code{\link[MeshAgreement]{read_mesh_one}},
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
