\name{read_mesh_obs}
\alias{read_mesh_obs}
\title{Read mesh files}
\description{Read a list of mesh objects from one observer from files.}
\usage{read_mesh_obs(x, name, fix_issues=TRUE,
              reconstruct=c("No", "AFS", "Poisson"),
              spacing=1)}
\arguments{
  \item{x}{\code{character} vector of file names.}
  \item{name}{\code{character} vector of names that should be assigned to the mesh objects. If missing, names are taken from the file names.}
  \item{fix_issues}{\code{logical}. Try to fix mesh issues, such as not being closed, not bounding a volume, or having self-intersections. Uses surface reconstruction if mesh is not closed.}
  \item{reconstruct}{\code{character}. Use \code{"No"} for no reconstruction. \code{"AFS"} uses the \code{AFSreconstruction} method of the \code{\link[cgalMeshes]{cgalMesh}} class to re-build the triangles from the vertices. Useful for meshes with defects such as not being closed or having self-intersections. \code{"Poisson"} currently not implemented and falls back to \code{AFS}.}
  \item{spacing}{Controls vertex spacing for \code{reconstruct="Poisson"}. Currently not implemented.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF.}
\value{List of mesh objects. Each object is a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[cgalMeshes]{cgalMesh}} object. Use \code{\link[MeshAgreement]{read_mesh}} to read a list of several mesh objects each for several observers. Use \code{\link[MeshAgreement]{read_mesh_one}} to read a single mesh file.}
\seealso{
\code{\link[cgalMeshes]{cgalMesh}},
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[MeshAgreement]{read_mesh_one}},
}
\examples{
if(interactive()) {
    src_dir <- system.file("extdata", package="MeshAgreement")
    ff      <- list.files(src_dir, pattern="^Obs01", full.names=TRUE)
    meshL   <- read_mesh_obs(ff, reconstruct="AFS")
    get_mesh_info_obs(meshL)
}
}
