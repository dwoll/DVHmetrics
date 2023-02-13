\name{read_mesh_obs}
\alias{read_mesh_obs}
\title{Read mesh files}
\description{Read a list of mesh objects from one observer from files.}
\usage{read_mesh_obs(x, name, fix_issues=TRUE,
              reconstruct=c("No", "AFS", "SSS", "Poisson"),
              jetSmoothing=NULL,
              scaleIterations=1,
              neighbors=12,
              samples=300,
              separateShells=FALSE,
              forceManifold=TRUE,
              borderAngle=45,
              normals="jet(12)",
              spacing="ave(12)",
              sm_angle=20,
              sm_radius=30,
              sm_distance=0.375)}
\arguments{
  \item{x}{\code{character} vector of file names.}
  \item{name}{\code{character} vector of names that should be assigned to the mesh objects. If missing, names are taken from the file names.}
  \item{fix_issues}{\code{logical}. Try to fix mesh issues, such as not being closed, not bounding a volume, or having self-intersections. Uses surface reconstruction if mesh is not closed.}
  \item{reconstruct}{\code{character}. Use \code{"No"} for no reconstruction. \code{"AFS"} uses \code{\link[cgalMeshes]{AFSreconstruction}}. \code{"SSS"} uses \code{\link[cgalMeshes]{SSSreconstruction}}. \code{"Poisson"} uses \code{PoissonReconstruction} if the development branch of \code{cgalMeshes} is installed, and falls back to \code{AFS} otherwise.}
  \item{jetSmoothing}{Used with \code{reconstruct="AFS"}. See \code{\link[cgalMeshes]{AFSreconstruction}}.}
  \item{scaleIterations}{Used with \code{reconstruct="SSS"}. See \code{\link[cgalMeshes]{SSSreconstruction}}.}
  \item{neighbors}{Used with \code{reconstruct="SSS"}. See \code{\link[cgalMeshes]{SSSreconstruction}}.}
  \item{samples}{Used with \code{reconstruct="SSS"}. See \code{\link[cgalMeshes]{SSSreconstruction}}.}
  \item{separateShells}{Used with \code{reconstruct="SSS"}. See \code{\link[cgalMeshes]{SSSreconstruction}}.}
  \item{forceManifold}{Used with \code{reconstruct="SSS"}. See \code{\link[cgalMeshes]{SSSreconstruction}}.}
  \item{borderAngle}{Used with \code{reconstruct="SSS"}. See \code{\link[cgalMeshes]{SSSreconstruction}}.}
  \item{normals}{Used with \code{reconstruct="Poisson"}. Currently only implemented in development branch of \code{cgalMeshes} which has function \code{PoissonReconstruction}.}
  \item{spacing}{Used with \code{reconstruct="Poisson"}. Currently only implemented in development branch of \code{cgalMeshes} which has function \code{PoissonReconstruction}.}
  \item{sm_angle}{Used with \code{reconstruct="Poisson"}. Currently only implemented in development branch of \code{cgalMeshes} which has function \code{PoissonReconstruction}.}
  \item{sm_radius}{Used with \code{reconstruct="Poisson"}. Currently only implemented in development branch of \code{cgalMeshes} which has function \code{PoissonReconstruction}.}
  \item{sm_distance}{Used with \code{reconstruct="Poisson"}. Currently only implemented in development branch of \code{cgalMeshes} which has function \code{PoissonReconstruction}.}
}
\details{Supported file formats: STL, PLY, OBJ, OFF.}
\value{List of mesh objects. Each object is a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[cgalMeshes]{cgalMesh}} object. Use \code{\link[MeshAgreement]{read_mesh}} to read a list of several mesh objects each for several observers. Use \code{\link[MeshAgreement]{read_mesh_one}} to read a single mesh file.}
\seealso{
\code{\link[cgalMeshes]{cgalMesh}},
\code{\link[cgalMeshes]{AFSreconstruction}},
\code{\link[cgalMeshes]{SSSreconstruction}},
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
