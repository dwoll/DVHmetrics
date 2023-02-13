\name{read_mesh_one}
\alias{read_mesh_one}
\title{Read one mesh file}
\description{Read one mesh object from a file.}
\usage{read_mesh_one(x, name, fix_issues=TRUE,
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
  \item{x}{\code{character}. File name.}
  \item{name}{\code{character}. Name that should be assigned to the mesh object. If missing, the name is taken from the file name.}
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
\value{Mesh object - a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[cgalMeshes]{cgalMesh}} object. Use \code{\link[MeshAgreement]{read_mesh}} to read a list of several mesh objects each for several observers. Use \code{\link[MeshAgreement]{read_mesh_obs}} to read a list of several mesh objects each for 1 observer.}
\seealso{
\code{\link[cgalMeshes]{cgalMesh}},
\code{\link[cgalMeshes]{AFSreconstruction}},
\code{\link[cgalMeshes]{SSSreconstruction}},
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[MeshAgreement]{read_mesh_obs}},
}
\examples{
if(interactive()) {
    src_dir <- system.file("extdata", package="MeshAgreement")
    ff      <- list.files(src_dir, pattern="^Obs01", full.names=TRUE)
    mesh    <- read_mesh_one(ff[1], reconstruct="AFS")
    get_mesh_info_one(mesh)
}
}
