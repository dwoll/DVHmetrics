\name{read_mesh}
\alias{read_mesh}
\title{Read mesh files}
\description{Read several mesh objects each for several observers from files.}
\usage{read_mesh(x, name, fix_issues=TRUE,
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
  \item{x}{\code{list} of \code{character} vectors of file names. One vector for the files from one observer.}
  \item{name}{\code{list} of \code{character} vectors of mesh names that should be assigned to the mesh objects. If missing, names are taken from the file names.}
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
\details{Supported file formats: STL, PLY, OBJ, OFF. For subsequent comparison, the first file from observer 1 will be paired with the first file from observer 2, and so on.}
\value{\code{list} (each component representing one observer) of \code{list}s of mesh objects. Each object is a \code{list} with components \code{name}, \code{mesh}, \code{volume}, \code{centroid}, where \code{mesh} is a \code{\link[cgalMeshes]{cgalMesh}} object. Sample object given in \code{\link[MeshAgreement]{data_heart_obsL}}. Use \code{\link[MeshAgreement]{read_mesh_obs}} to read several mesh files from one observer. Use \code{\link[MeshAgreement]{read_mesh_one}} to read a single mesh file.}
\seealso{
\code{\link[cgalMeshes]{cgalMesh}},
\code{\link[cgalMeshes]{AFSreconstruction}},
\code{\link[cgalMeshes]{SSSreconstruction}},
\code{\link[MeshAgreement]{data_heart_obsL}},
\code{\link[MeshAgreement]{read_mesh_obs}},
\code{\link[MeshAgreement]{read_mesh_one}},
}
\examples{
if(interactive()) {
    src_dir <- system.file("extdata", package="MeshAgreement")
    ff1     <- list.files(src_dir, pattern="^Obs01", full.names=TRUE)
    ff2     <- list.files(src_dir, pattern="^Obs02", full.names=TRUE)
    ff3     <- list.files(src_dir, pattern="^Obs03", full.names=TRUE)
    obsL    <- read_mesh(list(Obs01=ff1,
                              Obs02=ff2,
                              Obs03=ff3),
                         reconstruct="AFS")

    get_mesh_info(obsL)
}
}
