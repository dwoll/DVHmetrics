\name{get_mesh_info_one}
\alias{get_mesh_info_one}
\title{Get information about one mesh object}
\description{Get key information about one mesh object. This includes the number of vertices and triangles as well as the volume and the centroid.}
\usage{
get_mesh_info_one(x)
}
\arguments{
  \item{x}{A single mesh object as produced by \code{\link[MeshAgreement]{read_mesh_one}}.}
}
\details{Use \code{\link[MeshAgreement]{get_mesh_info}} for getting the information about a list of mesh observer/objects. Use \code{\link[MeshAgreement]{get_mesh_info_obs}} for getting the information about a list of meshes.}
\value{
\code{data.frame} with these variables
\item{name}{Name of the mesh.}
\item{n_verts}{Number of vertices.}
\item{n_faces}{Number of faces.}
\item{volume}{Mesh volume.}
\item{ctr_x}{Centroid x-coordinate.}
\item{ctr_y}{Centroid y-coordinate.}
\item{ctr_z}{Centroid z-coordinate.}
}
\seealso{
\code{\link[MeshAgreement]{read_mesh_one}},
\code{\link[MeshAgreement]{get_mesh_info_obs}},
\code{\link[MeshAgreement]{get_mesh_info}}
}
\examples{
heartL <- mesh3dL_to_cgalMeshL(data_heart_obsL)
get_mesh_info_one(heartL[["Obs01"]][["Obs01_HEART"]])
}
