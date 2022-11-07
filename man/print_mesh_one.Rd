\name{print_mesh_one}
\alias{print_mesh_one}
\title{Print information for one mesh object}
\description{Print the information stored in one mesh object. This includes the number of vertices and triangles as well as the volume and the centroid.}
\usage{
print_mesh_one(x)
}
\arguments{
  \item{x}{A single mesh object as produced by \code{\link[MeshAgreement]{read_mesh_one}}.}
}
\details{Use \code{\link[MeshAgreement]{print_mesh}} for printing information for a list of mesh observer/objects. Use \code{\link[MeshAgreement]{print_mesh_obs}} for printing information for a list of meshes.}
\value{Text output on the console.}
\seealso{
\code{\link[MeshAgreement]{read_mesh_one}},
\code{\link[MeshAgreement]{print_mesh_obs}},
\code{\link[MeshAgreement]{print_mesh}}
}
\examples{
print_mesh_one(data_heart_obsL[[1]][[1]])
}
