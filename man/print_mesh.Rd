\name{print_mesh}
\alias{print_mesh}
\title{Print information for a list of mesh objects}
\description{Print the information stored in a list of mesh objects. This includes their respective number of vertices and triangles as well as their volume and centroid.}
\usage{
print_mesh(x)
}
\arguments{
  \item{x}{A list of mesh objects as produced by \code{\link[MeshAgreement]{read_mesh}}.}
}
\details{Use \code{\link[MeshAgreement]{print_mesh_one}} for printing information for a single mesh object.}
\value{Text output on the console.}
\seealso{
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[MeshAgreement]{print_mesh_one}}
}
\examples{
print_mesh(data_heart_meshL)
}
