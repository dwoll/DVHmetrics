\name{print_mesh_one}
\alias{print_mesh_one}
\title{Print information for one mesh object}
\description{Print the information stored in one mesh object, i.e., the number of vertices and triangles as well as the volume and the centroid.}
\usage{
print_mesh_one(x, html=FALSE)
}
\arguments{
  \item{x}{A single mesh object as produced by \code{\link[MeshAgreement]{read_mesh_one}}.}
  \item{html}{\code{logical}. Format output in HTML code as is useful in a Shiny application.}
}
\details{Use \code{\link[MeshAgreement]{print_mesh}} for printing information for a list of mesh objects.}
\value{Text output.}
\seealso{
\code{\link[MeshAgreement]{read_mesh_one}},
\code{\link[MeshAgreement]{print_mesh}}
}
\examples{
print_mesh_one(data_heart_meshL[[1]])
}
