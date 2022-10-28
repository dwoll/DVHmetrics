\name{print_mesh}
\alias{print_mesh}
\title{Print information for a list of mesh objects}
\description{Print the information stored in a list of mesh objects, i.e., their respective number of vertices and triangles as well as their volume and centroid.}
\usage{
print_mesh(x, html=FALSE)
}
\arguments{
  \item{x}{A list of mesh objects as produced by \code{\link[MeshAgreement]{read_mesh}}.}
  \item{html}{\code{logical}. Format output in HTML code as is useful in a Shiny application.}
}
\details{Use \code{\link[MeshAgreement]{print_mesh}} for printing information for a list of mesh objects.}
\value{Text output.}
\seealso{
\code{\link[MeshAgreement]{read_mesh}},
\code{\link[MeshAgreement]{print_mesh_one}}
}
\examples{
print_mesh(data_heart_meshL)
}
