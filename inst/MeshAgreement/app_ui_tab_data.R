fluidPage(
    fluidRow(
        column(
            width=12,
            box(title="Load meshes",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                fileInput("file_sel", "Select files: (supported file formats: STL, PLY, OBJ, OFF)",
                          width="100%", multiple=TRUE),
                actionButton("apply_file_sel", "Apply")
            )
        )
    ),
    fluidRow(
        column(
            width=12,
            box(title="Information on selected meshes",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                p("Note: If information on mesh volume and mesh centroid are missing ('NA'), the mesh is not proper.",
                  "Likely, the mesh does not bound a volume. Some agreement measures will be unavailable."),
                htmlOutput("print_mesh_info")
            )
        )
    )
)
