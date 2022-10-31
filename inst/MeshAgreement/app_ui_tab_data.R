fluidPage(
    fluidRow(
        column(
            width=12,
            box(title="Load 3D mesh data",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                radioButtons("meshes_in",
                             label=NULL,
                             list("Use built-in data"=1,
                                  "Upload 3D mesh files"=2)),
                uiOutput("ui_select_files"),
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
                  "Likely, the mesh does not bound a volume, or has self-intersections. Some agreement measures will be unavailable.",
                  "Maybe, this can be fixed using surface reconstruction on import."),
                htmlOutput("print_mesh_info")
            )
        )
    )
)
