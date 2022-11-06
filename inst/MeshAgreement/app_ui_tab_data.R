fluidPage(
    fluidRow(
        column(
            width=12,
            box(title="Select 3D mesh data",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                radioButtons("meshes_input_source",
                             label="Data source",
                             list("Use built-in data"="builtin",
                                  "Upload 3D mesh files (Supported file formats: STL, PLY, OBJ, OFF)"="file"),
                             inline=TRUE),
                uiOutput("ui_surface_recon_method"),
                uiOutput("ui_surface_recon_spacing"),
                radioButtons("meshes_sel_mode",
                             label="Comparison mode",
                             list("All pairwise comparisons"="all_pairwise",
                                  "Define comparisons individually"="indiv"),
                             inline=TRUE),
                uiOutput("ui_select_comparisons"),
                uiOutput("ui_select_files"),
                actionButton("apply_file_sel", "Apply")
            )
        )
    ),
    fluidRow(
        column(
            width=12,
            box(title="Comparisons",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                uiOutput("ui_ranklist_files"),
                DT::dataTableOutput("ui_compare_table")
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
