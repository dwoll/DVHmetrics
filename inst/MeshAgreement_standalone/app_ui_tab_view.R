fluidPage(
    fluidRow(
        column(
            width=12,
            box(title="Settings",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                uiOutput("rgl_view_selection")
            )
        )
    ),
    fluidRow(
        column(
            width=6,
            box(title=htmlOutput("rgl_mesh1_name"),
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                rglwidgetOutput("rgl_mesh1",
                                width="512px",
                                height="512px")
            )
        ),
        column(
            width=6,
            box(title=htmlOutput("rgl_mesh2_name"),
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                rglwidgetOutput("rgl_mesh2",
                                width="512px",
                                height="512px")
            )
        )
    ),
    # fluidRow(
    #     column(
    #         width=6,
    #         box(title="Union",
    #             width=12,
    #             # height="650px",
    #             status=NULL,
    #             closable=FALSE,
    #             maximizable=FALSE,
    #             collapsible=FALSE,
    #             rglwidgetOutput("rgl_mesh_union",
    #                             width="512px",
    #                             height="512px")
    #         )
    #     ),
    #     column(
    #         width=6,
    #         box(title="Intersection",
    #             width=12,
    #             # height="650px",
    #             status=NULL,
    #             closable=FALSE,
    #             maximizable=FALSE,
    #             collapsible=FALSE,
    #             rglwidgetOutput("rgl_mesh_intersection",
    #                             width="512px",
    #                             height="512px")
    #         )
    #     )
    # ),
    fluidRow(
        column(
            width=6,
            box(title=htmlOutput("rgl_dist1_name"),
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                rglwidgetOutput("rgl_mesh_dist1",
                                width="512px",
                                height="512px")
            )
        ),
        column(
            width=6,
            box(title=htmlOutput("rgl_dist2_name"),
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                rglwidgetOutput("rgl_mesh_dist2",
                                width="512px",
                                height="512px")
            )
        )
    )
)
