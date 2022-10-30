fluidPage(
    fluidRow(
        column(
            width=12,
            box(title="Options for VCG metro distance calculations",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=TRUE,
                collapsed=TRUE,
                uiOutput("ui_mesh_agree_metro_options")
            )
        )
    ),
    fluidRow(
        column(
            width=12,
            box(title="Pairwise agreement",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                DT::dataTableOutput("table_agree_pairwise")
            )
        )
    ),
    fluidRow(
        column(
            width=12,
            box(title="Average agreement",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                DT::dataTableOutput("table_agree_aggr")
            )
        )
    ),
    fluidRow(
        column(
            width=12,
            box(title="Calculated agreement measures",
                width=12,
                p("See section",
                  actionButton(inputId="bttn_home_go_about2",
                               label="About"),
                  "for a definition of calculated measures")
            )
        )
    )
)
