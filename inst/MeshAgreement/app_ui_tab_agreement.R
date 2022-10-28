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
                checkboxInput("compare_limit", "Limit samples to save time", TRUE),
                numericInput("n_mesh_samples", "Number of samples", value=1000) #,
                # actionButton("apply_compare", "Compare meshes - takes a while...")
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
                  actionBttn(
                      inputId="btn_home_go_about2",
                      label="About",
                      size="sm",
                      no_outline=TRUE,
                      style = "minimal",
                      color = "primary"#,
                      #icon = icon("info-circle", lib="font-awesome") #, class="fa-lg")
                  ),
                  "for a definition of calculated measures")
            )
        )
    )
)
