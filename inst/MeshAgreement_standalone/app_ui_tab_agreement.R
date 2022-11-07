fluidPage(
    fluidRow(
        column(
            width=12,
            box(title="Options for agreement calculations",
                width=12,
                # height="650px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                collapsed=FALSE,
                p("Volume-overlap based metrics (DSC, JSC) are available, but take more time to compute than distance-based metrics."),
                checkboxInput("mesh_agree_do_ui", "Calculate DSC, JSC", FALSE)
            )
        )
    ),
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
                p("For the definition of the following options, see the",
                  tags$a(href="https://rdrr.io/github/zarquon42b/Rvcg/man/vcgMetro.html",
                         "documentation of the Rvcg::vcgMetro() function"),
                  "as well the",
                  tags$a(href="http://vcg.isti.cnr.it/vcglib/metro.html",
                         "VCG metro page")),
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
            box(title="Pairwise agreement",
                width=12,
                height="620px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                plotlyOutput("diag_agree_pairwise")
            )
        )
    ),
    fluidRow(
        column(
            width=12,
            box(title="Average agreement",
                width=12,
                height="620px",
                status=NULL,
                closable=FALSE,
                maximizable=FALSE,
                collapsible=FALSE,
                plotlyOutput("diag_agree_aggr")
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
