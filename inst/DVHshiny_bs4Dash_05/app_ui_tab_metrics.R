fluidPage(
    fluidRow(
        bs4Box(
            title="Define metrics",
            width=4,
            textInput("metrInput", "Metric(s):", value=c("DMEAN, D1cc, V10%")),
            #tags$textarea(id="defMetricsMult", rows=2, cols=10, ""),
            #actionButton('clearText_button','Clear metrics'),
            #radioButtons("metrInterp", label=h5("DVH interpolation"),
            #             list("Linear"=1,
            #                  "Monotone spline"=2,
            #                  "Local polynomial"=3)),
            checkboxInput("metrEUDparam", "Show EUD params ...", FALSE),
            conditionalPanel(condition="input.metrEUDparam == true",
                             textInput("metrEUDa", h5("exponent a"), value=""),
                             textInput("metrEUDfd", h5("fraction dose"), value=""),
                             textInput("metrEUDab", h5("alpha/beta ratio"), value="")
            ),
            checkboxInput("metrNTCPparam", "Show (N)TCP params ...", FALSE),
            conditionalPanel(condition="input.metrNTCPparam == true",
                             radioButtons("metrNTCPtype", h5("(N)TCP Model"),
                                          list("Probit (Lyman KB)"=1,
                                               "Logit (Niemierko)"=2,
                                               "Probit (Kaellman)"=3)),
                             textInput("metrNTCPtd50", h5("T(C)D50"), value=""),
                             textInput("metrNTCPn", h5("n (=1 / EUD-a)"), value=""),
                             conditionalPanel(condition="input.metrNTCPtype == '1'",
                                              textInput("metrNTCPm", h5("Lyman m"), value="")
                             ),
                             conditionalPanel(condition="input.metrNTCPtype != '1'",
                                              textInput("metrNTCPgamma50", h5("Logit/Poisson gamma50"), value="")
                             )
            ),
            uiOutput("metrSelPat"),
            actionButton("metrSelPatAll", label="(De)Select All"),
            uiOutput("metrSelStruct"),
            actionButton("metrSelStructAll", label="(De)Select All"),
            selectizeInput("metrSortBy", label=h5("Sort output table by:"),
                           choices=c("Value"=1,
                                     "Structure"=2,
                                     "Metric"=3,
                                     "Patient ID"=4),
                           multiple=TRUE)#,
            #options=c(placeholder='Click to select variables'))
        ),
        bs4Box(
            title="DVH metrics",
            width=8,
            DT::dataTableOutput("metrics"),
            downloadButton("saveMetrics", "Save as text file"),
            inputPanel(
                radioButtons("saveMetrDec", "Decimal separator:",
                             list("."=1, ","=2)),
                radioButtons("saveMetrSep", "Column separator:",
                             list("\\t (tab)"=1, "' ' (space)"=2, ", (comma)"=3, "; (semicolon)"=4)))
        )
    )
)
