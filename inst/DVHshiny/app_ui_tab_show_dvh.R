fluidPage(
    fluidRow(
        box(title="Plot options",
            width=4,
            radioButtons("plotByPat", label=h5("Plot by patient or by structure"),
                         list("By patient"=1,
                              "By structure"=2)),
            uiOutput("plotSelPat"),
            actionButton("plotSelPatAll", label="(De)Select All"),
            uiOutput("plotSelStruct"),
            actionButton("plotSelStructAll", label="(De)Select All"),
            radioButtons("plotPlotVol", label=h5("Plot relative/absolute volume"),
                         list("Relative volume"=1,
                              "Absolute volume"=2)),
            radioButtons("plotType", label=h5("DVH type"),
                         list("Cumulative"=1,
                              "Differential"=2)),
            checkboxInput("plotMSD", "Show M + SD areas", FALSE),
            sliderInput("plotThreshVol", label=h5("Threshold volume"),
                        min=0, max=100, value=1)
        ),
        box(title="Show DVH diagrams",
            width=8,
            downloadButton("saveDVHPDF", "Save as PDF"),
            downloadButton("saveDVHJPG", "Save as JPEGs (zipped to one file)"),
            #plotOutput("DVHplotOrg"),
            uiOutput("DVHplot"),
            downloadButton("saveDVHMSD", "Save DVH Mean + SD as text file"),
            inputPanel(
                radioButtons("saveDVHDec", "Decimal separator:",
                             list("."=1, ","=2)),
                radioButtons("saveDVHSep", "Column separator:",
                             list("\\t (tab)"=1, "' ' (space)"=2, ", (comma)"=3, "; (semicolon)"=4)))
        )
    )
)
