fluidPage(
    fluidRow(
        bs4Box(
            title="Constraint plot settings",
            width=4,
            radioButtons("constrByPat", label=h5("Plot by patient or by structure"),
                         list("By patient"=1,
                              "By structure"=2)),
            radioButtons("constrPlotVol", label=h5("Plot relative/absolute volume"),
                         list("Relative volume"=1,
                              "Absolute volume"=2)),
            sliderInput("constrThreshVol", label=h5("Threshold volume"),
                        min=0, max=100, value=1)
        ),
        bs4Box(
            title="Constraint plot",
            width=8,
            downloadButton("saveConstrPDF", "Save as PDF"),
            downloadButton("saveConstrJPG", "Save as JPEGs (zipped to one file)"),
            #plotOutput("constraintPlotOrg"),
            uiOutput("constraintPlot")
        )
    )
)
