fluidPage(
    fluidRow(
        box(title="BED/EQD2 settings",
            width=4,
            h4("Conversion type"),
            selectInput("BEDtype", label=NULL,
                        choices=list("BED"=1,
                                     "EQD2"=2,
                                     "Isoeffective dose"=3),
                        selected=1),
            conditionalPanel(condition="(input.BEDtype == '1') || (input.BEDtype == '2')",
                             h4("Input"),
                             textInput("BED_BED_D", h5("Total Dose"),
                                       value=c("50")),
                             textInput("BED_BED_FD", h5("Fractional Dose"),
                                       value=c("1.5 2 2.5")),
                             textInput("BED_BED_FN", h5("Number of fractions"),
                                       value=""),
                             textInput("BED_BED_AB", h5("alpha/beta ratio"),
                                       value="2")
            ),
            conditionalPanel(condition="input.BEDtype == '3'",
                             h4("Input"),
                             textInput("BED_IED_D1", h5("Total Dose 1"),
                                       value=c("50")),
                             textInput("BED_IED_D2", h5("Total Dose 2"),
                                       value=""),
                             textInput("BED_IED_FD1", h5("Fractional Dose 1"),
                                       value=c("1.5 2 2.5")),
                             textInput("BED_IED_FD2", h5("Fractional Dose 2"),
                                       value=c("2")),
                             textInput("BED_IED_FN1", h5("Number of fractions 1"),
                                       value=""),
                             textInput("BED_IED_FN2", h5("Number of fractions 2"),
                                       value=""),
                             textInput("BED_IED_AB", h5("alpha/beta ratio"),
                                       value="2")
            )
        ),
        box(title="BED / EQD2 / Isoeffective dose calculation",
            width=8,
            conditionalPanel(condition="input.BEDtype == '1'",
                             h5("BED")
            ),
            conditionalPanel(condition="input.BEDtype == '2'",
                             h5("EQD2")
            ),
            conditionalPanel(condition="input.BEDtype == '3'",
                             h5("Isoeffective dose")
            ),
            
            verbatimTextOutput("BED")
        )
    )
)
