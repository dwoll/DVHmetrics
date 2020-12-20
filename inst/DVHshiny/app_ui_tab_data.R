fluidPage(
    fluidRow(
        box(title="Enter data",
            width=4,
            radioButtons("DVHin", label=h4("Enter data"),
                         list("Use built-in data"=1,
                              "Upload DVH file"=2)),
            conditionalPanel(condition="input.DVHin == '2'",
                             h5("Upload DVH file: "),
                             radioButtons("DVHtype", "DVH file format:",
                                          list("Eclipse"=1, "CadPlan"=2, "MasterPlan"=3,
                                               "Pinnacle3"=4, "Monaco"=5, "HiArt"=6,
                                               "RayStation"=7, "ProSoma"=8, "PRIMO"=9)),
                             fileInput("DVHupload", "Select DVH file:", multiple=TRUE),
                             # radioButtons("fileenc", "File encoding:",
                             #              list("Default"=1, "UTF-8"=2, "UTF-8-BOM"=3)),
                             radioButtons("DVHplanInfo", "Information encoded in plan:",
                                          list("None"=1, "Prescribed dose"=2)),
                             checkboxGroupInput("DVHreadOpts", label=NULL,
                                                choices=c("Add to existing data"="DVHadd",
                                                          "Use Course for ID"="DVHcourse",
                                                          "Struct volume from DVH"="volume_from_dvh"))),
            actionButton("applyData", "Apply"),
            radioButtons("DVHverbose", "",
                         list("Short info on DVHs"=1,
                              "Detailed info on DVHs"=2))
        ),
        box(title="Information from imported DVH file(s)",
            width=8,
            verbatimTextOutput("DVHinfo")
        )
    )
)
