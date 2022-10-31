#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
## Daniel Wollschläger <wollschlaeger@uni-mainz.de>
## Agreement measures for 3D structures
## Pairwise distance-based and volume-overlap-based metrics
## DCOM:   Distance between centers of mass
## ASD:    Average surface distance
## RMSD:   Root mean squared surface distance
## HD_max: Hausdorff distance - max of both directed HDs
## HD_avg: Hausdorff distance - average of both directed HDs
## JSC:    Jaccard similarity coefficient
## DSC:    Dice similarity coefficient
#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------

library(shiny)
library(bs4Dash)
library(rgl)
library(MeshAgreement)
# library(Rvcg)
# library(SurfaceReconstruction)
# library(PolygonSoup)
# library(MeshesTools)
# library(Boov)

source("app_00_global.R")

#####---------------------------------------------------------------------------
## app code
#####---------------------------------------------------------------------------

shiny::shinyApp(
    ui=dashboardPage(
        title="Agreement measures for 3D structures",
        dark=NULL,
        help=FALSE,
        fullscreen=TRUE,
        scrollToTop=TRUE,
        header=dashboardHeader(
            title=NULL,
            fixed=FALSE,
            leftUI=tagList(tags$code(tags$h3("Agreement measures for 3D structures",
                                             style="display: flex; align-items: center; justify-content: center; margin-right: 20px;"))),
            rightUI=tagList()
        ),
        # controlbar=source("app_ui_controlbar.R", encoding="UTF8")$value,
        footer=dashboardFooter(
            fixed=FALSE,
            left=NULL,
            right=tagList(p(actionButton(
                inputId="bttn_footer_impressum",
                label="Impressum / Haftung / Urheberrecht",
                size="sm",
                no_outline=TRUE,
                style="minimal",
                color="primary")
            ))
        ),
        body=dashboardBody(
            tags$style(HTML(".irs-bar         { border-color: transparent; background-color: transparent; }",
                            ".irs-bar-edge    { border-color: transparent; background-color: transparent; }",
                            ".irs-bar--single { border-color: transparent; background-color: transparent; }")),
            tabItems(
                tabItem(
                    tabName="tab_home",
                    source("app_ui_tab_home.R", local=TRUE, encoding="UTF8")$value
                ),
                tabItem(
                    tabName="tab_data",
                    source("app_ui_tab_data.R", local=TRUE, encoding="UTF8")$value
                ),
                tabItem(
                    tabName="tab_agreement",
                    source("app_ui_tab_agreement.R", local=TRUE, encoding="UTF8")$value
                ),
                tabItem(
                    tabName="tab_view",
                    source("app_ui_tab_view.R", local=TRUE, encoding="UTF8")$value
                ),
                tabItem(
                    tabName="tab_about",
                    source("app_ui_tab_about.R", local=TRUE, encoding="UTF8")$value
                )
            )
        ),
        sidebar=source("app_ui_sidebar.R", encoding="UTF8")$value
    ),
    #####-----------------------------------------------------------------------
    ## server
    #####-----------------------------------------------------------------------
    server=function(input, output, session) {
        observeEvent(input$bttn_home_go_about1, {
            updateTabItems(session, inputId="sidebar_tabs", selected="tab_about")
        })
        observeEvent(input$bttn_home_go_about2, {
            updateTabItems(session, inputId="sidebar_tabs", selected="tab_about")
        })
        observeEvent(input$bttn_footer_impressum, {
            showModal(popup_info_impressum)
        })
        react_file_sel <- reactive({
            ## only change when explicitly applied
            input$apply_file_sel
            ## isolate against non-applied changes in data input UI elements
            isolate({
                if(input$meshes_in == 1) {
                    data_heart_meshL
                } else {
                    if(!is.null(input$file_sel)) {
                        read_mesh(input$file_sel$datapath,
                                  input$file_sel$name,
                                  reconstruct=input$read_mesh_reconstruct,
                                  spacing=input$read_mesh_reconstruct_pois_spacing)
                    } else {
                        NULL
                    }
                }
            })
        })
        react_mesh_ui <- reactive({
            meshL <- react_file_sel()
            if(!is.null(meshL)) {
                get_mesh_ui(meshL)
            } else {
                NULL
            }
        })
        react_mesh_metro <- reactive({
            meshL <- react_file_sel()
            if(!is.null(meshL)) {
                argL <- list(meshL,
                             chop        =TRUE,
                             silent      =TRUE,
                             colormeshes =TRUE,
                             nSamples    =input$vcgMetro_nSamples,
                             nSamplesArea=input$vcgMetro_nSamplesArea,
                             vertSamp    =input$vcgMetro_vertSamp,
                             edgeSamp    =input$vcgMetro_edgeSamp,
                             faceSamp    =input$vcgMetro_faceSamp,
                             unrefVert   =input$vcgMetro_unrefVert,
                             samplingTyp =input$vcgMetro_samplingTyp,
                             searchStruct=input$vcgMetro_searchStruct,
                             from        =input$vcgMetro_from,
                             to          =input$vcgMetro_to)
                
                do.call("get_mesh_metro", Filter(Negate(is.null), argL))
            } else {
                NULL
            }
        })
        ## list of pairwise meshes
        react_mesh_agree <- reactive({
            meshL  <- react_file_sel()
            uiL    <- react_mesh_ui()
            metroL <- react_mesh_metro()
            if(!is.null(meshL) && !is.null(uiL) && !is.null(metroL)) {
                mesh_pairL  <- get_mesh_pairs(meshL)
                agree_pairL <- Map(get_mesh_agree_pair,
                                   mesh_pairL,
                                   metro=metroL,
                                   ui   =uiL)
                
                d <- do.call("rbind", agree_pairL)
                rownames(d) <- NULL
                d
            } else {
                NULL
            }
        })
        output$ui_select_files <- renderUI({
            if(input$meshes_in == 2) {
                tagList(fileInput("file_sel",
                                  "Select files: (supported file formats: STL, PLY, OBJ, OFF)",
                                  width="100%", multiple=TRUE),
                        radioButtons("read_mesh_reconstruct",
                                     "Surface reconstruction on import",
                                     choices=c("none", "AFS", "Poisson"),
                                     selected="none"),
                        numericInput("read_mesh_reconstruct_pois_spacing",
                                     "Spacing parameter for Poisson reconstruction",
                                     value=0.05,
                                     min=0))
            } else {
                NULL
            }
        })
        output$print_mesh_info <- renderUI({
            meshL <- react_file_sel()
            if(!is.null(meshL)) {
                print_mesh_html(meshL)
            } else {
                NULL
            }
        })
        output$ui_mesh_agree_metro_options <- renderUI({
            tagList( numericInput("vcgMetro_nSamples",     "Number of samples", value=0L),
                     numericInput("vcgMetro_nSamplesArea", "Number of samples per area (overrides nSamples)", value=0L),
                    checkboxInput("vcgMetro_vertSamp",     "Vertex sampling",            value=TRUE),
                    checkboxInput("vcgMetro_edgeSamp",     "Edge sampling",              value=TRUE),
                    checkboxInput("vcgMetro_faceSamp",     "Face sampling",              value=TRUE),
                    checkboxInput("vcgMetro_unrefVert",    "Ignore unreferred vertices", value=FALSE),
                     radioButtons("vcgMetro_samplingTyp",  "Face sampling mode", choices=c("SS", "MC", "SD"), selected="SS"),
                     radioButtons("vcgMetro_searchStruct", "Search structure",   choices=c("SGRID", "AABB", "OCTREE", "HGRID"), selected="SGRID"),
                     numericInput("vcgMetro_from",         "Color mapping: minimum", value=0L),
                     numericInput("vcgMetro_to",           "Color mapping: maximum", value=0L)
            )
        })
        output$rgl_view_selection <- renderUI({
            meshL <- react_file_sel()
            if(!is.null(meshL)) {
                pairL <- get_mesh_pairs(meshL)
                selectInput("rgl_view_select",
                            "Select mesh pair",
                            choices=names(pairL),
                            multiple=FALSE)
            } else {
                NULL
            }
        })
        output$rgl_mesh1_name <- renderUI({
            p(get_name_elem(input$rgl_view_select, 1))
        })
        output$rgl_mesh2_name <- renderUI({
            p(get_name_elem(input$rgl_view_select, 2))
        })
        output$rgl_dist1_name <- renderUI({
            p(get_name_elem(input$rgl_view_select, 1))
        })
        output$rgl_dist2_name <- renderUI({
            p(get_name_elem(input$rgl_view_select, 2))
        })
        output$table_agree_pairwise <- DT::renderDataTable({
            # input$apply_compare
            # isolate({
                d_agree_pairW <- react_mesh_agree()
                
                if(!is.null(d_agree_pairW)) {
                    cols_numeric <- unname(which(vapply(d_agree_pairW, is.numeric, logical(1))))
                    DT_out <- DT::datatable(d_agree_pairW)
                    DT::formatRound(DT_out, columns=cols_numeric, digits=3)
                } else {
                    NULL
                }
            # })
        })
        output$table_agree_aggr <- DT::renderDataTable({
            # input$apply_compare
            # isolate({
                d_agree_pairW <- react_mesh_agree()
                if(!is.null(d_agree_pairW)) {
                    d_agree_aggr <- get_mesh_agree_aggr(d_agree_pairW)                    
                    cols_numeric <- unname(which(vapply(d_agree_aggr, is.numeric, logical(1))))
                    DT_out <- DT::datatable(d_agree_aggr)
                    DT::formatRound(DT_out, columns=cols_numeric, digits=3)
                } else {
                    NULL
                }
            # })
        })
        output$rgl_mesh1 <- renderRglwidget({
            meshL <- react_file_sel()
            view_select <- input$rgl_view_select
            if(!is.null(view_select)) {
                mesh_sel <- get_name_elem(view_select, 1)
                mesh     <- meshL[[mesh_sel]]
                if(!is.null(mesh)) {
                    try(rgl.close())
                    wire3d(mesh$mesh)
                    rglwidget()
                } else {
                    NULL
                }
            } else {
                NULL
            }
        })
        output$rgl_mesh2 <- renderRglwidget({
            meshL <- react_file_sel()
            view_select <- input$rgl_view_select
            if(!is.null(view_select)) {
                mesh_sel <- get_name_elem(view_select, 2)
                mesh     <- meshL[[mesh_sel]]
                if(!is.null(mesh)) {
                    try(rgl.close())
                    wire3d(mesh$mesh)
                    rglwidget()
                } else {
                    NULL
                }
            } else {
                NULL
            }
        })
        output$rgl_mesh_union <- renderRglwidget({
            mesh_uiL    <- react_mesh_ui()
            view_select <- input$rgl_view_select
            if(!is.null(mesh_uiL) && !is.null(view_select)) {
                ui <- mesh_uiL[[view_select]]
                if(!is.null(ui) && !is.null(ui$union)) {
                    try(rgl.close())
                    wire3d(ui$union)
                    rglwidget()
                } else {
                    NULL
                }
            } else {
                NULL
            }
        })
        output$rgl_mesh_intersection <- renderRglwidget({
            mesh_uiL    <- react_mesh_ui()
            view_select <- input$rgl_view_select
            if(!is.null(mesh_uiL) && !is.null(view_select)) {
                ui <- mesh_uiL[[view_select]]
                if(!is.null(ui) && !is.null(ui$intersection)) {
                    try(rgl.close())
                    wire3d(ui$intersection)
                    rglwidget()
                } else {
                    NULL
                }
            } else {
                NULL
            }
        })
        output$rgl_mesh_dist1 <- renderRglwidget({
            metroL      <- react_mesh_metro()
            view_select <- input$rgl_view_select

            if(!is.null(metroL) && !is.null(view_select)) {
                metro <- metroL[[view_select]]
                if(!is.null(metro)) {
                    try(rgl.close())
                    shade3d(metro$mesh1)
                    rglwidget()
                } else {
                    NULL
                }
            } else {
                NULL
            }
        })
        output$rgl_mesh_dist2 <- renderRglwidget({
            metroL      <- react_mesh_metro()
            view_select <- input$rgl_view_select

            if(!is.null(metroL) && !is.null(view_select)) {
                metro <- metroL[[view_select]]
                if(!is.null(metro)) {
                    try(rgl.close())
                    shade3d(metro$mesh2)
                    rglwidget()
                } else {
                    NULL
                }
            } else {
                NULL
            }
        })
    }
)
