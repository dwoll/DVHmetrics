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
## HD_95:  Hausdorff distance - average of both 95th percentiles of directed surface distances
## JSC:    Jaccard similarity coefficient
## DSC:    Dice similarity coefficient
#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------

library(shiny)
library(bs4Dash)
library(shinyWidgets)
library(dplyr)
library(rgl)
library(MeshAgreement)
# library(Rvcg)
# library(SurfaceReconstruction)
# library(PolygonSoup)
# library(MeshesTools)
library(Boov)

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
            right=tagList(p(actionBttn(
                inputId="btn_footer_impressum",
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
        observeEvent(input$btn_home_go_about1, {
            updateTabItems(session, inputId="sidebar_tabs", selected="tab_about")
        })
        observeEvent(input$btn_home_go_about2, {
            updateTabItems(session, inputId="sidebar_tabs", selected="tab_about")
        })
        observeEvent(input$btn_footer_impressum, {
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
                                  input$file_sel$name)
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
                n_samples <- if(input$compare_limit) {
                    input$n_mesh_samples
                } else {
                    0
                }
                
                get_mesh_metro(meshL, n_samples=n_samples, chop=TRUE)
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
                n_samples <- if(input$compare_limit) {
                    input$n_mesh_samples
                } else {
                    0
                }
                
                mesh_pairL  <- get_mesh_pairs(meshL)
                agree_pairL <- Map(get_mesh_agree_pair,
                                   mesh_pairL,
                                   metro    =metroL,
                                   ui       =uiL,
                                   n_samples=n_samples)
                
                bind_rows(agree_pairL)
            } else {
                NULL
            }
        })
        output$ui_select_files <- renderUI({
            if(input$meshes_in == 2) {
                fileInput("file_sel", "Select files: (supported file formats: STL, PLY, OBJ, OFF)",
                          width="100%", multiple=TRUE)
            } else {
                NULL
            }
        })
        output$print_mesh_info <- renderUI({
            meshL <- react_file_sel()
            if(!is.null(meshL)) {
                Map(print_mesh_one, meshL, html=TRUE)
            } else {
                NULL
            }
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
                    DT::datatable(d_agree_pairW) %>%
                        DT::formatRound(columns=cols_numeric, digits=3)
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
                    d_agree_pairL <- get_mesh_agree_pair_long(d_agree_pairW)
                    d_agree_aggr  <- get_mesh_agree_aggr(d_agree_pairL)                    
                    cols_numeric  <- unname(which(vapply(d_agree_aggr, is.numeric, logical(1))))
                    DT::datatable(d_agree_aggr) %>%
                        DT::formatRound(columns=cols_numeric, digits=3)
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
