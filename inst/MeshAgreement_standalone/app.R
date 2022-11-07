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
##
## TODO check input numeric ranges, types
#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------

library(shiny)
library(bs4Dash)
library(ggplot2)
library(plotly)
library(sortable)
library(rgl)
library(cgalMeshes)

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
                n_observers <- if(is.null(input$num_observers)) {
                    2L
                } else {
                    input$num_observers
                }
                
                meshL <- if(input$meshes_input_source == "builtin") {
                    ## use builtin data
                    ll <- data_heart_obsL
                    if(input$meshes_sel_mode == "indiv") {
                        ll[seq_len(min(c(length(ll), n_observers)))]
                    } else {
                        ll
                    }
                } else {
                    ## read data from file
                    ## number of file selection elements = number of observers
                    n_file_sel <- if(input$meshes_sel_mode == "all_pairwise") {
                        1L
                    } else {
                        n_observers
                    }
                    
                    ## for each file selection element -> read files
                    ll <- lapply(seq_len(n_file_sel), function(i) {
                        input_file_sel <- input[[sprintf("file_sel_%.2d", i)]]
                        if(!is.null(input_file_sel)) {
                            ## list of meshes
                            read_mesh_obs(input_file_sel$datapath,
                                          name=input_file_sel$name,
                                          fix_issues=input$read_mesh_fix_issues,
                                          reconstruct=input$read_mesh_reconstruct,
                                          spacing=input$read_mesh_reconstruct_pois_spacing)
                        } else {
                            NULL
                        }
                    })
                    
                    ll <- Filter(Negate(is.null), ll)
                    if(input$meshes_sel_mode == "all_pairwise") {
                        mesh_list_to_observer_list(unlist(ll, recursive=FALSE))
                    } else {
                        setNames(ll, sprintf("Observer_%.2d", seq_along(ll)))
                    }
                }
                
                meshL
            })
        })
        react_file_sel_sorted <- reactive({
            meshL <- react_file_sel()
            if(!is.null(meshL)) {
                lapply(seq_along(meshL), function(i) {
                    observer <- meshL[[i]]
                    ranklist <- input[[sprintf("ranklist_obs%.2d", i)]]
                    observer_sorted <- if(!is.null(ranklist)) {
                        observer[ranklist]
                    } else {
                        observer
                    }
                })
            } else {
                NULL
            }
        })
        react_mesh_ui <- reactive({
            meshL <- react_file_sel_sorted()
            if(!is.null(meshL)) {
                get_mesh_ui(meshL)
            } else {
                NULL
            }
        })
        react_mesh_metro <- reactive({
            meshL <- react_file_sel_sorted()
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
            meshL  <- react_file_sel_sorted()
            metroL <- react_mesh_metro()
            uiL    <- if(input$mesh_agree_do_ui) {
                react_mesh_ui()
            } else {
                list(NULL)
            }
            
            if(!is.null(meshL) && !is.null(metroL)) {
                mesh_pairL  <- get_mesh_pairs(meshL)
                agree_pairL <- Map(get_mesh_agree_pair,
                                   mesh_pairL,
                                   metro=metroL,
                                   ui   =uiL,
                                   do_ui=input$mesh_agree_do_ui)
                
                d <- do.call("rbind", agree_pairL)
                rownames(d) <- NULL
                d
            } else {
                NULL
            }
        })
        output$ui_select_comparisons <- renderUI({
            if(input$meshes_sel_mode == "indiv") {
                numericInput("num_observers", "Number of observers", value=2)
            } else {
                NULL
            }
        })
        output$ui_surface_recon_method <- renderUI({
            if(input$meshes_input_source == "file") {
                tagList(checkboxInput("read_mesh_fix_issues", "Try to fix mesh issues on import", TRUE),
                        radioButtons("read_mesh_reconstruct",
                                     "Surface reconstruction on import",
                                     choices=c("No", "AFS", "Poisson"),
                                     selected="No",
                                     inline=TRUE))
            } else {
                NULL
            }
        })
        output$ui_surface_recon_spacing <- renderUI({
            if(!is.null(input$meshes_input_source) &&
               !is.null(input$read_mesh_reconstruct) &&
               (input$meshes_input_source == "file") &&
               (input$read_mesh_reconstruct == "Poisson")) {
                numericInput("read_mesh_reconstruct_pois_spacing",
                             "Spacing parameter for Poisson reconstruction",
                             value=1)
            } else {
                NULL
            }
        })
        output$ui_select_files <- renderUI({
            if(input$meshes_input_source == "file") {
                if(input$meshes_sel_mode == "all_pairwise") {
                    n_observers <- 1L
                    sel_label   <- ""
                } else {
                    n_observers <- if(is.null(input$num_observers)) {
                        2L
                    } else {
                        input$num_observers
                    }
                    
                    sel_label <- paste0(" (Observer ", sprintf("%.2d", seq_len(n_observers)), ")")
                }
                
                file_selL <- lapply(seq_len(n_observers), function(i) {
                    finput_id    <- sprintf("file_sel_%.2d", i)
                    finput_label <- paste0("Select files", sel_label[i], ":")
                    column(width=12/n_observers,
                           fileInput(finput_id,
                                     finput_label,
                                     width="100%",
                                     multiple=TRUE))
                })
                
                ## weed out NULL components, convert the list to a tagList and return
                file_selL <- Filter(Negate(is.null), file_selL)
                fluidRow(do.call(tagList, file_selL))
            } else {
                NULL
            }
        })
        output$ui_ranklist_files <- renderUI({
            input$apply_file_sel
            isolate({
                if((input$meshes_sel_mode == "indiv")) {
                    n_observers <- if(is.null(input$num_observers)) {
                        2L
                    } else {
                        input$num_observers
                    }
                    
                    meshL <- react_file_sel()
                    if(!is.null(meshL)) {
                        n_obs_max <- min(c(length(meshL), n_observers))
                        ranklistL <- lapply(seq_len(n_obs_max), function(i) {
                            ranklist_ui_name <- sprintf("ui_ranklist_obs%.2d", i)
                            ranklist_inputid <- sprintf("ranklist_obs%.2d",    i)
                            meshL_obs <- meshL[[i]]
                            ranklist_labels <- vapply(meshL_obs, function(x) { x$name }, character(1))
                            column(width=12/n_obs_max,
                                   rank_list(sprintf("Files Observer %.2d", i),
                                             labels=ranklist_labels,
                                             input_id=ranklist_inputid))
                        })
                        
                        ## weed out NULL components, convert the list to a tagList and return
                        ranklistL <- Filter(Negate(is.null), ranklistL)
                        tagList(fluidRow(column(width=12,
                                                p("Drag-and-drop fifle names to define comparison sets.",
                                                  "All first elements are compared to each other between observers, and so on."))),
                                fluidRow(do.call(tagList, ranklistL)))
                    } else {
                        NULL
                    }
                } else {
                    NULL
                }
            })
        })
        output$ui_compare_table <- DT::renderDataTable({
            meshL <- react_file_sel_sorted()
            if(!is.null(meshL)) {
                pairL <- get_mesh_pairs(meshL, names_only=TRUE)
                DT::datatable(data.frame(Comparison=names(pairL)))
            } else {
                NULL
            }
        })
        output$print_mesh_info <- renderUI({
            meshL <- react_file_sel_sorted()
            if(!is.null(meshL)) {
                print_mesh_html(meshL)
            } else {
                NULL
            }
        })
        output$ui_mesh_agree_metro_options <- renderUI({
            tagList( numericInput("vcgMetro_nSamples",     "Number of samples", value=0),
                     numericInput("vcgMetro_nSamplesArea", "Number of samples per area (overrides nSamples)", value=0),
                    checkboxInput("vcgMetro_vertSamp",     "Vertex sampling",            value=TRUE),
                    checkboxInput("vcgMetro_edgeSamp",     "Edge sampling",              value=TRUE),
                    checkboxInput("vcgMetro_faceSamp",     "Face sampling",              value=TRUE),
                    checkboxInput("vcgMetro_unrefVert",    "Ignore unreferred vertices", value=FALSE),
                     radioButtons("vcgMetro_samplingTyp",  "Face sampling mode", choices=c("SS", "MC", "SD"), selected="SS", inline=TRUE),
                     radioButtons("vcgMetro_searchStruct", "Search structure",   choices=c("SGRID", "AABB", "OCTREE", "HGRID"), selected="SGRID", inline=TRUE),
                     numericInput("vcgMetro_from",         "Color mapping: minimum", value=0),
                     numericInput("vcgMetro_to",           "Color mapping: maximum", value=0)
            )
        })
        output$rgl_view_selection <- renderUI({
            meshL <- react_file_sel_sorted()
            if(!is.null(meshL)) {
                pairL <- get_mesh_pairs(meshL, names_only=TRUE)
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
                    DT_out <- DT::datatable(d_agree_pairW,
                                            extensions="Buttons",
                                            options=list(dom='Bfrtip',
                                                         buttons=c("csv", "excel")))
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
                    DT_out <- DT::datatable(d_agree_aggr,
                                            extensions="Buttons",
                                            options=list(dom='Bfrtip',
                                                         buttons=c("csv", "excel")))
                    DT::formatRound(DT_out, columns=cols_numeric, digits=3)
                } else {
                    NULL
                }
            # })
        })
        output$diag_agree_pairwise <- renderPlotly({
            d_agree_pairW <- react_mesh_agree()
            
            if(!is.null(d_agree_pairW)) {
                d_agree_pairL <- get_mesh_agree_long(d_agree_pairW)
                d_agree_pairL[["pair"]] <- paste(d_agree_pairL[["mesh1"]],
                                                 d_agree_pairL[["mesh2"]],
                                                 sep=" <-> ")
                p <- ggplot(d_agree_pairL, aes(x=pair, y=observed)) +
                    geom_point() +
                    facet_grid(metric ~ ., scales="free_y") +
                    xlab(NULL) +
                    ylab(NULL) +
                    theme_bw()
                
                ggplotly(p, height=600)
            } else {
                NULL
            }
            # })
        })
        output$diag_agree_aggr <- renderPlotly({
            d_agree_pairW <- react_mesh_agree()
            
            if(!is.null(d_agree_pairW)) {
                d_agree_aggrW <- get_mesh_agree_aggr(d_agree_pairW)
                d_agree_aggrL <- get_mesh_agree_aggr_long(d_agree_aggrW)
                # metric, statistic, observed
                p <- ggplot(d_agree_aggrL,
                            aes(x=statistic, y=observed,
                                group=group, color=group)) +
                    geom_point(position=position_dodge(w=0.2)) +
                    facet_grid(metric ~ ., scales="free_y") +
                    xlab(NULL) +
                    ylab(NULL) +
                    theme_bw()
                
                ggplotly(p, height=600)
            } else {
                NULL
            }
            # })
        })
        output$rgl_mesh1 <- renderRglwidget({
            meshL       <- react_file_sel_sorted()
            view_select <- input$rgl_view_select
            if(!is.null(meshL) && !is.null(view_select)) {
                pairL <- get_mesh_pairs(meshL, names_only=FALSE)
                mesh  <- pairL[[view_select]]$mesh1
                if(!is.null(mesh)) {
                    try(rgl.close())
                    wire3d(mesh$mesh$getMesh(rgl=TRUE))
                    rglwidget()
                } else {
                    NULL
                }
            } else {
                NULL
            }
        })
        output$rgl_mesh2 <- renderRglwidget({
            meshL       <- react_file_sel_sorted()
            view_select <- input$rgl_view_select
            if(!is.null(meshL) && !is.null(view_select)) {
                pairL <- get_mesh_pairs(meshL, names_only=FALSE)
                mesh  <- pairL[[view_select]]$mesh2
                if(!is.null(mesh)) {
                    try(rgl.close())
                    wire3d(mesh$mesh$getMesh(rgl=TRUE))
                    rglwidget()
                } else {
                    NULL
                }
            } else {
                NULL
            }
        })
        # output$rgl_mesh_union <- renderRglwidget({
        #     mesh_uiL    <- react_mesh_ui()
        #     view_select <- input$rgl_view_select
        #     if(!is.null(mesh_uiL) && !is.null(view_select)) {
        #         ui <- mesh_uiL[[view_select]]
        #         if(!is.null(ui) && !is.null(ui$union)) {
        #             try(rgl.close())
        #             wire3d(ui$union)
        #             rglwidget()
        #         } else {
        #             NULL
        #         }
        #     } else {
        #         NULL
        #     }
        # })
        # output$rgl_mesh_intersection <- renderRglwidget({
        #     mesh_uiL    <- react_mesh_ui()
        #     view_select <- input$rgl_view_select
        #     if(!is.null(mesh_uiL) && !is.null(view_select)) {
        #         ui <- mesh_uiL[[view_select]]
        #         if(!is.null(ui) && !is.null(ui$intersection)) {
        #             try(rgl.close())
        #             wire3d(ui$intersection)
        #             rglwidget()
        #         } else {
        #             NULL
        #         }
        #     } else {
        #         NULL
        #     }
        # })
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
