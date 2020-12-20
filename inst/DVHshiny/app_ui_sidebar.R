dashboardSidebar(
    id="sidebar",
    skin="light",
    fixed=TRUE,
    minified=TRUE,
    collapsed=FALSE,
    status="primary",
    sidebarMenu(
        flat = FALSE,
        compact = FALSE,
        childIndent = TRUE,
        menuItem(
            "DVH data",
            tabName="tab_data",
            icon=icon("file-upload", lib="font-awesome")
        ),
        menuItem(
            "Metrics",
            tabName="tab_metrics",
            icon=icon("table", lib="font-awesome")
        ),
        menuItem(
            "Show DVH",
            tabName="tab_show_dvh",
            icon = icon("chart-line", lib="font-awesome")
        ),
        menuItem(
            "Check constraints",
            tabName="tab_check_constraints",
            icon=icon("table", lib="font-awesome")
        ),
        menuItem(
            "Show constraints",
            tabName="tab_show_constraints",
            icon = icon("chart-line", lib="font-awesome")
        ),
        menuItem(
            "BED/EQD2",
            tabName="tab_bed_eqd2",
            icon=icon("table", lib="font-awesome")
        ),
        menuItem(
            "About",
            tabName="tab_about",
            icon=icon("lightbulb", lib="font-awesome")
        )
    )
)
