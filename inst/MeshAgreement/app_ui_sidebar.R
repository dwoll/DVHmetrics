dashboardSidebar(
    skin="light",
    fixed=TRUE,
    minified=TRUE,
    collapsed=FALSE,
    status="primary",
    sidebarMenu(
        id="sidebar_tabs",
        menuItem(
            "Home",
            tabName="tab_home",
            icon=icon("house", lib="font-awesome")
        ),
        menuItem(
            "Load data",
            tabName="tab_data",
            icon=icon("file-arrow-up", lib="font-awesome")
        ),
        menuItem(
            "Agreement",
            tabName="tab_agreement",
            icon=icon("code-compare", lib="font-awesome")
        ),
        menuItem(
            "View",
            tabName="tab_view",
            icon=icon("image", lib="font-awesome")
        ),
        menuItem(
            "About",
            icon=icon("info", lib="font-awesome"),
            tabName="tab_about"
        )
    ),
    tags$br(),
    tags$a(
        href="https://www.unimedizin-mainz.de/imbei/", 
        tags$img(src="img_logo_umm.png", 
                 title="Logo Universitätsmedizin Mainz", 
                 width="206",
                 height="55")
    )
)
