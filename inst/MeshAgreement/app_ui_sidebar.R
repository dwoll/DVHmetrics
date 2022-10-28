dashboardSidebar(
    skin="light",
    fixed=TRUE,
    minified=TRUE,
    collapsed=FALSE,
    status="primary",
    # brandColor="primary",
    # elevation=1,
    sidebarMenu(
        id="sidebar_tabs",
        menuItem(
            "Home",
            tabName="tab_home",
            # icon=icon("home", lib="font-awesome")
            icon=icon("house", lib="font-awesome")
        ),
        menuItem(
            "Load data",
            tabName="tab_data",
            # icon=icon("home", lib="font-awesome")
            icon=icon("file-arrow-up", lib="font-awesome")
        ),
        menuItem(
            "Agreement",
            tabName="tab_agreement",
            # icon=icon("map-marked-alt", lib="font-awesome")
            icon=icon("code-compare", lib="font-awesome")
        ),
        menuItem(
            "View",
            tabName="tab_view",
            icon=icon("image", lib="font-awesome")
        ),
        menuItem(
            "About",
            # icon=icon("info-circle", lib="font-awesome"),
            icon=icon("circle-info", lib="font-awesome"),
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
