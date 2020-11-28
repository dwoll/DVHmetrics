bs4DashSidebar(
    skin="light",
    status="primary",
    brandColor="primary",
    # elevation=1,
    bs4SidebarMenu(
        id="sidebar",
        bs4SidebarMenuItem(
            "DVH data",
            tabName="tab_data",
            icon="file-upload"
        ),
        bs4SidebarMenuItem(
            "Metrics",
            tabName="tab_metrics",
            icon="table"
        ),
        bs4SidebarMenuItem(
            "Show DVH",
            tabName="tab_show_dvh",
            icon = "chart-line"
        ),
        bs4SidebarMenuItem(
            "Check constraints",
            tabName="tab_check_constraints",
            icon="table"
        ),
        bs4SidebarMenuItem(
            "Show constraints",
            tabName="tab_show_constraints",
            icon = "chart-line"
        ),
        bs4SidebarMenuItem(
            "BED/EQD2",
            tabName="tab_bed_eqd2",
            icon="table"
        ),
        bs4SidebarMenuItem(
            "About",
            tabName="tab_about",
            icon="lightbulb"
        )
    )
)
