#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
## Global objects and functions for shiny app
#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------

options(rgl.useNULL=TRUE)
library(shiny)

## saved builtin data has class mesh3d, not cgalMesh
src_dir <- "www" # system.file("extdata", package="MeshAgreement")
f_name  <- paste0(src_dir, "/mesh3d_heart_obsL.rda")

cgalMesh_heart_obsL <- if(file.exists(f_name)) {
    mesh3d_heart_obsL <- readRDS(paste0(f_name))
    mesh3dL_to_cgalMeshL(mesh3d_heart_obsL)
} else {
    NULL
}

#####---------------------------------------------------------------------------
## popup
#####---------------------------------------------------------------------------

popup_info_impressum <- modalDialog(
    tagList(
        h3("Herausgeber"),
        p("Institut für Medizinische Biometrie, Epidemiologie und Informatik (IMBEI)", br(),
          "Universitätsmedizin der Johannes Gutenberg-Universität Mainz", br(),
          "Obere Zahlbacher Straße 69", br(),
          "55131 Mainz", br(),
          "Tel +49 6131 17-3252", br(),
          "Fax +49 6131 17-2968", br(), br(),
          "Die Universitätsmedizin Mainz ist eine Körperschaft des öffentlichen Rechts.", br(),
          "Vorsitzender des Aufsichtsrates:", br(),
          "Staatssekretär Dr. Denis Alt", br(),
          "Ministerium für Wissenschaft, Weiterbildung und Kultur", br(),
          "Mittlere Bleiche 61", br(),
          "55116  Mainz", br(), br(),
          "Umsatzsteuer-Identifikationsnummer: DE149065652"),
        h3("Autor / Redaktion"),
        p("Dieses Dashboard wurde von PD Dr. Daniel Wollschläger
          <wollschlaeger@uni-mainz.de> entwickelt."),
        h3("Haftung"),
        h4("Disclaimer"),
        p("Das IMBEI versucht die Richtigkeit und Aktualität der auf dieser
          Internetpräsenz bereitgestellten Informationen zu gewährleisten.
          Trotzdem können Fehler und Unklarheiten nicht vollständig
          ausgeschlossen werden. Das IMBEI übernimmt deshalb keine Gewähr für
          die Aktualität, Richtigkeit, Vollständigkeit oder Qualität der
          veröffentlichten Informationen. Für Schäden materieller oder
          immaterieller Art, die durch die Nutzung oder Nichtnutzung der
          dargebotenen Informationen unmittelbar oder mittelbar verursacht
          werden, haftet das IMBEI nicht, sofern ihm nicht vorsätzliches oder
          grob fahrlässiges Verschulden angelastet werden kann. Das IMBEI
          behält es sich vor, Teile des Internet-Angebots oder das gesamte
          Angebot ohne Vorankündigung zu verändern, zu ergänzen, zu löschen
          oder die Veröffentlichung zeitweise oder endgültig aus dem Internet
          zu entfernen."),
        h4("Erklärung zu Links auf fremde Webseiten"),
        p("Von den eigenen Inhalten sind Querverweise (Links) auf die von
          anderen Anbietern bereitgehaltenen Inhalte zu unterscheiden. Durch
          diese hält die Universitätsmedizin der Johannes Gutenberg-Universität
          Mainz insofern fremde Inhalte zur Nutzung bereit. Für diese fremden
          Inhalte ist die Universitätsmedizin der Johannes Gutenberg-Universität
          Mainz nur dann verantwortlich, wenn von ihnen (d.h. auch von einem
          rechtswidrigen bzw. strafbaren Inhalt) positive Kenntnis vorliegt
          und es technisch möglich und zumutbar ist, deren Nutzung zu verhindern."),
        p("Bei Links handelt es sich allerdings stets um dynamische Verweise,
          deren Inhalte sich stetig ändern können. Die Universitätsmedizin der
          Johannes Gutenberg-Universität Mainz hat bei der erstmaligen
          Verknüpfung zwar den fremden Inhalt daraufhin überprüft, ob durch
          ihn eine mögliche zivilrechtliche oder strafrechtliche
          Verantwortlichkeit ausgelöst wird. Der Inhaltsanbieter ist aber
          nicht dazu verpflichtet, die Inhalte, auf die er in seinem Angebot
          verweist, ständig auf Veränderungen zu überprüfen, die eine
          Verantwortlichkeit neu begründen könnten. Erst wenn er feststellt
          oder von anderen darauf hingewiesen wird, dass ein konkretes Angebot,
          zu dem er einen Link bereitgestellt hat, eine zivil- oder
          strafrechtliche Verantwortlichkeit auslöst, wird der Verweis auf
          dieses Angebot umgehend aufgehoben, soweit dies technisch möglich
          und zumutbar ist."),
        h3("Urheberrecht"),
        p("Copyright (c) 2022, IMBEI - Universitätsmedizin der Johannes
          Gutenberg-Universität Mainz.", br(),
          "Quelle und Copyright der verwendeten und hier wiedergegebenen Daten
              sind jeweils kenntlich gemacht. Insbesondere liegt das Copyright
              für Daten zu Sterbefällen und zur Bevölkerung beim Statistischen
              Landesamt Rheinland-Pfalz.", br(),
          "Icons: Font Awesome by Dave Gandy - http://fontawesome.io/", br(), br(),
          "Alle Rechte vorbehalten. Alle Inhalte der Internetpräsenz des IMBEI
          sind urheberrechtlich geschützt. Für die Vervielfältigung, Bearbeitung,
          Übersetzung, Einspeicherung, Verarbeitung und Wiedergabe von Inhalten
          in Datenbanken oder anderen elektronischen Medien und Systemen muss
          die Zustimmung des Urhebers eingeholt werden. Wir erlauben das
          Fotokopieren und Herunterladen unserer Internetseiten für private,
          wissenschaftliche und nicht kommerzielle Zwecke. Das IMBEI ahndet
          Verletzungen seiner Urheberrechte. Wir gestatten ausdrücklich Zitate
          unserer Dokumente und Internetseiten und freuen uns, wenn Sie auf
          unsere Seiten verlinken."),
        h3("Datenschutz"),
        p("Es gelten die Datenschutz-Bestimmungen der",
          tags$a(href="https://www.unimedizin-mainz.de/footer/datenschutz.html",
                 "Universitätsmedizin Mainz"),
          ".")
    ),
    title=NULL,
    footer=modalButton("Ok"),
    size="l",
    easyClose=TRUE,
    fade=FALSE
)
