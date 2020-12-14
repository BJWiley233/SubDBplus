# Already inside server
output$pageStub <- renderUI(fluidPage(theme = "slate.min.css",
                                      tags$style(HTML("
                                                      .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter {
                                                      color: #a9a8ae;
                                                      }
                                                      #Error {
                                                        color: #da3910;
                                                      }
                                                      .action-button {
                                                        border-bottom-width: 1px !important;
                                                        border-bottom-style: dotted !important;
                                                        text-decoration:none !important;
                                                      }
                                                      .action-button:hover {
                                                        border-bottom-width: 0.5px !important;
                                                        border-bottom-style: solid !important;
                                                        background-color: #535160;
                                                        text-decoration:none;
                                                      }
                                                      ")
                                      ),
                                      
  titlePanel(title = "Search"),
  
  fluidRow(column(4, h4("Filters")), column(4), column(2, h6("\tFor substrates direction=down & length=1", 
                                                             style="color:#fdfd08"))),
           
  fluidRow(column(4, pickerInput("headProteinFamily", "Protein Family:",
                                choices = protein.families$headProteinFamily,
                                selected = protein.families$headProteinFamily,
                                multiple = TRUE,
                                options = list(`actions-box` = TRUE))),
   
          column(3, autocomplete_input("geneNamePreferred", "Gene Name: (type to choose IDs)", 
                                       options = unique(dat$geneNamePreferred[!is.na(dat$geneNamePreferred)]),
                                       # https://stackoverflow.com/questions/57798381/how-to-reduce-space-between-label-and-choices-in-selectinput
                                       create = T), div(style = "margin-top:-15px"),
                 # https://rdrr.io/cran/shinybrms/src/inst/shinybrms_app/app.R
                 p("Examples: ",
                 HTML(paste(actionLink("geneNamePreferredExampleCDC7", 
                                       "CDC7", style = "font-size:85%")), .noWS = "after"), ",",
                 HTML(paste0(actionLink("geneNamePreferredExamplePRKACA", 
                                        "PRKACA", style = "font-size:85%")), .noWS = "after"), ",",
                 HTML(paste0(actionLink("geneNamePreferredExampleCasp3", 
                                        "Casp3", style = "font-size:85%")), .noWS = "after"), ",",
                 HTML(paste0(actionLink("geneNamePreferredExampleBACE1", 
                                        "BACE1", style = "font-size:85%")), .noWS = "after"),
                 style = "font-size:85%; color:#fff")
                 ),
          column(1),
          column(2, selectInput("direction", "Path Direction", 
                                choices = c("down", "up", "both"), 
                                selected = "down")),
          column(2, textOutput("text"))),
  
  fluidRow(column(4, pickerInput("proteinOrganism", "Organism:",
                                choices = organism$organism,
                                selected = organism$organism,
                                multiple = TRUE,
                                options = list(`actions-box` = TRUE))),
           column(3, autocomplete_input("uniProtID", "UniProt ID:",
                                       options = dat$uniProtID,
                                       create=T), div(style = "margin-top:-15px"),
                  p("Examples: ",
                    HTML(paste(actionLink("UPIDExampleCDC7", 
                                          "O00311", style = "font-size:85%")), .noWS = "after"), ",",
                    HTML(paste0(actionLink("UPIDExamplePRKACA", 
                                           "P17612", style = "font-size:85%")), .noWS = "after"), ",",
                    HTML(paste0(actionLink("UPIDExampleCasp3", 
                                           "P70677", style = "font-size:85%")), .noWS = "after"), ",",
                    HTML(paste0(actionLink("UPIDExampleBACE1", 
                                           "P56817", style = "font-size:85%")), .noWS = "after"), ",",
                    HTML(paste0(actionLink("UPIDExampleACE2", 
                                           "Q9BYF1", style = "font-size:85%")), .noWS = "after"),
                    style = "font-size:85%; color:#fff"),
           
                    
           ),
           column(1),
           column(2, selectInput("length", "Path Length", 
                                choices = 1:5, selected = 1)),
           column(2, textOutput("text2"))),
  
  fluidRow(column(7), column(1, textOutput("upID")),
           column(3, sliderInput("limit", "Limit Relationship Results to:", min=1, max=200, value=10),
                  bsTooltip("limit", paste0("The limit of relationships to return. ",
                                            "Note since there can be multiple entries ",
                                            "per relationship, i.e. a kinase may phosphorylate ",
                                            "more than one site on a substrate, you may get more ",
                                            "results than your limit in the table."),
                            placement = "bottom", trigger = "hover",
                            options = list(width='200px'))
                  # bsTooltip("limit", "The limit of relationships to return.\n", 
                  #           placement = "right", trigger = "hover",
                  #           options = NULL)
                  )),
  
  fluidRow(column(7, DT::dataTableOutput("tb_chosen")), column(1),
           column(1, actionButton("getGraph", "Get Graph")),
           column(3, textOutput("Error")))
  )
)

#output$geneNamePreferredExamples <- renderText({ "Examples:" })
#output$comma <- renderText({ "," })

observe({
  update_autocomplete_input(session, "geneNamePreferred", "Gene Name: (type to choose IDs)",
                            options = unique(dat[dat$headProteinFamily %in% input$headProteinFamily &
                                          dat$organism %in% input$proteinOrganism,
                                          "geneNamePreferred"][!is.na(dat[dat$headProteinFamily %in% input$headProteinFamily &
                                                                      dat$organism %in% input$proteinOrganism,
                                                                      "geneNamePreferred"])])
                            )


  update_autocomplete_input(session, "uniProtID", "UniProt ID:",
                            options = dat[dat$headProteinFamily %in% input$headProteinFamily &
                                          dat$organism %in% input$proteinOrganism,
                                          "uniProtID"])
  
  # testing selections
  # output$text <- renderText({ nchar(input$geneNamePreferred) })
  # output$text2 <- renderText({ length(input$tb_chosen_cells_selected) })
  # output$upID <- renderText({ input$uniProtID })
  
  req(input$geneNamePreferred)
  ## Anything below here requires gene name to be entered for observing
  
  # if delete Gene Name
  if (nchar(input$geneNamePreferred) == 0) {
    update_autocomplete_input(session, "uniProtID", "UniProt ID:",
                              options = dat[dat$headProteinFamily %in% input$headProteinFamily &
                                            dat$organism %in% input$proteinOrganism, "uniProtID"],
                              value = "", create=T)
  } else {
    UP.ids <- dat[dat$headProteinFamily %in% input$headProteinFamily &
                  dat$organism %in% input$proteinOrganism &
                  dat$geneNamePreferred %in% input$geneNamePreferred,
                  "uniProtID"]
    if (length(UP.ids)==1) { UP.ids = list(UP.ids) }

    update_autocomplete_input(session, "uniProtID", "UniProt ID:",
                              options = UP.ids,
                              value = "", create=T)
  }
  
  ## if selecting from table
  if(length(input$tb_chosen_cells_selected)==0) {
    
  } 
  else {
    indices <- input$tb_chosen_cells_selected
    UP.ids <- dat[dat$headProteinFamily %in% input$headProteinFamily &
                  dat$organism %in% input$proteinOrganism &
                  dat$geneNamePreferred %in% input$geneNamePreferred,
                  "uniProtID"]
    if (length(UP.ids)==1) {
      value <- UP.ids
      UP.ids = list(UP.ids)
    } else {
      value <- UP.ids[indices[1]]
    }
    
    update_autocomplete_input(session, "uniProtID", "UniProt ID:",
                              options = UP.ids, create = F,
                              value = value)  
    
  }
  #output$upID <- renderText({ input$uniProtID })
  
  
})

## submit button click
observeEvent(input$getGraph, {
  submit.click()
})


submit.click <- reactive({
  ## direction and length with always have entry really
  if (input$uniProtID == "" || input$direction == "" || input$length < 1) {
    output$Error <- renderText({ "Error: Please enter UniProt ID\n or search Gene Name to select an ID." })
  } else {
    updateNavbarPage(session=session, "graph_and_data")
    fname = paste0("Graph", ".R") # remove leading "?", add ".R"
    cat(paste0("Session filename: ", fname, ".\n"))      # print the URL for this session
    source(fname, local=TRUE)
  }
  
})

## There is a better JavaScript way to handle actionLinks I think
observeEvent(input$geneNamePreferredExampleCDC7, {
  update_autocomplete_input(session, "geneNamePreferred", value = "CDC7")
})
observeEvent(input$geneNamePreferredExamplePRKACA, {
  update_autocomplete_input(session, "geneNamePreferred", value = "PRKACA")
})
observeEvent(input$geneNamePreferredExampleCasp3, {
  update_autocomplete_input(session, "geneNamePreferred", value = "Casp3")
})
observeEvent(input$geneNamePreferredExampleBACE1, {
  update_autocomplete_input(session, "geneNamePreferred", value = "BACE1")
})

observeEvent(input$UPIDExampleCDC7, {
  update_autocomplete_input(session, "uniProtID", value = "O00311")
})
observeEvent(input$UPIDExamplePRKACA, {
  update_autocomplete_input(session, "uniProtID", value = "P17612")
})
observeEvent(input$UPIDExampleCasp3, {
  update_autocomplete_input(session, "uniProtID", value = "P70677")
})
observeEvent(input$UPIDExampleBACE1, {
  update_autocomplete_input(session, "uniProtID", value = "P56817")
})
observeEvent(input$UPIDExampleACE2, {
  update_autocomplete_input(session, "uniProtID", value = "Q9BYF1")
})


# https://community.rstudio.com/t/change-the-color-of-column-headers-in-dt-table/77343
output$tb_chosen <- DT::renderDataTable(
                        datatable(style = "bootstrap", class = "compact",
                                  subset(dat, 
                                         dat$headProteinFamily %in% input$headProteinFamily &
                                         dat$geneNamePreferred %in% input$geneNamePreferred &
                                         dat$organism %in% input$proteinOrganism),
                                  selection=list(mode = "single", target = "cell"),
                                  options = list(
                                    initComplete = JS(
                                      "function(settings, json) {",
                                      "$(this.api().table().header()).css({'color': '#fff'});",
                                      "}")),
                                  colnames = c("Head Protein Family", "Organism", 
                                               "Gene Name", "UniProt ID"))
                        )
                
