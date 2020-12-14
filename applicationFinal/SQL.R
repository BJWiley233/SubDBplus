output$pageStub <- renderUI(
  fluidPage(theme = "slate.min.css",
            tags$style(
              HTML("
               .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter {
                 color: #a9a8ae;
               ")
            ),
            # Application title
            titlePanel("Search MySQL Database"),
            navbarPage("SQL", id="sqlNav",
                       tabPanel(title="Table Schema & SQL", value = "query",
                          fluidRow(column(8, DT::dataTableOutput("schema")),
                                   column(4, textAreaInput("sql", "SQL Query",
                                            placeholder = 
                                    sprintf(paste0(
                                            "/* Which kinase has most unique substrates */\n",
                                            "SELECT *\n",
                                            "FROM kinasePhosphoSitePlus\n",
                                            "WHERE uniProtIDKin = \n",
                                            "   (SELECT sel1.uniProtIDKin\n",
                                            "    FROM (SELECT DISTINCT uniProtIDKin, uniProtIDSub\n",
                                            " \t\t FROM kinasePhosphoSitePlus) as sel1\n",
                                            "    GROUP BY uniProtIDKin\n",
                                            "    ORDER BY COUNT(uniProtIDKin) DESC\n",
                                            "    LIMIT 1);")), width = "400px", height = "300px"),
                                   actionButton("sqlSubmit", "Submit")))
                          ),
                       tabPanel(title="Results", value="Results",
                                fluidRow(column(11, DT::dataTableOutput("results")),
                                         column(1, downloadButton("downloadData", "Download")))
                       )
            )
  )
  
)

getData <- reactive({
  req(input$sql)
  data <- loadData(input$sql)
  return(data)
})

observeEvent(input$sqlSubmit, {
  updateTabsetPanel(session, "sqlNav",
                    selected = "Results")
})

get.data.on.submit <- eventReactive(input$sqlSubmit, {
  data.sql <- getData()
  if (any(c('uniProtIDKin','geneNamePreferredKin','kinTaxid',
           'kinOrganism', 'subModSite', 'sitePlusMinus7AA') %in% colnames(data.sql))) {
    "<a href='https://www.phosphosite.org/homeAction' target='_blank'>PSP</a>"
    data.sql$source <- "<a href='https://www.phosphosite.org/homeAction' target='_blank'>PSP</a>"
  }
  
  ## if kinase substrate mod site with substrate ID link both mod site
  ## and substrate ID to PSP
  if (all(c('subModSite', 'uniProtIDSub') %in% colnames(data.sql))) {
    data.sql$subModSite <- sprintf("<a href='https://www.phosphosite.org/uniprotAccAction?id=%s' target='_blank'>%s</a>",
                                   data.sql$uniProtIDSub, data.sql$subModSite)
    data.sql$uniProtIDSub <- sprintf("<a href='https://www.phosphosite.org/uniprotAccAction?id=%s' target='_blank'>%s</a>",
                                     data.sql$uniProtIDSub, data.sql$uniProtIDSub)
    
  } 
  
  ## if kinase ID link to that kinase in PSP
  if ('uniProtIDKin' %in% colnames(data.sql)) {
    data.sql$uniProtIDKin <- sprintf("<a href='https://www.phosphosite.org/uniprotAccAction?id=%s' target='_blank'>%s</a>",
                                     data.sql$uniProtIDKin, data.sql$uniProtIDKin)
  }
  
  ## factor organism and headProteinFamily if selected
  if ('organism' %in% colnames(data.sql)) {
    data.sql$organism <- factor(data.sql$organism)
  }
  
  if ('headProteinFamily' %in% colnames(data.sql)) {
    data.sql$headProteinFamily <- factor(data.sql$headProteinFamily)
  }

  data.sql
 
})
 
  
  
output$results <- DT::renderDataTable(datatable(
  get.data.on.submit(),  style = "bootstrap", class = "compact",
  filter = "top",
  options = list(
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'color': '#fff'});",
      "}"),
    autoWidth = T,
    scrollX='800px', 
    scrollY = '500px',
    pageLength = 50), 
  callback = JS('table.page(3).draw(false);'),
  escape = F
  )
)

output$downloadData <- downloadHandler(
  filename = function() {
    paste0(file.prefix(), "_", gsub(" ", "_", date()), "_SQL.csv")
  },
  content = function(file) {
    write.csv(getData(), file, row.names = FALSE)
  }
)

output$schema <- DT::renderDataTable(datatable(
                      schema.mysql,  style = "bootstrap", class = "compact",
                      filter = "top",
                      options = list(
                        initComplete = JS(
                          "function(settings, json) {",
                          "$(this.api().table().header()).css({'color': '#fff'});",
                          "}"),
                        autoWidth = T,
                        columnDefs = list(
                          list(
                            targets = c(1,3,5,7,8),
                            render = JS(
                              "function(data, type, row, meta) {",
                              "return type === 'display' && data.length > 15 ?",
                              "'<span title=\"' + data + '\">' + data.substr(0, 15) + '...</span>' : data;",
                              "}"))),
                        scrollX='400px', 
                        scrollY = '500px',
                        pageLength = 50), 
                      callback = JS('table.page(3).draw(false);'),
                      colnames = c("Element Name", "Table Name", "Description",
                                   "Source","Data Type", "Node or Edge",
                                   "Neo4j Elem. Name", "Neo4j Data Type", "Neo4j Parent Data Type"),
                      escape = F
  )
)
