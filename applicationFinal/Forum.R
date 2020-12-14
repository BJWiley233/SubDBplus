output$pageStub <- renderUI(
  fluidPage(theme = "slate.min.css",
            tags$style(
              HTML("
               .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter {
                 color: #a9a8ae;
               ")
            ),
            # Application title
            titlePanel(sprintf("FORUM (Coming Soon. Check my IntAct Table %s)",
                               emoji(search_emoji('smile'))[1])),
            fluidRow(column(12, DT::dataTableOutput("schema")))
  )
  
)


#intact[1:1000,]
output$schema <- DT::renderDataTable(datatable(
  intact[1:1000,],  style = "bootstrap", class = "compact",
  filter = "top",
  options = list(
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'color': '#fff'});",
      "}"),
    autoWidth = T,
    # takes to long to load
    # columnDefs = list(
    #   list(
    #     targets = c(6:9, 19),
    #     render = JS(
    #       "function(data, type, row, meta) {",
    #       "return type === 'display' && data.length > 30 ?",
    #       "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
    #       "}"))),
    scrollX='400px'), 
  callback = JS('table.page(3).draw(false);'),
  escape = F
)
)
