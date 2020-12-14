output$pageStub <- renderUI(
  fluidPage(theme = "slate.min.css",
            br(),
            #https://community.rstudio.com/t/how-to-embed-videos-in-shiny/38937/6
            HTML('<iframe width="760" height="445" src="https://www.youtube.com/embed/8ew1xJZf7tY" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
  )
)
