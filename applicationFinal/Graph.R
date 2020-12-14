library(RMySQL)
library(DBI)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(dqshiny)
library(DT)
library(shinyjs)
library(plotly)
library(dplyr)
library(jsonlite)
library(qdapRegex)
library(emojifont)
library(shinyBS)
library(stringr)

output$pageStub <- renderUI(fluidPage(theme = "slate.min.css",
                                      tags$style(HTML("
                                                      .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter {
                                                        color: #a9a8ae;
                                                      }
                                                      #SwissModel a {
                                                        color: #5b33ff;
                                                      }
                                                      .has-feedback .form-control {
                                                        padding-right: 0px;
                                                      }
                                                      ")
                                      ),
  navbarPage("Graph and Data",
             tabPanel(title = "Graph",
                      fluidRow(column(2), plotOutput("legend", width="1200px", height="100px")),
                      fluidRow(plotlyOutput("plot", width="1200px", height="800px"))
             ),
             tabPanel(title="Graph Data Table",
                      fluidRow(column(12, 
                                        div(DT::dataTableOutput("data_table"),
                                            style = "font-size:95%; width:1200px")
                                      )
                        )
             ),
             tabPanel(title="Substrates",
                      fluidRow(column(11, 
                                      div(DT::dataTableOutput("substrates_table"),
                                          style = "font-size:95%; width:1200px")
                                      ),
                               column(1, downloadButton("downloadSubstrateData", "Download"))
                      )
             ),
             tabPanel(title="PDB Structures",
                      mainPanel(
                        fluidRow(column(7, DT::dataTableOutput("PDB_structures")),
                                 ## For 3D model by U of Pitt
                                 shinydashboard::box(title="3D Structure", width = 5, status="primary", solidHeader =TRUE, 
                                     uiOutput("structure_3d")),
                                 column(2, textOutput("text2"))),
                        fluidRow(column(8, htmlOutput("SwissModel"))),
                      )
             ),
             tabPanel(title="PDB Binding Sites and Drugs",
                      mainPanel(
                        fluidRow(column(8, DT::dataTableOutput("binding_drug"),
                                        style = "width:1200px"))
                      )
             )

    )
  )
)


observe({
  ## set name for graph title
  req(input$uniProtID)
  name <- dat[dat$uniProtID==input$uniProtID, "geneNamePreferred"]

  #' \link[app.R]{getGraph}
  G <- getGraph(input$uniProtID, input$direction,
                as.numeric(input$length), limit=as.numeric(input$limit))
  
  ## main graph app
  if ("neo" %in% class(G)) {
    
    ## Needs at least 1 protein to have a proteinName attribute needs fixing.
    # G$nodes$proteinName <- apply(G$nodes, 1, function(x){
    #   ifelse(is.na(x[['proteinName']]),
    #          strsplit(x[['altProtNames']], "|", fixed=T)[[1]],
    #          x[['proteinName']])
    # })
    if ("proteinName" %in% colnames(G$nodes)) {
      G$nodes$proteinName <- apply(G$nodes, 1, function(x) {
        ifelse(is.na(x[['proteinName']]),
               strsplit(x[['altProtNames']], "|", fixed=T)[[1]][1],
               x[['proteinName']])
      })
      ## might need to check if alt names is in there too but should really
      ## always be there with the Merge in Neo4j by Python class.
    } else {
      G$nodes$proteinName <- apply(G$nodes, 1, function(x) {
        strsplit(x[['altProtNames']], "|", fixed=T)[[1]][1]
        
      })
    }
    
    final <- G$relationships %>%
      group_by(id, startNode, endNode) %>% 
      summarise(
        startNode=startNode,
        endNode=endNode,
        type=type,
        id=id,
        entries=list(unique(entries)),
        name=name) %>%
      relocate(id, .after = type) %>%
      distinct()
    
    graph_object <- igraph::graph_from_data_frame(
      d = final, 
      directed = TRUE, 
      vertices = G$nodes
    )
    index.protein.searched <- which(G$nodes$uniprotID == input$uniProtID)

    L.g <- layout.circle(graph_object)
    vs.g<- V(graph_object)
    es.g <- get.edgelist(graph_object)
    
    Nv.g <- length(vs.g) #number of nodes
    Ne.g <- length(es.g[,1]) #number of edges
    
    L.g <- layout.fruchterman.reingold(graph_object)
    Xn.g <- L.g[,1]
    Yn.g <- L.g[,2]
    
    v.colors <- c("dodgerblue", "#0afb02", "#fcf51c")
    # for different interactions types 
    e.colors <- c("orchid", "orange", "#4444fb", "#5cf61d", "#6b6bae", 
                  "#0cf3fa", "#f20e42", "#cdafb6", "#8bd0f8", "#b40fb9", "#fdfbfd")
    
    
    v.attrs <- vertex_attr(graph_object)
    edge_attr(graph_object, "color", index = E(graph_object)) <- 
      e.colors[as.factor(edge_attr(graph_object)$name)]
    e.attrs <- edge_attr(graph_object)
    
    output$plot <- renderPlotly({
      
      ## set color of your protein to red, all others color of molecule type
      ## this factoring becomes issue with list vs. vectors with more than 1 factor
      #colors <- v.colors[as.factor(v.attrs$label)]
      ## change to list instead and Protein factor comes before Molecule using forcats::fct_rev
      colors <- v.colors[forcats::fct_rev(as.factor(data.frame(v.attrs$label)))]
      colors[index.protein.searched] <- "red"
      sizes <- rep(20, Nv.g)
      sizes[index.protein.searched] <- 30
      
      # Creates the nodes (plots the points)
      network.g <- plot_ly(x = ~Xn.g, y = ~Yn.g, #Node points
                           mode = "text+markers", 
                           text = vs.g$name, 
                           hoverinfo = "text",
                           hovertext = paste0("Gene Name: ", v.attrs$name, "\n",
                                              "Protein Name: ", v.attrs$proteinName, "\n",
                                              "UniProt ID: ", v.attrs$uniprotID, "\n",
                                              "Organism: ", v.attrs$organism, "\n",
                                              "TaxID: ", v.attrs$taxid),
                           marker = list(
                             color = colors,
                             size = sizes),
                           textfont = list(color = '#efeff5', size = 16, layer="above"),
      )
      
      #Create edges
      edge_shapes.g <- list()
      names(Xn.g) <- names(vs.g)
      names(Yn.g) <- names(vs.g)
      for(i in 1:Ne.g) {
        v0.g <- as.character(es.g[i,1])
        v1.g <- as.character(es.g[i,2])
        
        dir <- c(Xn.g[v1.g], Yn.g[v1.g]) - c(Xn.g[v0.g], Yn.g[v0.g])
        ## if self make small arrow
        if (all(dir == 0)) {
          new.p1 <- c(Xn.g[v1.g], Yn.g[v1.g])*(.9999)
          new.p2 <- c(Xn.g[v1.g], Yn.g[v1.g])*(1.0001)
        } else {
          new.p1 <- c(Xn.g[v0.g], Yn.g[v0.g]) + .2*normalize(dir)
          new.p2 <- c(Xn.g[v1.g], Yn.g[v1.g]) + -.1*normalize(dir)
        }
        
        edge_shape.g = list(
          type = "line",
          line = list(color = e.attrs$color[i], width = 2, layer="below"),
          opacity = 0.7,
          x0 = new.p1[1],
          y0 = new.p1[2],
          x1 = new.p2[1],
          y1 = new.p2[2]
        )
        
        edge_shapes.g[[i]] <- edge_shape.g
      }
      
      
      axis.g <- list(title = "", showgrid = FALSE, 
                     showticklabels = FALSE, zeroline = FALSE)
      title <- ifelse(input$length==1 & input$direction=="down",
                      sprintf("<b>%s Substrates", name),
                      sprintf("<b>%s Paths", name))
      
      p.g <- plotly::layout(
        network.g,
        title = list(text=title, 
                     font=list(size=30, style="italic", color="#c7c7df")
        ),
        shapes = edge_shapes.g,
        xaxis = axis.g,
        yaxis = axis.g,
        showlegend=FALSE,
        margin = list(l=50, r=50, b=100, t=100, pad=4),
        plot_bgcolor = "#19191f",
        paper_bgcolor = "#19191f"
      )
      
      arrow.x.start <- lapply(edge_shapes.g, function(x) x$x0) 
      arrow.x.end <- lapply(edge_shapes.g, function(x) x$x1) 
      arrow.y.start <- lapply(edge_shapes.g, function(x) x$y0) 
      arrow.y.end <- lapply(edge_shapes.g, function(x) x$y1) 
      
      ## edge properties
      ent <- lapply(e.attrs$entries, function(x) {
        string = ""
        t <- list()
        for (i in 1:length(x)) {
          string <- paste0(string, i, ". ")
          t[[i]] = jsonlite::fromJSON(x[i][[1]])
          for (j in paste(names(t[[i]]), ":", t[[i]], "\n")) {
            string = paste0(string, j)
          }
        }
        return(string)
      })
      
      p.g %>% add_trace(type = 'scatter') %>%
        
        add_annotations( x = ~arrow.x.end,
                         y = ~arrow.y.end,
                         xref = "x", yref = "y",
                         axref = "x", ayref = "y",
                         text = "",
                         hoverinfo = c(~arrow.x.end, ~arrow.y.end),
                         hovertext = paste(ent),
                         opacity = 0.7,
                         ax = ~arrow.x.end,
                         ay = ~arrow.y.end,
                         layer="below") %>%
        
        add_annotations( x = ~arrow.x.end,
                         y = ~arrow.y.end,
                         xref = "x", yref = "y",
                         axref = "x", ayref = "y",
                         text = "",
                         showarrow = T,
                         arrowcolor = ~e.attrs$color,
                         opacity = 0.7,
                         ax = ~arrow.x.start,
                         ay = ~arrow.y.start, 
                         layer="below")
      
    })
    
    output$legend <- renderPlot({
      par(mar=c(1,1,1.8,1))
      plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
      legend("topleft", legend = levels(as.factor(e.attrs$name)), lty = 1, lwd = 3, 
             col = c(unique(e.colors)), box.lty = 0, ncol = 5, cex=1.2, text.col = "#c7c7df")
      mtext("Reaction type", at=0.2, cex=2, col = "#c7c7df")
    }, bg = "#19191f")
    
    ## get table from graph
    relation <- data.frame(G$relationships)
    relation <- relation %>% rename(reaction=name)

    
    df1 <- left_join(relation, G$nodes, by = c("startNode"="id"))
    df2 <- left_join(df1, G$nodes, by = c("endNode"="id"))
    df3 <- df2[c("name.x", "uniprotID.x", "proteinName.x", "taxid.x", "organism.x",
                 "name.y", "uniprotID.y", "proteinName.y", "taxid.y", "organism.y",
                 "reaction", "entries")]
    
    ## list the entries
    df3$entries <- apply(df3, 1, function(x) {
      l = jsonlite::fromJSON(x['entries'][[1]])
      string = ""
      if (!is.null(l$interactionID)) {
        if (grepl("EBI-[0-9]+", l$interactionID)) {
          l$interactionID = sprintf("<a href='https://www.ebi.ac.uk/intact/interaction/%s' target='_blank'>%s</a>",
                                    l$interactionID, l$interactionID)
        } else if (grepl("CLE[0-9]+", l$interactionID)) {
          if (!is.null(l$`publicationID(s)`)) {
            link.id <- ex_between(l$`publicationID(s)`, "[", "]")[[1]]
            publication <- ex_between(l$`publicationID(s)`, "<%", "[")[[1]]
            publication <- gsub("%", "", publication)
            l$`publicationID(s)` <- sprintf("<a href='https://www.ebi.ac.uk/merops/cgi-bin/refs?id=%s' target='_blank'>%s</a>",
                                            link.id, publication)
            l$interactionID <- sprintf("<a href='https://www.ebi.ac.uk/merops/cgi-bin/show_substrate?SpAcc=%s' target='_blank'>%s</a>",
                                       x['uniprotID.y'],  l$interactionID)
          } else {
            l$interactionID <- sprintf("<a href='https://www.ebi.ac.uk/merops/cgi-bin/show_substrate?SpAcc=%s' target='_blank'>%s</a>",
                                       x['uniprotID.y'],  l$interactionID)
          }
        }
      }
      for (j in paste(names(l), ":", l, "<br>"))
        string = paste(string, j)
      return(trimws(string, which="both"))
    })
    
    ## According to PSP download agreement must make link to their site
    ## if displaying modification site information derived by PSP
    ##
    links <- apply(df3, 1, function(row) {
      if (grepl('PhosphoSitePlus', row[["entries"]])) {
        row[["uniprotID.x"]] <- sprintf("<a href='https://www.phosphosite.org/uniprotAccAction?id=%s' target='_blank'>%s</a>", 
                                        row[["uniprotID.x"]], row[["uniprotID.x"]])
        row[["uniprotID.y"]] <- sprintf("<a href='https://www.phosphosite.org/uniprotAccAction?id=%s' target='_blank'>%s</a>", 
                                        row[["uniprotID.y"]], row[["uniprotID.y"]])
      } ## CHEBI
      else if (grepl('CHEBI', row[["uniprotID.y"]])) {
        row[["uniprotID.x"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s' target='_blank'>%s</a>", 
                                      row[["uniprotID.x"]], row[["uniprotID.x"]])
        row[["uniprotID.y"]] <- sprintf("<a href='https://www.ebi.ac.uk/chebi/searchId.do;?chebiId=%s' target='_blank'>%s</a>", 
                                      row[["uniprotID.y"]], row[["uniprotID.y"]])
      } else {
        row[["uniprotID.x"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s' target='_blank'>%s</a>", 
                                        row[["uniprotID.x"]], row[["uniprotID.x"]])
        row[["uniprotID.y"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s' target='_blank'>%s</a>", 
                                        row[["uniprotID.y"]], row[["uniprotID.y"]])
      }
      c(row[["uniprotID.x"]], row[["uniprotID.y"]])
    })
    df3[,c("uniprotID.x", "uniprotID.y")] <- t(links)
      
    
    colnames(df3) <- c("Prot Gene Name", "Prot UniProt ID", "Prot Protein Name", "Prot taxid", "Prot organism",
                       "Sub Gene Name", "Sub UniProt ID", "Sub Protein Name", "Sub taxid", "Sub organism",
                       "Reaction type", "Reaction info")
    
    
    output$data_table <- DT::renderDataTable(
                            datatable(df3, style = "bootstrap", class = "compact",
                                      filter = "top",
                                      options = list(
                                        initComplete = JS(
                                          "function(settings, json) {",
                                          "$(this.api().table().header()).css({'color': '#fff'});",
                                          "}"),
                                        # https://github.com/rstudio/DT/issues/171
                                        autoWidth = T,
                                        width = "100%",
                                        scrollX=T,
                                        bSortClasses = TRUE,
                                        targets = 12,
                                        render = JS(
                                          "function(data, type, row, meta) {",
                                          "return type === 'display' && data.length > 60 ?",
                                          "'<span title=\"' + data + '\">' + data.substr(0, 60) + '...</span>' : data;",
                                          "}"),
                                        LengthMenu = c(5, 30, 50),
                                        columnDefs = list(
                                          list(className = 'dt-body-left'),
                                          list(width='325px', targets=12)),
                                        scrollY = '500px',
                                        pageLength = 50
                                      ),
                                      escape = F
                              )
    )
    ## no neo4j graph
  } else {
    output$plot <- renderPlotly({ empty_plot("No interaction data in neo4j for your protein!
There may be structure data on other tabs.
Check 'PDB Structures' or 'PDB Binding Sites and Drugs' tabs
or select a different UniProt ID.")
    })
  }
  
  ################################################################################################
  ## NEW!!! strictly a substrates table even if user selects some sort of pathway
  ## Ne4j query returns "row" type instead of "graph" type
  #' \link[app.R]{getSubstrates}
  Rows <- getSubstrates(input$uniProtID)
  if ("neo" %in% class(Rows) & length(Rows) > 0) {
    
    rowDF <- dplyr::bind_cols(Rows)
    colnames(rowDF) <- names(Rows)
    rowDF$Prot1protNameAlt <- ifelse(rowDF$Prot1protNameAlt=="null", NA, rowDF$Prot1protNameAlt)
    rowDF$Prot2protNameAlt <- ifelse(rowDF$Prot2protNameAlt=="null", NA, rowDF$Prot2protNameAlt)
    
    ## if Protein name is blank get first alternative
    rowDF$Prot1protName <- apply(rowDF, 1, function(x) {
      if(is.na(x['Prot1protName'])) {
        if (!is.na(x['Prot1protNameAlt'])) {
          return(jsonlite::fromJSON(x['Prot1protNameAlt'][[1]])[1])
        } else {
          return(NA)
        }
      } else {
        return(x['Prot1protName'])
      }
    }
    )
    
    ## if Sub name is blank get first alternative
    rowDF$Prot2protName <- apply(rowDF, 1, function(x) {
      if(is.na(x['Prot2protName'])) {
        if (!is.na(x['Prot2protNameAlt'])) {
          return(jsonlite::fromJSON(x['Prot2protNameAlt'][[1]])[1])
        } else {
          return(NA)
        }
      } else {
        return(x['Prot2protName'])
      }
    })
    
    ## for download without hyperlinks
    downloadSubs <- rowDF[,c(1:3, 5:9, 11:14)]
    colnames(downloadSubs) <- c("Prot Gene Name", "Prot UniProt ID", "Prot Protein Name", "Prot taxid", "Prot organism",
                          "Sub Gene Name", "Sub UniProt ID", "Sub Protein Name", "Sub taxid", "Sub organism",
                          "Reaction type", "Reaction info")
    
    
    ## Add links to Relationships
    rowDF$`Relationship details` <- apply(rowDF, 1, function(x) {
      string = ""
      l <- jsonlite::fromJSON(x["Relationship details"][[1]])
      for (i in 1:length(l)) {
        string <- paste0(string, '<b>', i, ". ", "</b>")
        l.i <- jsonlite::fromJSON(l[[i]])
        if (!is.null(l.i$interactionID)) {
          ## intact
          if (grepl("EBI-[0-9]+", l.i$interactionID)) {
            l.i$interactionID = sprintf("<a href='https://www.ebi.ac.uk/intact/interaction/%s' target='_blank'>%s</a>",
                                        l.i$interactionID, l.i$interactionID)
            ## merops
          } else if (grepl("CLE[0-9]+", l.i$interactionID)) {
            if (!is.null(l.i$`publicationID(s)`)) {
              link.id <- ex_between(l.i$`publicationID(s)`, "[", "]")[[1]]
              publication <- ex_between(l.i$`publicationID(s)`, "<%", "[")[[1]]
              publication <- gsub("%", "", publication)
              l.i$`publicationID(s)` <- sprintf("<a href='https://www.ebi.ac.uk/merops/cgi-bin/refs?id=%s' target='_blank'>%s</a>",
                                                link.id, publication)
              l.i$interactionID <- sprintf("<a href='https://www.ebi.ac.uk/merops/cgi-bin/show_substrate?SpAcc=%s' target='_blank'>%s</a>",
                                           x['Prot2UPID'],  l.i$interactionID)
            } else {
              l.i$interactionID <- sprintf("<a href='https://www.ebi.ac.uk/merops/cgi-bin/show_substrate?SpAcc=%s' target='_blank'>%s</a>",
                                           x['Prot2UPID'],  l.i$interactionID)
            }
          }
        }
        
        for (j in paste(names(l.i), ":", l.i, "<br>")) {
          string = paste0(string, j)
        }
      }
      return(string)
    })
    
    # Links to PSP, UniProt, and CHEBI
    linksRow <- apply(rowDF, 1, function(row) {
      if (grepl('PhosphoSitePlus', row[["Relationship details"]])) {
        row[["Prot1UPID"]] <- sprintf("<a href='https://www.phosphosite.org/uniprotAccAction?id=%s' target='_blank'>%s</a>",
                                      row[["Prot1UPID"]], row[["Prot1UPID"]])
        row[["Prot2UPID"]] <- sprintf("<a href='https://www.phosphosite.org/uniprotAccAction?id=%s' target='_blank'>%s</a>",
                                      row[["Prot2UPID"]], row[["Prot2UPID"]])
      } ## CHEBI
      else if (grepl('CHEBI', row[["Prot2UPID"]])) {
        row[["Prot1UPID"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s' target='_blank'>%s</a>", 
                                      row[["Prot1UPID"]], row[["Prot1UPID"]])
        row[["Prot2UPID"]] <- sprintf("<a href='https://www.ebi.ac.uk/chebi/searchId.do;?chebiId=%s' target='_blank'>%s</a>", 
                                      row[["Prot2UPID"]], row[["Prot2UPID"]])
      } else {
        row[["Prot1UPID"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s' target='_blank'>%s</a>",
                                      row[["Prot1UPID"]], row[["Prot1UPID"]])
        row[["Prot2UPID"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s' target='_blank'>%s</a>",
                                      row[["Prot2UPID"]], row[["Prot2UPID"]])
      }
      c(row[["Prot1UPID"]], row[["Prot2UPID"]])
    })
    
    rowDF[,c("Prot1UPID", "Prot2UPID")] <- t(linksRow)
    #rowDF2 <- rowDF[,c(1:8, 10:13)]
    rowDF2 <- rowDF[,c(1:3, 5:9, 11:14)]
    
    colnames(rowDF2) <- c("Prot Gene Name", "Prot UniProt ID", "Prot Protein Name", "Prot taxid", "Prot organism",
                         "Sub Gene Name", "Sub UniProt ID", "Sub Protein Name", "Sub taxid", "Sub organism",
                         "Reaction type", "Reaction info")
    
    output$substrates_table <- DT::renderDataTable(
                                datatable(rowDF2, style = "bootstrap", class = "compact",
                                          filter = "top",
                                          options = list(
                                            initComplete = JS(
                                              "function(settings, json) {",
                                              "$(this.api().table().header()).css({'color': '#fff'});",
                                              "}"),
                                            # https://github.com/rstudio/DT/issues/171
                                            autoWidth = T,
                                            width = "100%",
                                            scrollX=T,
                                            bSortClasses = TRUE,
                                            targets = 12,
                                            render = JS(
                                              "function(data, type, row, meta) {",
                                              "return type === 'display' && data.length > 60 ?",
                                              "'<span title=\"' + data + '\">' + data.substr(0, 60) + '...</span>' : data;",
                                              "}"),
                                            LengthMenu = c(5, 30, 50),
                                            columnDefs = list(
                                              list(className = 'dt-body-left'),
                                              list(width='325px', targets=12)),
                                            scrollY = '500px',
                                            pageLength = 50
                                          ),
                                          escape = F
                                )
    )
    
    output$downloadSubstrateData <- downloadHandler(
      filename = function() {
        paste0(file.prefix(), "_", gsub(" ", "_", date()), "_", input$uniProtID, "_SUBSTRATES.csv")
      },
      content = function(file) {
        
        write.csv(downloadSubs, file, row.names = FALSE)
      }
    )
    
  }
  ################################################################################################
  ## pdb structures tab
  pdb.data <- loadData(structures.query(input$uniProtID))
  pdb.data$pdbID <- sprintf("<a href='https://www.rcsb.org/structure/%s' target='_blank'>%s</a>", 
                            pdb.data$pdbID, pdb.data$pdbID)
  scrolly = "500px"
  if (nrow(pdb.data) == 0) {
    url <- a(input$uniProtID, href=sprintf("https://swissmodel.expasy.org/repository/uniprot/%s", 
                                           input$uniProtID), target='_blank')
    scrolly = "0px"
    output$SwissModel <- renderUI({  
      HTML(paste0("There are no structures for your UniProt protein.","<br>",
                     "Click link for Swiss-Model model of ", url))
    })
  }
 
  
  output$PDB_structures <- DT::renderDataTable(
                              datatable(pdb.data,  style = "bootstrap", class = "compact",
                                        filter = "top",
                                        selection=list(mode = "single", target = "cell"),
                                        options = list(
                                          initComplete = JS(
                                            "function(settings, json) {",
                                            "$(this.api().table().header()).css({'color': '#fff'});",
                                            "}"),
                                          scrollY = scrolly,
                                          pageLength = 25),
                                        escape = F
                              )
  )

  ################################################################################################
  ## binding sites and drugs
  #O00311
  pdb.drug.bind.data <- loadData(drugBankBinding.query(input$uniProtID))
  
  ## link to DrugBank
  pdb.drug.bind.data$drugBankID <- apply(pdb.drug.bind.data, 1, function(x) {
    if (is.na(x['drugBankID']) & !(is.na(x['ligandShort']))) {
      sprintf("<a href='https://go.drugbank.com/unearth/q?utf8=%%E2%%9C%%93&searcher=drugs&query=%s' target='_blank'>DB Search<a/>", 
              x['ligandShort'])
    } else if (is.na(x['drugBankID']) & (is.na(x['ligandShort']))) {
      sprintf("<a href='https://go.drugbank.com/unearth/q?utf8=%%E2%%9C%%93&searcher=drugs&query=' target='_blank'>DB Search<a/>", "")
    } else {
      sprintf("<a href='https://go.drugbank.com/drugs/%s' target='_blank'>%s<a/>", x['drugBankID'], x['drugBankID'])
    }
  })
  
  ## link to RCSB ligands
  pdb.drug.bind.data$ligandShort <- ifelse(is.na(pdb.drug.bind.data$ligandShort), NA,
                                           sprintf("<a href='https://www.rcsb.org/ligand/%s' target='_blank'>%s<a/>", 
                                                   pdb.drug.bind.data$ligandShort,
                                                   pdb.drug.bind.data$ligandShort))
  
  
  
  #pdb.drug.bind.data$ligandShort <- factor(pdb.drug.bind.data$ligandShort)
  
  pdb.drug.bind.data$pdbID <- sprintf("<a href='https://www.rcsb.org/structure/%s' target='_blank'>%s</a>",
                                      pdb.drug.bind.data$pdbID, pdb.drug.bind.data$pdbID)

  
  output$binding_drug <- DT::renderDataTable(
                              datatable(pdb.drug.bind.data, style = "bootstrap", class = "compact",
                                        filter = "top",
                                        options = list(
                                          initComplete = JS(
                                            "function(settings, json) {",
                                            "$(this.api().table().header()).css({'color': '#fff'});",
                                            "}"),
                                          #https://rstudio.github.io/DT/options.html
                                          autoWidth = T,
                                          width = "100%",
                                          scrollX=T,
                                          targets = 10,
                                          render = JS(
                                            "function(data, type, row, meta) {",
                                            "return type === 'display' && data.length > 10 ?",
                                            "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                                            "}")
                                          ,
                                          scrollY = '500px',
                                          pageLength = 50),
                                        colnames = c("UniProt Protein Chain", "PDB ID", "PDB Site ID", "Structure Residue #",
                                                     "UniProt Residue #", "Residue", "Residue Chain", "Ligand Residue #",
                                                     "Ligand Short", "Ligand Long", "Ligand Chain", "DrugBank ID"),
                                        escape = F
                              )
  )
  
  ################################################################################################
  ## 3D images from PITT javascript script, see works cited
  observe ({
    req(input$PDB_structures_cells_selected)
    
    if (length(input$PDB_structures_cells_selected)>0) {
      pdb <- ex_between(pdb.data[input$PDB_structures_cells_selected],">","</a")[[1]]
    } else {
      pdb=""
    }
      output$structure_3d <- renderUI({
        
        tabPanel("3D Structure",
                 tags$head(tags$script(src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js")),
                 tags$div(
                   style="height: 400px; width: 700px; position: relative;",
                   class='viewer_3Dmoljs',
                   'data-pdb'=pdb,
                   'data-backgroundcolor'='0xffffff',
                   'data-style'='cartoon'))
        
      })
    
  })

})

            