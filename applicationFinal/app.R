library(RMySQL)
library(DBI)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(dqshiny)
library(DT)
library(data.table)
library(shinyjs)
library(plotly)
library(dplyr)
library(jsonlite)
library(qdapRegex)
library(emojifont)
library(shinyBS)
library(stringr)

loadData <- function(query) {
  db <- RMySQL::dbConnect(RMySQL::MySQL(),
                          db = "protTest",
                          username = "root",
                          password = "**",
                          host = "127.0.0.1")
  dat <- dbGetQuery(db, query)
  dbDisconnect(db)
  dat
}

protein.families.query <- sprintf("SELECT DISTINCT headProteinFamily 
                                   FROM proteinsUniprot")
protein.families <- loadData(protein.families.query)

organisms.query <- sprintf("SELECT DISTINCT organism 
                            FROM proteinsUniprot")
organism <- loadData(organisms.query)

dat.query <- sprintf("SELECT headProteinFamily, organism, geneNamePreferred, uniProtID
                      FROM proteinsUniprot")
dat <- loadData(dat.query)

structures.query <- function(id) {
  sprintf("SELECT * FROM uniprotPdbJoin
           WHERE uniProtID = '%s'", id)
}

## for the Forum.R placeholder
intact <- loadData("SELECT * FROM intact")

cols <- c("organismA", "organismB", "detectionMethod", 
          "interactionType", "isNegative")
intact[cols] <- lapply(intact[cols], factor)


drugBankBinding.query <- function(id) {
  ## if drugBankID is blank
  ## https://go.drugbank.com/unearth/q?utf8=%E2%9C%93&searcher=drugs&query=%s, ligandShort
  sprintf("SELECT DISTINCT j.`chain` as 'UniProt Protein Chain', pb.*, sel.drugBankID
           FROM uniprotPdbJoin j
           INNER JOIN 
          	 pdbBindingSites pb
          	 ON j.pdbID = pb.pdbID
           LEFT JOIN 
          	 (SELECT DISTINCT d.drugShort, d.drugBankID FROM pdbDrugBank d) sel
              ON pb.ligandShort = sel.drugShort
           WHERE j.uniProtID = '%s'", id)
}

## Data Dictionary
schema.mysql <- read.table("~/JHU_Fall_2020/Biological_DBs/Project/data_dict.csv",
                           sep = ",", header = T, quote = "\"")
#https://stackoverflow.com/questions/33180058/coerce-multiple-columns-to-factors-at-once
cols <- c("MySQL.table.name", "Data.type",
          "Neo4j.node.or.edge..or.property.thereof.", "Neo4j.data.type")
schema.mysql[cols] <- lapply(schema.mysql[cols], factor)

# levels(schema.mysql$Source) <- unlist(lapply(levels(schema.mysql$Source), function (x) {
#   rvest::html_text(xml2::read_html(x))
# }))

## for shortening arrows in graph
normalize <- function(x) {x / sqrt(sum(x^2))}

## for download sql query results to file
#https://stackoverflow.com/questions/42734547/generating-random-strings
file.prefix <- function() {
  a <- do.call(paste0, replicate(7, sample(LETTERS, 1, TRUE), FALSE))
  b <- do.call(paste0, replicate(7, sample(LETTERS, 1, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(999999, 1, TRUE)), b)
}



empty_plot <- function(title = NULL){
  p <- plotly_empty(type = "scatter", mode = "markers") %>%
    plotly::config(
      displayModeBar = FALSE
    ) %>%
    layout(
      title = list(
        text = title,
        yref = "paper",
        y = 1.,
        x=0.02,
        font = list(color="red"),
        xanchor = "left"
      )
    )
  return(p)
} 

############################################################################
## Neo4j stuff
############################################################################
library(neo4r)
library(httr)
library(igraph)
library(rlist)

con <- neo4j_api$new(url = "http://localhost:7474", 
                     db = "protTest", user = "neo4j", 
                     password = "Swimgolf1212**", isV4 = TRUE)
status_code(GET("http://localhost:7474"))


## Function for Graph
## Used to get Graph Table as well
getGraph <- function(UP.id, direction, length, limit=10) {
  
  if (direction=="down") {
    direct.str <- sprintf("-[:INTERACTION*1..%d]->", length)
  } else if (direction=="up") {
    direct.str <- sprintf("<-[:INTERACTION*1..%d]-", length)
  } else {
    direct.str <- sprintf("<-[:INTERACTION*1..%d]->", length)
  }
  
  G <- sprintf(
    # "MATCH p = (n:Protein {uniprotID:'%s'})%s() 
    #   WITH *, relationships(p) AS rs
    #   UNWIND rs AS relType
    #   RETURN 
    #   	startNode(last(rs)).name AS Protein1, 
    #       relType['name'] AS `Interaction type`,
    #   	endNode(last(rs)).name AS Protein2, 
    #       rs[0] AS `Relationship details`, p AS Path
    #   LIMIT %d",  UP.id, direct.str, limit
    
    ## Above ne04j query is incorrect, this will
    ## now limit return of relationships correctly at least 
    ## instead of trying to limit entries or nodes which doesn't work
    "MATCH p = (n:Protein {uniprotID:'%s'})%s() 
      WITH *, relationships(p) AS rs
      RETURN 
        p
      LIMIT %d",  UP.id, direct.str, limit
  ) %>%
    call_neo4j(con, type="graph")
  
  if (length(G)==0) {
    return("No results in graph. Try another protein.")
  } else {
 
    ## tibble needs to be same length
    for (i in 1:length(G$nodes$properties)) {
      G$nodes$properties[[i]][['altGeneNames']] <- paste(unlist(G$nodes$properties[[i]][['altGeneNames']]),
                                                         collapse = "|")
      G$nodes$properties[[i]][['altProtNames']] <- paste(unlist(G$nodes$properties[[i]][['altProtNames']]),
                                                         collapse = "|")
    }

    G$nodes <- G$nodes %>%
      unnest_nodes(what = "properties")
    
    G$relationships <- G$relationships %>%
      unnest_relationships() %>%
      select(startNode, endNode, type, everything())
    
    return(G)
  }
}

## Function to get all substrates from selected UniProt ID
## Main SIMPLE goal from database
getSubstrates <- function(UP.id) {
  Rows <- sprintf(
    "MATCH p = (n:Protein {uniprotID:'%s'})-[:INTERACTION]->() 
     WITH *, relationships(p) as rs
     UNWIND rs as relationship
     RETURN 
    	  startNode(relationship).name as Protein1, 
        startNode(relationship).uniprotID as Prot1UPID,
        startNode(relationship).proteinName as Prot1protName,
        apoc.convert.toJson(startNode(relationship).altProtNames) as Prot1protNameAlt,
        startNode(relationship).taxid as Prot1tax,
        startNode(relationship).organism as Prot1org,
        endNode(relationship).name as Protein2, 
        endNode(relationship).uniprotID as Prot2UPID,
        endNode(relationship).proteinName as Prot2protName,
        apoc.convert.toJson(endNode(relationship).altProtNames) as Prot2protNameAlt,
        endNode(relationship).taxid as Prot2tax,
        endNode(relationship).organism as Prot2org,
        relationship['name'] as `Interaction type`,
        apoc.convert.toJson(relationship['entries']) AS `Relationship details`
    ",  UP.id
  ) %>%
    call_neo4j(con, type="row", output = "r")
  
  Rows
}



# https://github.com/open-meta/uiStub/blob/master/app/old-faithful.R
ui <- function(request) { uiOutput("uiStub") }

server <- function(input, output, session) {
  cat("Session started.\n")                               # this prints when a session starts
  onSessionEnded(function() {cat("Session ended.\n\n")})  # this prints when a session ends
  
  output$uiStub <- renderUI(
    fluidPage(
      theme = shinytheme("slate"),
      fluidRow(column(12, h1("Protein Substrates"))),
      fluidRow(column(12,
                      HTML("<h3><a href='?RoadMap'>RoadMap</a> | ",
                           "<a href='?Search'>Search</a> |",
                           "<a href='?SQL'>SQL</a> |",
                           "<a href='?Forum'>Forum</a> |",
                           "<a href='?Demonstration'>Demonstration</a> |",
                           "<a href='?page3'>Nothing</a>",
                           "</h3>"))
      ),
      uiOutput("pageStub")
    )
  )
  
  validFiles = c("RoadMap.R", "Search.R", "Graph.R", "SQL.R",
                 "Forum.R", "Demonstration.R")
  
  fname = isolate(session$clientData$url_search)
  if (nchar(fname)==0) { fname = "?RoadMap"}
  fname = paste0(substr(fname, 2, nchar(fname)), ".R") # remove leading "?", add ".R"
  
  cat(paste0("Session filename: ", fname, ".\n"))      # print the URL for this session
      
  if(!fname %in% validFiles){                          # is that one of our files?
    output$pageStub <- renderUI(tagList(              # 404 if no file with that name
      fluidRow(
        column(5,
               HTML("<h2>404 Not Found Error:</h2><p>That URL doesn't exist. Use the",
                    "menu above to navigate to the page you were looking for.</p>")
        )
      )
    ))
    return()    # to prevent a "file not found" error on the next line after a 404 error
  }
  source(fname, local=TRUE)                            # load and run server code for this page
}

# Run the application
shinyApp(ui = ui, server = server, enableBookmarking = "url")

