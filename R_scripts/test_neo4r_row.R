
library(neo4r) ## from github
library(magrittr)
library(httr)
library(igraph)
library(dplyr)
library(purrr)
library(rlist)
library(igraph)
library(plotly)

con <- neo4j_api$new(url = "http://localhost:7474", 
                     db = "protTest", user = "neo4j", 
                     password = "**", isV4 = TRUE)
status_code(GET("http://localhost:7474"))



getGraph <- function(UP.id, direction, length, limit=10) {
  # direction <- "down"
  # length <- 3
  if (direction=="down") {
    direct.str <- sprintf("-[:INTERACTION*1..%d]->", length)
  } else if (direction=="up") {
    direct.str <- sprintf("<-[:INTERACTION*1..%d]-", length)
  } else {
    direct.str <- sprintf("<-[:INTERACTION*1..%d]->", length)
  }
  
  G <- sprintf(
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
    
    # G$nodes$properties <- lapply(G$nodes$properties, 
    #                              function(x) list.remove(x, c("altProtNames", "altGeneNames")))
    # G$nodes$properties[[6]]
    for (i in 1:length(G$nodes$properties)) {
      G$nodes$properties[[i]][['altGeneNames']] <- paste(unlist(G$nodes$properties[[i]][['altGeneNames']]),
                                                         collapse = "|")
      G$nodes$properties[[i]][['altProtNames']] <- paste(unlist(G$nodes$properties[[i]][['altProtNames']]),
                                                         collapse = "|")
    }
    # lapply(G$nodes$properties, 
    #        function(x) unlist(x[['altGeneNames']]))
    # x=G$nodes$properties[[6]]
    
    G$nodes <- G$nodes %>%
      unnest_nodes(what = "properties")
    
    G$relationships <- G$relationships %>%
      unnest_relationships() %>%
      select(startNode, endNode, type, everything())

    
    return(G)
  }
  
}

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

normalize <- function(x) {x / sqrt(sum(x^2))}

direction <- "both"
length <- 4
#name <- "Setd2"
#UP.id <- "E9Q5F9"
#name <- "PRKCD"
name <- "CTSB"
#UP.id <- "Q05655"
UP.id <- "P17612"
limit <- 100


G <- getGraph(UP.id, direction=direction, length=length, limit=limit)


Rows <- getSubstrates(UP.id)

rowDF <- dplyr::bind_cols(Rows)
colnames(rowDF) <- names(Rows)

rowDF$Prot1protNameAlt <- ifelse(rowDF$Prot1protNameAlt=="null", NA, rowDF$Prot1protNameAlt)
rowDF$Prot2protNameAlt <- ifelse(rowDF$Prot2protNameAlt=="null", NA, rowDF$Prot2protNameAlt)

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
  }
)

downloadSubs <- rowDF[,c(1:3, 5:9, 11:14)]
colnames(downloadSubs) <- c("Prot Gene Name", "Prot UniProt ID", "Prot Protein Name", "Prot taxid", "Prot organism",
                            "Sub Gene Name", "Sub UniProt ID", "Sub Protein Name", "Sub taxid", "Sub organism",
                            "Reaction type", "Reaction info")

downloadSubs$`Reaction info` <- apply(downloadSubs, 1, function(x) {

  l <- jsonlite::fromJSON(x["Reaction info"][[1]])
  
})
  
# rowDF$`Relationship details` <- lapply(rowDF$`Relationship details`, function(x) {
#   string = ""
#   l <- jsonlite::fromJSON(x[[1]])
#   for (i in 1:length(l)) {
#     string <- paste0(string, '<b>', i, ". ", "</b>")
#     l.i <- jsonlite::fromJSON(l[[i]])
#     if (!is.null(l.i$interactionID)) {
#       if (grepl("EBI-[0-9]+", l.i$interactionID)) {
#         l.i$interactionID = sprintf("<a href=https://www.ebi.ac.uk/intact/interaction/%s>%s</a>",
#                                     l.i$interactionID, l.i$interactionID)
#       } else if (grepl("CLE[0-9]+", l.i$interactionID)) {
#         link.id <- ex_between(l.i$`publicationID(s)`, "[", "]")[[1]]
#         publication <- ex_between(l.i$`publicationID(s)`, "<%", "[")[[1]]
#         publication <- gsub("%", "", publication)
#         l.i$`publicationID(s)` = sprintf("<a href=https://www.ebi.ac.uk/merops/cgi-bin/refs?id=%s>%s</a>",
#                                        link.id, publication)
#         l.i$interactionID
#       }
#     }
#     
#     for (j in paste(names(l.i), ":", l.i, "<br>")) {
#       string = paste0(string, j)
#     }
#   }
#   return(string)
# })

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



  

linksRow <- apply(rowDF, 1, function(row) {
  if (grepl('PhosphoSitePlus', row[["Relationship details"]])) {
    row[["Prot1UPID"]] <- sprintf("<a href=https://www.phosphosite.org/uniprotAccAction?id=%s>%s</a>", 
                                    row[["Prot1UPID"]], row[["Prot1UPID"]])
    row[["Prot2UPID"]] <- sprintf("<a href=https://www.phosphosite.org/uniprotAccAction?id=%s>%s</a>", 
                                    row[["Prot2UPID"]], row[["Prot2UPID"]])
  } ## CHEBI
    else if (grepl('CHEBI', row[["Prot2UPID"]])) {
    row[["Prot1UPID"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s'>%s</a>", 
                                  row[["Prot1UPID"]], row[["Prot1UPID"]])
    row[["Prot2UPID"]] <- sprintf("<a href='https://www.ebi.ac.uk/chebi/searchId.do;?chebiId=%s'>%s</a>", 
                                  row[["Prot2UPID"]], row[["Prot2UPID"]])
  } else {
    row[["Prot1UPID"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s'>%s</a>", 
                                    row[["Prot1UPID"]], row[["Prot1UPID"]])
    row[["Prot2UPID"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s'>%s</a>", 
                                    row[["Prot2UPID"]], row[["Prot2UPID"]])
  }
  c(row[["Prot1UPID"]], row[["Prot2UPID"]])
})
# dim(links)
# dim(t(links))
rowDF[,c("Prot1UPID", "Prot2UPID")] <- t(linksRow)
colnames(rowDF)
rowDF2 <- rowDF[,c(1:3, 5:9, 11:14)]



colnames(rowDF2) <- c("Prot Gene Name", "Prot UniProt ID", "Prot Protein Name", "Prot taxid", "Prot organism",
                   "Sub Gene Name", "Sub UniProt ID", "Sub Protein Name", "Sub taxid", "Sub organism",
                   "Reaction type", "Reaction info")

# df3$entries <- sapply(df3$entries, function(x) {
#   l = jsonlite::fromJSON(x[[1]])
#   string = ""
#   if (!is.null(l$interactionID)) {
#     if (grepl("EBI-[0-9]+", l$interactionID)) {
#       l$interactionID = sprintf("<a href=https://www.ebi.ac.uk/intact/interaction/%s>%s</a>",
#                                 l$interactionID, l$interactionID)
#     } else if (grepl("CLE[0-9]+", l$interactionID)) {
#       link.id <- ex_between(l$`publicationID(s)`, "[", "]")[[1]]
#       publication <- ex_between(l$`publicationID(s)`, "<%", "[")[[1]]
#       publication <- gsub("%", "", publication)
#       l$`publicationID(s)` = sprintf("<a href=https://www.ebi.ac.uk/merops/cgi-bin/refs?id=%s>%s</a>",
#                                      link.id, publication)
#     }
#   }
#   for (j in paste(names(l), ":", l, "\n"))
#     string = paste(string, j)
#   return(trimws(string, which="both"))
# })



##################################################################
# graph_object <- igraph::graph_from_data_frame(
#   d = G2$relationships, 
#   directed = TRUE, 
#   vertices = G2$nodes
# )

# plot(graph_object, edge.label=G2$relationships$name,
#      vertex.label=G2$nodes$name, edge.label.cex=0.8,
#      vertex.label.cex=0.7) 
##################################################################

direction <- "both"
length <- 3
#name <- "Setd2"
#UP.id <- "E9Q5F9"
#name <- "PRKCD"
name <- "CTSB"
UP.id <- "Q05655"
#UP.id <- "P08254"
limit <- 100


G <- getGraph(UP.id, direction=direction, length=length, limit=limit)


relation <- data.frame(G$relationships)
relation <- relation %>% dplyr::rename(reaction=name)
#G2$nodes <- G2$nodes[c("id", "organism", "taxid", "name", "uniprotID", "altProtNames", "proteinName", "label")] 



if ("proteinName" %in% colnames(G$nodes)) {
  G$nodes$proteinName <- apply(G$nodes, 1, function(x) {
    ifelse(is.na(x[['proteinName']]),
           strsplit(x[['altProtNames']], "|", fixed=T)[[1]][1],
           x[['proteinName']])
  })
} else {
  G$nodes$proteinName <- apply(G$nodes, 1, function(x) {
           strsplit(x[['altProtNames']], "|", fixed=T)[[1]][1]
           
  })
}
  

df1 <- left_join(relation, G$nodes, by = c("startNode"="id"))
df2 <- left_join(df1, G$nodes, by = c("endNode"="id"))
df3 <- df2[c("name.x", "uniprotID.x", "proteinName.x", "taxid.x", "organism.x",
             "name.y", "uniprotID.y", "proteinName.y", "taxid.y", "organism.y",
             "reaction", "entries")]


df3$entries <- apply(df3, 1, function(x) {
  l = jsonlite::fromJSON(x['entries'][[1]])
  string = ""
  if (!is.null(l$interactionID)) {
    if (grepl("EBI-[0-9]+", l$interactionID)) {
      l$interactionID = sprintf("<a href=https://www.ebi.ac.uk/intact/interaction/%s>%s</a>",
                                     l$interactionID, l$interactionID)
    } else if (grepl("CLE[0-9]+", l$interactionID)) {
      if (!is.null(l$`publicationID(s)`)) {
        link.id <- ex_between(l$`publicationID(s)`, "[", "]")[[1]]
        publication <- ex_between(l$`publicationID(s)`, "<%", "[")[[1]]
        publication <- gsub("%", "", publication)
        l$`publicationID(s)` <- sprintf("<a href=https://www.ebi.ac.uk/merops/cgi-bin/refs?id=%s>%s</a>",
                                        link.id, publication)
        l$interactionID <- sprintf("<a href=https://www.ebi.ac.uk/merops/cgi-bin/show_substrate?SpAcc=%s>%s</a>",
                                   x['uniprotID.y'],  l$interactionID)
      } else {
        l$interactionID <- sprintf("<a href=https://www.ebi.ac.uk/merops/cgi-bin/show_substrate?SpAcc=%s>%s</a>",
                                   x[[1]]['uniprotID.y'],  l$interactionID)
      }
    }
  }
  for (j in paste(names(l), ":", l, "<br>"))
    string = paste(string, j)
  return(trimws(string, which="both"))
})

# for (i in 1:length(df3$entries)) {
#   l = jsonlite::fromJSON(df3$entries[i][[1]])
#   string = ""
#   if (!is.null(l$interactionID)) {
#     if (grepl("EBI-[0-9]+", l$interactionID)) {
#       l$interactionID = sprintf("<a href=https://www.ebi.ac.uk/intact/interaction/%s>%s</a>",
#                                 l$interactionID, l$interactionID)
#     } else if (grepl("CLE[0-9]+", l$interactionID)) {
#       if (!is.null(l$`publicationID(s)`)) {
#         cat("yes")
#       } else {
#         "just link the CLE to merops"
#       }
#       link.id <- ex_between(l$`publicationID(s)`, "[", "]")[[1]]
#       publication <- ex_between(l$`publicationID(s)`, "<%", "[")[[1]]
#       publication <- gsub("%", "", publication)
#       l$`publicationID(s)` = sprintf("<a href=https://www.ebi.ac.uk/merops/cgi-bin/refs?id=%s>%s</a>",
#                                      link.id, publication)
#     }
#   }
#   for (j in paste(names(l), ":", l, "<br>"))
#     string = paste(string, j)
# }



Q7BU69
#grepl('PhosphoSitePlus', df3$entries)
links <- apply(df3, 1, function(row) {
  if (grepl('PhosphoSitePlus', row[["entries"]])) {
    row[["uniprotID.x"]] <- sprintf("<a href=https://www.phosphosite.org/uniprotAccAction?id=%s>%s</a>", 
                                    row[["uniprotID.x"]], row[["uniprotID.x"]])
    row[["uniprotID.y"]] <- sprintf("<a href=https://www.phosphosite.org/uniprotAccAction?id=%s>%s</a>", 
                                    row[["uniprotID.y"]], row[["uniprotID.y"]])
  } else {
    row[["uniprotID.x"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s'>%s</a>", 
                                    row[["uniprotID.x"]], row[["uniprotID.x"]])
    row[["uniprotID.y"]] <- sprintf("<a href='https://www.uniprot.org/uniprot/%s'>%s</a>", 
                                    row[["uniprotID.y"]], row[["uniprotID.y"]])
  }
  c(row[["uniprotID.x"]], row[["uniprotID.y"]])
})
# dim(links)
# dim(t(links))
df3[,c("uniprotID.x", "uniprotID.y")] <- t(links)

# df3$uniprotID.y <- sprintf("<a href='https://www.uniprot.org/uniprot/%s'>%s</a>", 
#                            df3$uniprotID.y, df3$uniprotID.y)

colnames(df3) <- c("Prot Gene Name", "Prot UniProt ID", "Prot Protein Name", "Prot taxid", "Prot organism",
                   "Sub Gene Name", "Sub UniProt ID", "Sub Protein Name", "Sub taxid", "Sub organism",
                   "Reaction type", "Reaction info")


#df3$`Reaction info` <- lapply(df3$`Reaction info`, function(x) gsub(',"', ', "', x))

# library(plotly)
# library(igraphdata)
# library(dplyr)
# library(jsonlite)
## this gets a nice dataframe for the edges
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
index.protein.searched <- which(G$nodes$uniprotID == UP.id)

L.g <- layout.circle(graph_object)
vs.g <- V(graph_object)
es.g <- get.edgelist(graph_object)

Nv.g <- length(vs.g) #number of nodes
Ne.g <- length(es.g[,1]) #number of edges

L.g <- layout.fruchterman.reingold(graph_object)
Xn.g <- L.g[,1]
Yn.g <- L.g[,2]

v.colors <- c("dodgerblue", "green", "yellow")
# for different interactions types 
e.colors <- c("orchid", "orange", "turquoise", "grey", "purple", "black", "darkblue", "")

#vertex_attr(graph_object, "UniProt ID", index = V(graph_object)) <- node.data.g$uniprotID
v.attrs <- vertex_attr(graph_object)
edge_attr(graph_object, "color", index = E(graph_object)) <- 
  e.colors[as.factor(edge_attr(graph_object)$name)]
e.attrs <- edge_attr(graph_object)
as.factor(e.attrs$name)


# Creates the nodes (plots the points)
# as.factor(data.frame(v.attrs$label))
#colors <- v.colors[as.factor(v.attrs$label)]
colors <- v.colors[forcats::fct_rev(as.factor(data.frame(v.attrs$label)))]
colors[index.protein.searched] <- "red"
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
                     size = 20),
                   textfont = list(color = '#000000', size = 16, layer="above"),
              )

#Create edges
edge_shapes.g <- list()
names(Xn.g) <- names(vs.g)
names(Yn.g) <- names(vs.g)


for(i in 1:Ne.g) {
  v0.g <- as.character(es.g[i,1])
  v1.g <- as.character(es.g[i,2])
  
  dir <- c(Xn.g[v1.g], Yn.g[v1.g]) - c(Xn.g[v0.g], Yn.g[v0.g])
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
    y1 = new.p2[2],
    name = e.attrs$name[i]
    # x0 = Xn.g[v0.g],
    # y0 = Yn.g[v0.g],
    # x1 = Xn.g[v1.g],
    # y1 = Yn.g[v1.g]
  )
  
  edge_shapes.g[[i]] <- edge_shape.g
}


axis.g <- list(title = "", showgrid = FALSE, 
             showticklabels = FALSE, zeroline = FALSE)
title <- ifelse(length==1 & direction=="down",
                sprintf("%s Substrates", name),
                sprintf("%s Paths", name))

p.g <- layout(
  network.g,
  title = list(text=title, 
               font=list(size=25)
               ),
  shapes = edge_shapes.g,
  xaxis = axis.g,
  yaxis = axis.g,
  showlegend=TRUE
) %>% add_trace(type="scatter")



arrow.x.start <- lapply(edge_shapes.g, function(x) x$x0) 
arrow.x.end <- lapply(edge_shapes.g, function(x) x$x1) 
arrow.y.start <- lapply(edge_shapes.g, function(x) x$y0) 
arrow.y.end <- lapply(edge_shapes.g, function(x) x$y1) 


ent <- lapply(e.attrs$entries, function(x) {
  string = ""
  l <- list()
  for (i in 1:length(x)) {
    string <- paste0(string, i, ". ")
    l[[i]] = jsonlite::fromJSON(x[i][[1]])
    for (j in paste(names(l[[i]]), ":", l[[i]], "\n")) {
      string = paste0(string, j)
    }
  }
  return(string)
})

x=e.attrs$entries[[1]]

p.g %>% plotly::add_trace(type = 'scatter') %>%
  
  plotly::add_annotations( x = ~arrow.x.end,
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
  
  plotly::add_annotations( x = ~arrow.x.end,
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

par(mar=c(1,1,1.4,1))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = levels(as.factor(e.attrs$name)), lty = 1, lwd = 3, 
       col = c(unique(e.colors))[1:nlevels(as.factor(e.attrs$name))], ncol=5,
       box.lty = 0, cex=1.2)
mtext("Reaction type", at=0.2, cex=2)
  
