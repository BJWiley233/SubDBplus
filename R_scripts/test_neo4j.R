library(neo4r)
library(magrittr)
library(httr)
library(igraph)
library(dplyr)
library(purrr)
library(rlist)

con <- neo4j_api$new(url = "http://localhost:7474", 
                     db = "protTest", user = "neo4j", 
                     password = "Swimgolf1212**", isV4 = TRUE)
status_code(GET("http://localhost:7474"))



getGraph3 <- function(UP.id, direction, length, limit=10) {
  
  if (direction=="down") {
    direct.str <- sprintf("-[:INTERACTION*0..%d]->", length)
  } else if (direction=="up") {
    direct.str <- sprintf("<-[:INTERACTION*0..%d]-", length)
  } else {
    direct.str <- sprintf("<-[:INTERACTION*0..%d]->", length)
  }
  
  G <- sprintf(
    "MATCH p = (n:Protein {uniprotID:'%s'})%s() 
      WITH *, relationships(p) AS rs
      UNWIND rs AS relType
      RETURN 
      	startNode(last(rs)).name AS Protein1, 
          relType['name'] AS `Interaction type`,
      	endNode(last(rs)).name AS Protein2, 
          rs[0] AS `Relationship details`, p AS Path
      LIMIT %d",  UP.id, direct.str, limit
  ) %>%
    call_neo4j(con, type="graph")
  
  if (length(G)==0) {
    return("No results in graph. Try another protein.")
  } else{
    #G2 <- G
    G$nodes$properties <- lapply(G$nodes$properties, 
                                 function(x) list.remove(x, c("altProtNames", "altGeneNames")))
    
    G$nodes <- G$nodes %>%
      unnest_nodes(what = "properties")
    
    G$relationships <- G$relationships %>%
      unnest_relationships() %>%
      select(startNode, endNode, type, everything())
    
    graph_object <- igraph::graph_from_data_frame(
      d = G$relationships, 
      directed = TRUE, 
      vertices = G$nodes
    )
    
    plot(graph_object, edge.label=G$relationships$name,
         vertex.label=G$nodes$name, edge.label.cex=0.8,
         vertex.label.cex=0.7) 
    return(G)
  }
  
}

UP.id="Q05655"
direction="down" 
length=1
limit=10
library(dplyr)
G <- getGraph3(UP.id, direction=direction, length=length, limit=10)
relation <- data.frame(G$relationships)
relation <- relation %>% rename(reaction=name)
G2 <- G$nodes[c("id", "organism", "taxid", "name", "uniprotID")] 

df1 <- left_join(relation, G2, by = c("startNode"="id"))
df2 <- left_join(df1, G2, by = c("endNode"="id"))
df3 <- df2[c("name.x", "uniprotID.x", "taxid.x", "organism.x",
             "name.y", "uniprotID.y", "taxid.y", "organism.y",
             "reaction", "entries")]
colnames(df3) <- c("Protein", "Protein UniProt ID", "Protein taxid", "Protein organism",
                   "Substrate", "Substrate UniProt ID", "Substrate taxid", "Substrate organism",
                   "Reaction type", "Reaction info")
df3$`Reaction info` <- lapply(df3$`Reaction info`, function(x) gsub(',"', ', "', x))
df3$`Reaction info`[[1]][1]
class(df3$`Reaction info`[[1]])
gsub(',"', ', "', df3$`Reaction info`[[1]])
