library(neo4r)
library(magrittr)

con <- neo4j_api$new(url = "http://localhost:7474", 
                     db = "protTest", user = "neo4j", 
                     password = "**", isV4 = TRUE)
con$ping()

library(httr)
status_code(GET("http://localhost:7474"))

sprintf(
"MATCH p = (n:Protein { name:'%s' })<-[:REGULATES*1..%d]->(b:Protein)
WITH *, relationships(p) AS rs
RETURN
  startNode(last(rs)).name AS Protein1,
  last(rs).direction AS Regulates,
  endNode(last(rs)).name AS Protein2,
  length(p) AS pathLength", 'C', 2) %>%
  call_neo4j(con)

G <- sprintf(
"MATCH p = (n:Protein { name:'%s' })<-[:REGULATES*1..%d]->(b:Protein)
WITH *, relationships(p) AS rs
RETURN
  p as PATH", 'C', 2) %>%
  call_neo4j(con, type="graph")

########################################################################
UP.id <- 'Q05655'
direction <- 'down'
length <- 1
limit <- 50
# Q05655
G <- sprintf(
"MATCH p = (n:Protein {uniprotID:'%s'})-[:INTERACTION*0..%d]->() 
WITH *, relationships(p) as rs
UNWIND rs as relType
RETURN 
	startNode(last(rs)).name as Protein1, 
    relType['name'] AS `Interaction type`,
	endNode(last(rs)).name as Protein2, 
    rs AS `Relationship details`, p as Path
 LIMIT %d",  UP.id, length, limit
) %>%
  call_neo4j(con, type="graph")
########################################################################

G
class(G$nodes$properties)
G$nodes
#sapply(G$nodes$properties, function(x) x$name)
library(ggraph)
library(igraph)
library(dplyr)
library(purrr)
library(rlist)
# Create a dataframe with col 1 being the ID, 
# And columns 2 being the names

G2 <- G
G2$nodes$properties <- lapply(G2$nodes$properties, function(x) list.remove(x, c("altProtNames", "altGeneNames")))
unnest_nodes(G2$nodes, what = "properties")

#what = "properties"
G2$nodes <- G2$nodes %>%
  unnest_nodes(what = "properties") #%>% 
  # We're extracting the first label of each node, but 
  # this column can also be removed if not needed
  #mutate(label = map_chr(label, 1))
head(G2$nodes)



G2$relationships <- G2$relationships %>%
  unnest_relationships() %>%
  select(startNode, endNode, type, everything())
head(G2$relationships)



graph_object <- igraph::graph_from_data_frame(
  d = G2$relationships, 
  directed = TRUE, 
  vertices = G2$nodes
)

labels <- lapply(G2$relationships$properties, function(x) x$name)

plot(graph_object, edge.label=G2$relationships$name,
     vertex.label=G2$nodes$name, edge.label.cex=0.8,
     vertex.label.cex=0.7) 
# TODO change color of edge depending on interactionType
# Add negative maybe for the inhibitor reactions in IntAct?
     
     #edge.color=ifelse(G$relationships$direction=="+", "red", "blue"),
     #edge.label.color=ifelse(G$relationships$direction=="+", "red", "blue"))











# library(ggraph)
# graph_object %>%
#   ggraph() + 
#   geom_node_label(aes(label = label)) +
#   geom_edge_link() + 
#   theme_graph()
# 
# class(graph_object)
# 
# "MATCH (n)-[r]->(m)
#     RETURN n,r,m" %>%
#   call_neo4j(con)
