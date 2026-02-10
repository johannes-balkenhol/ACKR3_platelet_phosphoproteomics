BiocManager::install("OmnipathR")
BiocManager::install("dnet")
BiocManager::install("gprofiler2")
install.packages('igraph')
BiocManager::install("NetIndices")

library(OmnipathR)
library(tidyr)
library(dnet)
library(gprofiler2)
library(igraph)
library(biomaRt)
library(NetIndices)
library(ggraph)

get_interaction_resources()

#### The interactions are stored into a data frame.
interactions <-
    import_omnipath_interactions()

interactions <-
    import_omnipath_interactions(resources=c("Adhesome","SignaLink3","PhosphoSite",
    "SIGNOR","connectomeDB2020"))



#### filter interactions with Uniprot IDs present in the dataset

uniprot_id <- sapply(strsplit(rownames(top.all), ";"), "[[", 1)

filtered_interactions <- interactions[interactions$source %in% uniprot_id &
                                        interactions$target %in% uniprot_id,]

## Build the directed protein interaction network
network <- graph_from_data_frame(filtered_interactions, directed = TRUE)

## Visualize the network
dev.new()
plot(network)



ksh[["tree_row"]][["labels"]]

#### Find the shortest path between two proteins
## tranlate unirpto to gene symbol for visualization
protein_data <- unique(data.frame(c(interactions$source, interactions$target), c(interactions$source_genesymbol, interactions$target_genesymbol)))
colnames(protein_data) <- c("uniprot_id","symbol")


protein_data[grep(pattern = "ARRB", x = protein_data$symbol),]
protein_data[grep(pattern = "Q9Y2L6", x = protein_data$uniprot_id),]
top.all[grep(pattern = "PTK2", x = top.all$symbol),]

protein1 <- "P49407" # Example protein 1 (ARRB)
protein2 <- "P62993" # Example protein 2 (GRB2)
protein3 <- "P42345" # Example protein 3 (MTOR)
protein4 <- "Q03135" # Example protein 4 (CAV1)
protein5 <- "Q7Z460" # Example protein 4 (CLASP1)
protein6 <- "O75122" # Example protein 4 (CLASP2)
protein7 <- "Q15942" # Example protein 4 (ZYX)
protein8 <- "Q05209" # Example protein 4 (...)

# shortest.paths(graph, v=igraph.vs.all(graph), mode = "all")
# get.shortest.paths(graph, from, mode = "all")
# average.path.length(graph, directed=TRUE, unconnected=TRUE)

#https://igraph.org/r/doc/distances.html


shortest_path <- get.shortest.paths(network, from = protein1, to = protein8, mode="out")

print(shortest_path)

shortest_path_id <- shortest_path$vpath[[1]]


#### Visualize shortest graph
## Filter the network to only include nodes and edges in the shortest path
uniprot_ids <- names(as.list(shortest_path_id))
subnetwork <- subgraph.edges(network, eids = as.numeric(E(network)[.from(uniprot_ids[-length(uniprot_ids)]) & .to(uniprot_ids[-1]) ]))


V(subnetwork)$name <- protein_data$symbol[match(V(subnetwork)$name, protein_data$uniprot_id)]

## aethetics 
V(subnetwork)$color <- "#D3D3D3"
V(subnetwork)$label.color <- "black"
E(subnetwork)$color <-"#000000"
E(subnetwork)$width <-3

dev.new()
par(mar=c(.1,.1,.1,.1))
plot.igraph(subnetwork,
layout=layout.fruchterman.reingold,
xlim = c(-1, 1),
ylim = c(-1, 1),
vertex.label.cex=1.5,
edge.arrow.size=1.5)



### or
g <- ggraph(subnetwork, layout = "fr") + 
  geom_edge_link(aes(width = weight, color = "red")) +
  geom_node_point(aes(color = "red"), size = 5) +
  scale_color_manual(values = c("black", "red")) +
  theme_void() +
  theme(legend.position = "none")
g



#### Visualize full network

# Convert Uniprot IDs to gene symbols
#gene_symbols <- get_genes(shortest_path)

# Print the shortest path with gene symbols
#cat(paste(gene_symbols, collapse = " -> "))

# Highlight the shortest path in the network plot
V(network)$color <- "#D3D3D3"
#V(network)[gene_symbols]$color <- "red"
V(network)[shortest_path_id]$color <- "#0077ff"
#V(network)[shortest_path_id]$size <- -10
V(network)[shortest_path_id]$label.color <- "black"

E(network,path=shortest_path_id)$color <-"darkgreen"
E(network,path=shortest_path_id)$width <-2


#E(network)[from(gene_symbols[-length(gene_symbols)]) & to(gene_symbols[-1]) ]$color <- "red"
#E(network)[.from(shortest_path_id[-length(shortest_path_id)]) & .to(shortest_path_id[-1]) ]$color <- "red"
#plot(network)

dev.new()
par(mar=c(.1,.1,.1,.1))
plot.igraph(network,
layout=layout.fruchterman.reingold,
xlim = c(-1, 1),
ylim = c(-1, 1),
vertex.label.cex=.5,
edge.arrow.size=.5)




#### other network properties

# https://assemblingnetwork.wordpress.com/2013/06/10/network-basics-with-r-and-igraph-part-ii-of-iii/
dev.new()
par(mar=c(.1,.1,.1,.1))
plot.igraph(network,
layout=layout.fruchterman.reingold,
vertex.size=7,
vertex.label.cex=.5,
edge.arrow.size=.5)

network      # Tells me that it is an IGRAPH object with 100 nodes and 197 links,
# made with the Barabasi algorithm
V(network)   # gives the vertex sequence
E(network)   # gives the edge sequence (edge list)
 
# The "GenInd()" function requires an input of an adjacency matrix
network.adj <- get.adjacency(network,sparse=F)

network.properties <- GenInd(network.adj)

in.deg.network<-degree(network,v=V(network),mode="in")
out.deg.network<-degree(network,v=V(network),mode="out")
all.deg.network<-degree(network,v=V(network),mode="all")
 
# Degree distribution is the cumulative frequency of nodes with a given degree
# this, like degree() can be specified as "in", "out", or "all"
deg.distr<-degree.distribution(network,cumulative=T,mode="all")

power<-power.law.fit(all.deg.network)

# Then I can plot the degree distribution
plot(deg.distr,log="xy",
ylim=c(.01,10),
bg="black",pch=21,
xlab="Degree",
ylab="Cumulative Frequency")

# And the expected power law distribution
lines(1:200,10*(1:200)^((-power$alpha)+1))

diameter(network)
nodes.diameter<-get.diameter(network)


