BiocManager::install("OmnipathR")
BiocManager::install("dnet", force = TRUE)
BiocManager::install("gprofiler2")

library(OmnipathR)
library(tidyr)
library(dnet)
library(gprofiler2)

setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/")
#Tutorial on: https://r.omnipathdb.org/articles/omnipath_intro.html

get_interaction_resources()

## The interactions are stored into a data frame.
interactions <- import_omnipath_interactions(resources=c("SignaLink3","PhosphoSite","SIGNOR"))
platelet_proteins <- read.csv("../interactions/platelet_prots_uniprot.txt", col.names = "uniprotid")

UniprotID.top.10 <- sapply(strsplit(rownames(top.10.sign), ";"), "[[", 1)
UniprotID.top.10 <- unique(UniprotID.top.10)


platelet_network <- interactions[interactions$source %in% platelet_proteins$uniprotid & 
                                   interactions$target %in% platelet_proteins$uniprotid,]

total.interactions.top.10 <- interactions[interactions$source %in% UniprotID.top.10 |
                                            interactions$target %in% UniprotID.top.10,]
platelet.interactions.top.10 <- platelet_network[platelet_network$source %in% UniprotID.top.10 | 
                                      platelet_network$target %in% UniprotID.top.10,]
top.10_primary.platelet <- unique(c(unique(platelet.interactions.top.10$source), 
                           unique(platelet.interactions.top.10$target)))

top.10_primary.total <- unique(c(unique(total.interactions.top.10$source), 
                                    unique(total.interactions.top.10$target)))
PKA_regulators <- platelet_network[platelet_network$target_genesymbol == "PRKACA",]

active_kinases <- c("CDK2", "MAPK3")

find_paths <- function(x,y) {
  print_path_vs(all_shortest_paths(OPI_g,from = x,
                                   to = y)$res,OPI_g)
}

for (i in 1:length(active_kinases)) {
  for (j in 1:length(PKA_regulators$source_genesymbol)) {
    find_paths(active_kinases[i],PKA_regulators$source_genesymbol[j])
  }
}
find_paths("SRC", "PRKACA")
  
library(ggVennDiagram)

ggVennDiagram(list(platelet_proteins$uniprotid, UniprotID))

## We visualize the first interactions in the data frame.
print_interactions(head(interactions))

OPI_g <- interaction_graph(interactions = platelet_network)

##It is to note that the functions print_path\_es and print_path\_vs display
#very similar results, but the first one takes as an input an edge sequence 
#and the second one a node sequence.

## Find and print shortest paths on the directed network between proteins
## of interest:
print_path_es(shortest_paths(OPI_g,from = "CDK2",to = "PRKACA",
    output = 'epath')$epath[[1]],OPI_g)

## Find and print all shortest paths between proteins of interest:
print_path_vs(all_shortest_paths(OPI_g,from = "CDK2",
    to = "PRKACA")$res,OPI_g)






## We apply a clustering algorithm (Louvain) to group proteins in
## our network. We apply here Louvain which is fast but can only run
## on undirected graphs. Other clustering algorithms can deal with
## directed networks but with longer computational times,
## such as cluster_edge_betweenness. These cluster methods are directly
## available in the igraph package.
OPI_g_undirected <- as.undirected(OPI_g, mode=c("mutual"))
OPI_g_undirected <- simplify(OPI_g_undirected)
cl_results <- cluster_fast_greedy(OPI_g_undirected)
## We extract the cluster where a protein of interest is contained
cluster_id <- cl_results$membership[which(cl_results$names == "ERBB2")]
module_graph <- induced_subgraph(OPI_g_undirected,
                                 V(OPI_g)$name[which(cl_results$membership == cluster_id)])






