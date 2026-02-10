BiocManager::install("OmnipathR")
BiocManager::install("dnet")
BiocManager::install("gprofiler2")

library(OmnipathR)
library(tidyr)
library(dnet)
library(gprofiler2)

get_interaction_resources()

## The interactions are stored into a data frame.
interactions <-
    import_omnipath_interactions()

interactions <-
    import_omnipath_interactions(resources=c("SignaLink3","PhosphoSite",
    "SIGNOR"))

## We visualize the first interactions in the data frame.
print_interactions(head(interactions))

OPI_g <- interaction_graph(interactions = interactions)

## Find and print shortest paths on the directed network between proteins
## of interest:
print_path_es(shortest_paths(OPI_g,from = "TYRO3",to = "STAT3",
    output = 'epath')$epath[[1]],OPI_g)

    ## Find and print all shortest paths between proteins of interest:
print_path_vs(all_shortest_paths(OPI_g,from = "DYRK2",
    to = "MAPKAPK2")$res,OPI_g)


    OPI_g_undirected <- as.undirected(OPI_g, mode=c("mutual"))
OPI_g_undirected <- simplify(OPI_g_undirected)
cl_results <- cluster_fast_greedy(OPI_g_undirected)
## We extract the cluster where a protein of interest is contained
cluster_id <- cl_results$membership[which(cl_results$names == "ERBB2")]
module_graph <- induced_subgraph(OPI_g_undirected,
    V(OPI_g)$name[which(cl_results$membership == cluster_id)])



    ## We query and store the interactions into a dataframe
interactions <-
    import_pathwayextra_interactions(resources=c("BioGRID","STRING"),
    organism = 10090)



## We query and store the interactions into a dataframe
interactions <-
    import_kinaseextra_interactions(resources=c("PhosphoPoint",
    "PhosphoSite"), organism = 10116)

## We select the interactions in which Dpysl2 gene is a target
interactions_TargetDpysl2 <- dplyr::filter(interactions,
    target_genesymbol == "Dpysl2")

## We print these interactions:
print_interactions(interactions_TargetDpysl2)
