###################################################################################################################################################
###################################################################################################################################################
################# Ideas
### 1. the plan is to extract the central cascade and combine it in string with the significant differential abundant psites 
### repeat this for all timepints (problem: other timepoints have to many proteins)

### 2. the plan is to use the central cascade up to degree 3 neighbor and map the differential abundant sites on the network 
### repeat this process for all timepoints (problem: not all psites are well represented in ths network and entrez gene to uniprot conversion creates multiples)

### 3. the plan is to extract the central cascade and combine it in string with the significant differential abundant psites 
### that are of interest from other analysis (e.g. kinase target sites ) 
### repeat this for all timepints (problem: selective representation)

### 4. use string only to project only selected psites for each timepoint (prefred for small subnetworks)





################# export undirected network
#### get the uniprot ID

###### sif file
SELECT DISTINCT UPPER(gene_name1) AS gene_name1, inttype, UPPER(gene_name2) AS gene_name2
FROM cyto_edges
ORDER BY gene_name1;
# full 11527
# filtered: 1959

##### edges table
SELECT DISTINCT CONCAT( UPPER(gene_name1), " (", 'interacts with', ") ", UPPER(gene_name2) ) AS shared_name,
MAX(predicted_interaction_score) AS predicted_interaction_score, hpi_score, inttype,
iscentral, ismischnik, isstated, 
isswapped, isorganism, interactions1, interactions2,
`action`, a_is_acting, action_score,
pubmed_id, mischnik_publication, interactiontype
FROM cyto_edges
GROUP BY gene_name1, gene_name2
ORDER BY gene_name1;
# full 11527
# filtered 1959

###### nodes table
SELECT DISTINCT UPPER(gene_name) AS gene_name, rpkm_hsa, rpkm_mmu, rpkm_delta,
isdesi, isweyrich, ismischnik, iscentral, isnonortholog, 
neighbour_degree, connector, `TYPE`, description,
GROUP_CONCAT(species SEPARATOR " ; ") AS species,
interesting_proteins, isorganism
FROM cyto_nodes
#where rpkm_hsa = 0 and rpkm_mmu = 0 and iscentral = 0 and isdesi = 0
GROUP BY gene_name
ORDER BY gene_name;
# full 1811
# filtered 414