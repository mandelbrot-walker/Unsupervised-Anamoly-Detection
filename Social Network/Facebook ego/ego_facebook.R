library(igraph)
library(centiserve)
library(tidyverse)
library(factoextra)
library(Rtsne)
library(Plasmidprofiler)
library(MASS)
library(data.table)
library(corrplot)
library(tibble)
library(caret)
library(plyr)
library(gridExtra) 
library(CINNA)
library(ClusterR)
library(mclust)
library(kohonen)
library(kernlab)
library(Rdimtools)
library(scatterplot3d)
library(rgl)
library(fpc)

#------------------------------------------Data loader and centrality calculation start-------------------------------------# 
#edges<-read_csv("facebook_combined.txt") #Reading in ego_facebook data
edges<-read.delim("facebook_combined.txt",header = TRUE, sep = " ")

g<-graph.data.frame(edges) #graph data frame for igraph

transitivity(g)
clique_num(g)
library(gRbase)
tg<-igraph.to.graphNEL(g)
gc()
gclq<-getCliques(tg)

dg<-degree(g) # Calculation of Degre centrality
write.csv(dg, "dg_ego_facebook.csv")

btn<-betweenness(g) # Calculation of Betweeness centrality
write.csv(btn, "btn_ego_facebook.csv")

eig<-evcent(g)$vector # Calculation of Eigenvector centrality
write.csv(eig, "eig_ego_facebook.csv")

clsn<-closeness(g) # Calculation of Closeness centrality
write.csv(clsn, "clsn_ego_facebook.csv")

pgr<-page_rank(g)$vector #Calculation of Page Rank centrality
write.csv(pgr, "pgr_ego_facebook.csv")

katz<-katzcent(g) # Error: Graph is not loop free #need to recheck for this dataset

crsc<-crossclique(g) # Calculation of Cross-Clique centrality
write.csv(crsc, "crsc_ego_facebook.csv")

frm<-closeness.freeman(g) # Not calculatable as graphis not strongly connected

edge_connectivity(g) # Outputs 0
        

