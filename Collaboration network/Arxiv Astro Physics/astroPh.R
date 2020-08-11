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
edges<-read.delim("CA-AstroPh.txt",header = TRUE, sep = "\t")

g<-graph.data.frame(edges) #graph data frame for igraph

transitivity(g) #  Check for cross clique calculation

# cent<-proper_centralities(g)
# 
# calculate_centralities(g, include = cent[1:50])%>%
#   pca_centralities(scale.unit = TRUE, ncp = 50) # takes indefinite time

dg<-degree(g) # Calculation of Degree centrality
write.csv(dg, "dg_astroPh.csv")

btn<-betweenness(g) # Calculation of Betweenness centrality
write.csv(btn, "btn_astroPh.csv")

eig<-evcent(g)$vector # Calculation of Eigenvector centrality
write.csv(eig, "eig_astroPh.csv")

clsn<-closeness(g) # Calculation of Closeness centrality
write.csv(clsn, "clsn_astroPh.csv")

pgr<-page_rank(g)$vector #Calculation of Page Rank centrality
write.csv(pgr, "pgr_astroPh.csv")

#katz<-katzcent(g) # Error in alpha >= maxEigenvalue : invalid comparison with complex values

crsc<-crossclique(g) # Calculation of Cross-Clique centrality
write.csv(crsc, "crsc_astroPh.csv")

frm<-closeness.freeman(g) # Not calculatable as graphis not strongly connected

edge_connectivity(g) # Outputs 0