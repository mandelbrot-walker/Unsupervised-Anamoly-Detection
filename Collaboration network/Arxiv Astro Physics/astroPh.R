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

dg<-degree(g) # Calculation of Degree centrality
write.csv(dg, "dg_CA-AstroPh.csv")

btn<-betweenness(g) # Calculation of Betweenness centrality
write.csv(btn, "btn_CA-AstroPh.csv")

eig<-evcent(g)$vector # Calculation of Eigenvector centrality
write.csv(eig, "eig_CA-AstroPh.csv")

clsn<-closeness.latora(g) # Calculation of Closeness centrality
write.csv(clsn, "clsn_CA-AstroPh.csv")

pgr<-page_rank(g)$vector #Calculation of Page Rank centrality
write.csv(pgr, "pgr_CA-AstroPh.csv")

#katz<-katzcent(g) # Error in alpha >= maxEigenvalue : invalid comparison with complex values

#crsc<-crossclique(g) # Calculation of Cross-Clique centrality
#write.csv(crsc, "crsc_CA-AstroPh.csv")

#cntr<-centroid(g) #  Error: Graph is not strongly connected.
radial<-radiality(g) #  takes awhile
#clstrnk<-clusterrank(g)
denmnc<-dmnc(g)
lby<-lobby(g)
#library(linkcomm)
#comm<-communitycent(g) #  takes too long therefore stopped
lvg<-leverage(g)
#subg<-subgraph.centrality(g) #  takes a long time
topol<-topocoefficient(as.undirected(g))
ecc<-eccentricity(g)
gkp<-geokpath(g)
library(sna)
str<-calculate_centralities(as.directed(g), include = "Stress Centrality") #  takes a lot of memory 
infc<-calculate_centralities(as.directed(g), include = "Information Centrality") #  takes a long time
#mkc<-markovcent(g) #  takes a lot of memory and time therefore stopped
#entc<-entropy(g) #  takes too long therefore stopped
lbc<-local_bridging_centrality(g)

frm<-closeness.freeman(g) # Not calculatable as graphis not strongly connected

edge_connectivity(g) # Outputs 0

clstrnk[is.na(clstrnk)] <- 0

gr<-g # temporary variable gr

V(gr)$degree <- dg                               #  Degree centrality
V(gr)$eig <- eig                                 #  Eigenvector centrality
V(gr)$closeness <- clsn                          #  Closeness centrality
V(gr)$pagerank <- pgr                            #  Pagerank centrality
V(gr)$betweenness <- btn                         #  Vertex betweenness centrality
#V(gr)$crossclique <- crsc                        #  Crossclique centrality
V(gr)$hubs <- hub.score(g)$vector                #  Hub centrality
V(gr)$authorities <- authority.score(g)$vector   #  Authority centrality
V(gr)$radial<-radial
V(gr)$clusterrank<-clstrnk
V(gr)$dmnc<-denmnc
V(gr)$lobby<-lby
V(gr)$leverage<-lvg
V(gr)$subgraph<-subg
V(gr)$topologicalcoeff<-topol
V(gr)$eccentricity<-ecc
V(gr)$gdkpath<-gkp
V(gr)$stress<-unlist(str)
V(gr)$informationcent<-unlist(infc)
V(gr)$localbrigdecent<-lbc


centrality <- data.frame(row.names   = V(gr)$name,
                         degree      = V(gr)$degree,
                         eigenvector = V(gr)$eig,
                         closeness   = V(gr)$closeness,
                         pagerank    = V(gr)$pagerank,
                         betweenness = V(gr)$betweenness,
                         hubscore    = V(gr)$hubs,
                         authorities = V(gr)$authorities,
                         radiality   = V(gr)$radial,
                         clusterrank = V(gr)$clusterrank,
                         densitymnc  = V(gr)$dmnc,
                         lobby       = V(gr)$lobby,
                         leverage    = V(gr)$leverage,
                         #subgraph    = V(gr)$subgraph,
                         topologicalcoeff = V(gr)$topologicalcoeff,
                         eccentricity = V(gr)$eccentricity,
                         geodkpath   = V(gr)$gdkpath,
                         stress      = V(gr)$stress,
                         informationcent = V(gr)$informationcent,
                         localbridge = V(gr)$localbrigdecent
) #  Non-normalized centrality values

centrality <- centrality[order(row.names(centrality)),] #  These values are not normalized


head(centrality) #  check centrality variables

ndegree      = normalize(dg) 
neigenvector = normalize(eig) 
ncloseness   = normalize(clsn)
npagerank    = normalize(pgr)
#ncrossclique = normalize(crsc)
nbetweenness = normalize(btn)
nhubscore    = normalize(V(gr)$hubs)
nauthorities = normalize(V(gr)$authorities)
nradiality   = normalize(radial)
nclusterrank = normalize(clstrnk)
ndmnc        = normalize(denmnc)
nlobby       = normalize(lby)
nleverage    = normalize(lvg)
nsubgraph    = normalize(abs(subg))
ntopologicalcoeff = normalize(topol)
neccentricity = normalize(ecc)
ngeodkpath   = normalize(gkp)
nstress      = normalize(unlist(str))
ninformationcent = normalize(unlist(infc))
nlocalbridge = normalize(lbc)


ncentrality  <- data.frame(degree      = ndegree,
                           eigenvector = neigenvector,
                           closeness   = ncloseness,
                           pagerank    = npagerank,
                           #crossclique = ncrossclique,
                           betweenness = nbetweenness,
                           hubscore    = nhubscore,
                           authorities = nauthorities,
                           radiality   = nradiality,
                           clusterrank = nclusterrank,
                           dmnc        = ndmnc,
                           lobby       = nlobby,
                           leverage    = nleverage,
                           #subgraph    = nsubgraph,
                           topologicalcoeff = ntopologicalcoeff,
                           eccentricity = neccentricity,
                           geodkpath   = ngeodkpath,
                           stress      = nstress,
                           informationcent = ninformationcent,
                           localbridge = nlocalbridge
) #  normalized values 8 variables

rm(ndegree,neigenvector,ncloseness,npagerank,nbetweenness,nhubscore,nauthorities,nradiality,nclusterrank,ndmnc,
   nlobby,nleverage,nsubgraph,ntopologicalcoeff,neccentricity,ngeodkpath,nstress,ninformationcent,nlocalbridge)

#------------------------------------------Data loader and centrality calculation End-------------------------------------#

var8<-ncentrality
