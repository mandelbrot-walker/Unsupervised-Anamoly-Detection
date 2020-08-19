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
edges<-read.delim("Email-EuAll.txt",header = TRUE, sep = "\t")

g<-graph.data.frame(edges) #  graph data frame for igraph

transitivity(g) #  Check for cross clique calculation

cent<-proper_centralities(g)
# 
# calculate_centralities(g, include = cent[1:50])%>%
#   pca_centralities(scale.unit = TRUE, ncp = 50) # takes indefinite time

c<-c("Page Rank","Closeness centrality (Latora)","Degree Centrality","eigenvector centralities",
     "Kleinberg's authority centrality scores", "Kleinberg's hub centrality scores","Shortest-Paths Betweenness Centrality", 
     "Centroid value", "Radiality Centrality", "ClusterRank", "DMNC - Density of Maximum Neighborhood Component",
     "clustering coefficient", "Lobby Index (Centrality)", "Community Centrality", "Leverage Centrality",
     "Load Centrality", "subgraph centrality scores", "Topological Coefficient", "Diffusion Degree",
     "Eccentricity Centrality", "Geodesic K-Path Centrality", "Stress Centrality", "Information Centrality", 
     "Markov Centrality", "Entropy Centrality", "Local Bridging Centrality")
#cent[3,29,11,16,20,21,27,7,26,9,24,31,36,43,12,14,18,42,45,25,33,49])

calculate_centralities(g, include = c)%>%
  pca_centralities(scale.unit = TRUE, ncp = 50) # takes indefinite time


dg<-degree(g) # Calculation of Degree centrality
write.csv(dg, "dg_p2p_Gnutella04.csv")

btn<-betweenness(g) # Calculation of Betweenness centrality
write.csv(btn, "btn_p2p_Gnutella04.csv")

eig<-evcent(g)$vector # Calculation of Eigenvector centrality
write.csv(eig, "eig_p2p_Gnutella04.csv")

clsn<-closeness.latora(g) # Calculation of Closeness centrality
write.csv(clsn, "clsn_p2p_Gnutella04.csv")

pgr<-page_rank(g)$vector #Calculation of Page Rank centrality
write.csv(pgr, "pgr_p2p_Gnutella04.csv")

#katz<-katzcent(g) # Error in alpha >= maxEigenvalue : invalid comparison with complex values

#crsc<-crossclique(g) # Calculation of Cross-Clique centrality
#write.csv(crsc, "crsc_p2p_Gnutella04.csv")

#cntr<-centroid(g) #  Error: Graph is not strongly connected.
radial<-radiality(g) #  takes awhile
save.image(".Rdata",safe =F)
#clstrnk<-clusterrank(g)
denmnc<-dmnc(g)
save.image(".Rdata",safe =F)
lby<-lobby(g)
save.image(".Rdata",safe =F)
#library(linkcomm)
#comm<-communitycent(g) #  takes too long therefore stopped
lvg<-leverage(g)
save.image(".Rdata",safe =F)
#subg<-subgraph.centrality(g) #  takes a long time
topol<-topocoefficient(as.undirected(g))
save.image(".Rdata",safe =F)
ecc<-eccentricity(g)
save.image(".Rdata",safe =F)
gkp<-geokpath(g)
save.image(".Rdata",safe =F)
#library(sna)
#str<-calculate_centralities(g, include = "Stress Centrality") #  takes a lot of memory 
#infc<-calculate_centralities(g, include = "Information Centrality") #  takes a long time
#mkc<-markovcent(g) #  takes a lot of memory and time therefore stopped
#entc<-entropy(g) #  takes too long therefore stopped
lbc<-local_bridging_centrality(g)
save.image(".Rdata",safe =F)
frm<-closeness.freeman(g) # Not calculatable as graphis not strongly connected

edge_connectivity(g) # Outputs 0

clstrnk[is.na(clstrnk)] <- 0

cmp<-decompose.graph(g)
gkp<-c(0)

t<-cmp[[1]]
tmp<-geokpath(t, weights = NULL, mode = c("out"))
gkp[(length(gkp) + 1):(length(gkp) + length(tmp))]<-tmp

for (i in 2:length(cmp)) {
  t<-cmp[[i]]
  tmp<-geokpath(t)
  gkp[(length(gkp) + 1):(length(gkp) + length(tmp))]<-tmp
}
gkp<-gkp[-1]

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

#------------------------------------------Data loader and centrality calculation End-------------------------------------#

#------------------------------------------Write data start-------------------------------------#
write.csv(centrality,"Centrality_p2p_Gnutella04.csv")
write.csv(ncentrality,"NCentrality_p2p_Gnutella04.csv")

#------------------------------------------Write data end-------------------------------------#

#------------------------------------------Boxplot and corelation matrix start-------------------------------------#
M <- cor(ncentrality)

#plot correlation matrix
bmp("corrplot.bmp", width = 1280, height = 720)
c<-corrplot(M, method = "circle") # correlation matrix 
dev.off()

#plot boxplot
bmp("boxplot.bmp", width = 1280, height = 720)
boxplot(ndegree,neigenvector,ncloseness,npagerank,nbetweenness,nhubscore,nauthorities,nradiality,nclusterrank,ndmnc,
        nlobby,nleverage,nsubgraph,ntopologicalcoeff,neccentricity,ngeodkpath,nstress,ninformationcent,nlocalbridge,
        main = "Multiple boxplots for comparision",
        at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
        names = c("degree", "eigenvector", "closeness", "pagerank","betweenness","hubscore","authorities",
                  "radiality","clusterrank","dmnc","lobby","leverage","subgraph","topologicalcoeff","eccentricity",
                  "geodkpath","stress","informationcent","localbridge"),
        las = 2,
        col = c("orange","black"),
        border = "brown",
        horizontal = TRUE,
        notch = FALSE
) #multiple boxplot 
dev.off()
#------------------------------------------Boxplot and corelation matrix end-------------------------------------#
