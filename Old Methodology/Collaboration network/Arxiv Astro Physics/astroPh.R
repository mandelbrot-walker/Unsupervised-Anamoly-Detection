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
cmp<-decompose.graph(g)
str<-c(0)

t<-cmp[[1]]
tmp<-calculate_centralities(as.directed(t), include = "Information Centrality")
str[(length(str) + 1):(length(str) + length(unlist(tmp)))]<-unlist(tmp)

str2<-str[-1]

for (i in 2:length(cmp)) {
        if(toString(which_loop(cmp[[i]])) == "TRUE"){
                print(i)
                tmp <- 0
                str2[(length(str2) + 1):(length(str2) + length(tmp))]<-tmp
                i = i+1
                print(i)
        } else{
                t<-cmp[[i]]
                tmp<-calculate_centralities(as.directed(t), include = "Information Centrality")
                str2[(length(str2) + 1):(length(str2) + length(unlist(tmp)))]<-unlist(tmp)
        }
}
str<-str[-1]

for (j in 1:length(cmp)) {
        result<-which_loop(cmp[[j]])
        # if (result[i]) {
        #         print(j)
        # } else{
        #         which_loop(cmp[[j]])
        # }
}


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
#V(gr)$clusterrank<-clstrnk
V(gr)$dmnc<-denmnc
V(gr)$lobby<-lby
V(gr)$leverage<-lvg
#V(gr)$subgraph<-subg
V(gr)$topologicalcoeff<-topol
V(gr)$eccentricity<-ecc
V(gr)$gdkpath<-gkp
#V(gr)$stress<-unlist(str)
#V(gr)$informationcent<-unlist(infc)
V(gr)$localbrigdecent<-lbc


centrality <- data.frame(row.names   = V(gr)$name,
                         degree      = V(gr)$degree,
                         eigenvector = V(gr)$eig,
                         closeness   = V(gr)$closeness,
                         pagerank    = V(gr)$pagerank,
                         betweenness = V(gr)$betweenness,
                         hubscore    = V(gr)$hubs,
                         authorities = V(gr)$authorities,
                         #radiality   = V(gr)$radial,
                         #clusterrank = V(gr)$clusterrank,
                         densitymnc  = V(gr)$dmnc,
                         lobby       = V(gr)$lobby,
                         leverage    = V(gr)$leverage,
                         #subgraph    = V(gr)$subgraph,
                         #topologicalcoeff = V(gr)$topologicalcoeff,
                         eccentricity = V(gr)$eccentricity,
                         geodkpath   = V(gr)$gdkpath,
                         #stress      = V(gr)$stress,
                         #informationcent = V(gr)$informationcent,
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
#nradiality   = normalize(radial)
#nclusterrank = normalize(clstrnk)
ndmnc        = normalize(denmnc)
nlobby       = normalize(lby)
nleverage    = normalize(lvg)
#nsubgraph    = normalize(abs(subg))
#ntopologicalcoeff = normalize(topol)
neccentricity = normalize(ecc)
ngeodkpath   = normalize(gkp)
#nstress      = normalize(unlist(str))
#ninformationcent = normalize(unlist(infc))
nlocalbridge = normalize(lbc)


ncentrality  <- data.frame(degree      = ndegree,
                           eigenvector = neigenvector,
                           closeness   = ncloseness,
                           pagerank    = npagerank,
                           #crossclique = ncrossclique,
                           betweenness = nbetweenness,
                           hubscore    = nhubscore,
                           authorities = nauthorities,
                           #radiality   = nradiality,
                           #clusterrank = nclusterrank,
                           dmnc        = ndmnc,
                           lobby       = nlobby,
                           leverage    = nleverage,
                           #subgraph    = nsubgraph,
                           #topologicalcoeff = ntopologicalcoeff,
                           eccentricity = neccentricity,
                           geodkpath   = ngeodkpath,
                           #stress      = nstress,
                           #informationcent = ninformationcent,
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
boxplot(ndegree,neigenvector,ncloseness,npagerank,nbetweenness,nhubscore,nauthorities,ndmnc,
        nlobby,nleverage,neccentricity,ngeodkpath,nlocalbridge,
        main = "Multiple boxplots for comparision",
        at = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
        names = c("degree", "eigenvector", "closeness", "pagerank","betweenness","hubscore","authorities",
                  "dmnc","lobby","leverage","eccentricity",
                  "geodkpath","localbridge"),
        las = 2,
        col = c("orange","black"),
        border = "brown",
        horizontal = TRUE,
        notch = FALSE
) #multiple boxplot 
dev.off()
#------------------------------------------Boxplot and corelation matrix end-------------------------------------#

#------------------------------------------PCA start-------------------------------------#

#principal component analysis

res.pca8<-prcomp(scale(var8),center=TRUE) #  8 variables

#show pca values
print(res.pca8)

p1<-fviz_eig(res.pca8,addlabels = TRUE)#Scree plot 8 var

#plot screeplot
bmp("screeplot.bmp", width = 1920, height = 1080)
grid.arrange(p1, nrow = 1)
dev.off()

p1<-fviz_pca_biplot(res.pca8, label ="var")#biplot 8 var

#plot biplot
bmp("biplot.bmp", width = 1920, height = 1080)
grid.arrange(p1, nrow = 1)
dev.off()


#plot contribution to dimensions
p5d1<-fviz_contrib(res.pca8, choice = "var", axes = 1)
p5d2<-fviz_contrib(res.pca8, choice = "var", axes = 2)

bmp("dimension contribution.bmp", width = 1920, height = 1080)
grid.arrange(p5d1,p5d2, nrow = 1, ncol = 2)
dev.off()

rm(p5d1,p5d2) #  remove contribution plot variables


# Color by contributions to the PC uses cos^2 value

p1<-fviz_pca_var(res.pca8,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
) #  8 variables

#plot contribution wise pca
#plot biplot
bmp("cos2 contribution.bmp", width = 1920, height = 1080)
grid.arrange(p1, nrow = 1)
dev.off()

eig.val8 <- get_eigenvalue(res.pca8) 
eig.val8


#  Results for Variables
res.var <- get_pca_var(res.pca8)#  8 variables 
res.var$coord          #  Coordinates
res.var$contrib        #  Contributions to the PCs
res.var$cos2           #  Quality of representation 

#------------------------------------------PCA end-------------------------------------#
