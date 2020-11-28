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
library(uwot)
library(clusternor)

#------------------------------------------Data loader and centrality calculation start-------------------------------------# 
edges<-read.delim("Email-EuAll.txt",header = TRUE, sep = "\t")


g<-graph.data.frame(edges) #graph data frame for igraph

transitivity(g) #  Check for cross clique calculation

#dg<-degree(g) # Calculation of Degree centrality
#btn<-betweenness(g) # Calculation of Betweenness centrality
#eig<-evcent(g)$vector # Calculation of Eigenvector centrality
clsn<-closeness.latora(g) # Calculation of Closeness centrality
#pgr<-page_rank(g)$vector #  Calculation of Page Rank centrality
auth<-authority.score(g)$vector 
hubs<-hub.score(g)$vector              #  Hub centrality
denmnc<-dmnc(g)
lby<-lobby(g)
lvg<-leverage(g)
ecc<-eccentricity(g)
lbc<-local_bridging_centrality(g)

cmp<-decompose.graph(g)
#infc<-c(0)
#infc[-1]


#cent<-proper_centralities(g)
# 
# calculate_centralities(g, include = cent[1:50])%>%
#   pca_centralities(scale.unit = TRUE, ncp = 50) # takes indefinite time

#c<-c("Page Rank","Closeness centrality (Latora)","Degree Centrality","eigenvector centralities",
# "Kleinberg's authority centrality scores", "Kleinberg's hub centrality scores","Shortest-Paths Betweenness Centrality", 
# "DMNC - Density of Maximum Neighborhood Component","Lobby Index (Centrality)", "Leverage Centrality",
# "Eccentricity Centrality", "Information Centrality","Local Bridging Centrality")
#cent[3,29,11,16,20,21,27,7,26,9,24,31,36,43,12,14,18,42,45,25,33,49])

#calculate_centralities(g, include = c)%>%
#  pca_centralities(scale.unit = TRUE, ncp = 50) # takes indefinite time

#katz<-katzcent(g) # Error in alpha >= maxEigenvalue : invalid comparison with complex values
#crsc<-crossclique(g) # Calculation of Cross-Clique centrality
#cntr<-centroid(g) #  Error: Graph is not strongly connected.
#radial<-radiality(g) #  takes awhile
#clstrnk<-clusterrank(g)
#library(linkcomm)
#comm<-communitycent(g) #  takes too long therefore stopped
#subg<-subgraph.centrality(g) #  takes a long time
#topol<-topocoefficient(as.undirected(g))
#gkp<-geokpath(g)
#library(sna)
#str<-calculate_centralities(g, include = "Stress Centrality") #  takes a lot of memory 
#mkc<-markovcent(g) #  takes a lot of memory and time therefore stopped
#entc<-entropy(g) #  takes too long therefore stopped
#frm<-closeness.freeman(g) # Not calculatable as graphis not strongly connected
#write.csv(dg, "dg_p2p_Gnutella04.csv")
#write.csv(btn, "btn_p2p_Gnutella04.csv")
#write.csv(eig, "eig_p2p_Gnutella04.csv")
#write.csv(clsn, "clsn_p2p_Gnutella04.csv")
#write.csv(pgr, "pgr_p2p_Gnutella04.csv")
#write.csv(crsc, "crsc_p2p_Gnutella04.csv")

#edge_connectivity(g) # Outputs 0
#clstrnk[is.na(clstrnk)] <- 0

#V(gr)$crossclique <- crsc                       #  Crossclique centrality
#V(gr)$radial<-radial
#V(gr)$clusterrank<-clstrnk
#V(gr)$subgraph<-subg
#V(gr)$topologicalcoeff<-topol
#V(gr)$gdkpath<-gkp
#V(gr)$stress<-unlist(str)

#radiality   = V(gr)$radial,
#clusterrank = V(gr)$clusterrank,
#subgraph    = V(gr)$subgraph,
#topologicalcoeff = V(gr)$topologicalcoeff,
#geodkpath   = V(gr)$gdkpath,
#stress      = V(gr)$stress,

#ncrossclique = normalize(crsc)
#nradiality   = normalize(radial)
#nclusterrank = normalize(clstrnk)
#nsubgraph    = normalize(abs(subg))
#ntopologicalcoeff = normalize(topol)
#ngeodkpath   = normalize(gkp)
#nstress      = normalize(unlist(str))

#crossclique = ncrossclique,
#radiality   = nradiality,
#clusterrank = nclusterrank,
#subgraph    = nsubgraph,
#topologicalcoeff = ntopologicalcoeff,
#geodkpath   = ngeodkpath,
#stress      = nstress,

gr<-g # temporary variable gr

V(gr)$degree <- dg                               #  Degree centrality
V(gr)$eig <- eig                                 #  Eigenvector centrality
V(gr)$pagerank <- pgr                            #  Pagerank centrality
V(gr)$authorities <- auth                        #  Authority centrality
V(gr)$hubs <- hubs                               #  Hub centrality
V(gr)$betweenness <- btn                         #  Vertex betweenness centrality
V(gr)$closeness <- clsn                          #  Closeness centrality
V(gr)$eccentricity<-ecc
V(gr)$dmnc<-denmnc
V(gr)$lobby<-lby
V(gr)$leverage<-lvg
V(gr)$localbrigdecent<-lbc


centrality <- data.frame(row.names   = V(gr)$name,
                         degree      = V(gr)$degree,
                         eigenvector = V(gr)$eig,
                         pagerank    = V(gr)$pagerank,
                         authorities = V(gr)$authorities,
                         hubscore    = V(gr)$hubs,
                         betweenness = V(gr)$betweenness,
                         closeness   = V(gr)$closeness,
                         eccentricity = V(gr)$eccentricity,
                         densitymnc  = V(gr)$dmnc,
                         lobby       = V(gr)$lobby,
                         leverage    = V(gr)$leverage,
                         localbridge = V(gr)$localbrigdecent
) #  Non-normalized centrality values

centrality <- centrality[order(row.names(centrality)),] #  These values are not normalized


head(centrality) #  check centrality variables

ndegree      = normalize(dg) 
neigenvector = normalize(eig) 
npagerank    = normalize(pgr)
nauthorities = normalize(auth)
nhubscore    = normalize(hubs)
nbetweenness = normalize(btn)
ncloseness   = normalize(clsn)
neccentricity = normalize(ecc)
ndmnc        = normalize(denmnc)
nlobby       = normalize(lby)
nleverage    = normalize(lvg)
nlocalbridge = normalize(lbc)


ncentrality  <- data.frame(degree      = ndegree,              #  1
                           eigenvector = neigenvector,         #  2
                           pagerank    = npagerank,            #  3
                           authorities = nauthorities,         #  4
                           hubscore    = nhubscore,            #  5
                           betweenness = nbetweenness,         #  6
                           closeness   = ncloseness,           #  7
                           eccentricity = neccentricity,       #  8
                           dmnc        = ndmnc,                #  9
                           lobby       = nlobby,               #  10
                           leverage    = nleverage,            #  11
                           localbridge = nlocalbridge          #  12
) #  normalized values 12 variables

#------------------------------------------Datasets-----------------------------------------------------------------------#

m1_all_var13<-ncentrality
m2_without_dg_var12<-within(ncentrality, rm(degree))
m3_without_comm_var7<-within(ncentrality, rm(eccentricity,dmnc,lobby,leverage,localbridge))
m4_without_nodes_var6<-within(ncentrality, rm(degree,eigenvector,pagerank,authorities,hubscore,betweenness,closeness))
m5_without_ranks_var9<-within(ncentrality, rm(eigenvector,pagerank,authorities,hubscore))
m6_without_dist_var9<-within(ncentrality, rm(betweenness,closeness,eccentricity))
m7_mix_match1_var6<-within(ncentrality, rm(eigenvector,closeness,dmnc,lobby,leverage,localbridge))
m8_mix_match2_var6<-within(ncentrality, rm(degree,authorities,hubscore,betweenness,eccentricity,leverage))
m9_mix_match3_var6<-within(ncentrality, rm(closeness,eigenvector,authorities,hubscore,betweenness,leverage))

m1<-m1_all_var13
m2<-m2_without_dg_var12
m3<-m3_without_comm_var7
m4<-m4_without_nodes_var6
m5<-m5_without_ranks_var9
m6<-m6_without_dist_var9
m7<-m7_mix_match1_var6
m8<-m8_mix_match2_var6
m9<-m9_mix_match3_var6

rm(m1_all_var13,m2_without_dg_var12,m3_without_comm_var7,m4_without_nodes_var6,m5_without_ranks_var9,m6_without_dist_var9,
   m7_mix_match1_var6,m8_mix_match2_var6,m9_mix_match3_var6)
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
boxplot(ndegree,neigenvector,npagerank,nauthorities,nhubscore,nbetweenness,ncloseness,neccentricity,
        ndmnc,nlobby,nleverage,nlocalbridge,
        main = "Multiple boxplots for comparision",
        at = c(1,2,3,4,5,6,7,8,9,10,11,12),
        names = c("degree","eigenvector","pagerank","authorities","hubscore","betweenness","closeness",
                  "eccentricity","dmnc","lobby","leverage","localbridge"),
        las = 2,
        col = c("orange","black"),
        border = "brown",
        horizontal = TRUE,
        notch = FALSE
) #multiple boxplot 
dev.off()

rm(ndegree,neigenvector,ncloseness,npagerank,nbetweenness,nhubscore,nauthorities,ndmnc,
   nlobby,nleverage,neccentricity,nlocalbridge)
#------------------------------------------Boxplot and corelation matrix end-------------------------------------#

#------------------------------------------PCA start-------------------------------------#

#principal component analysis

res.m1<-prcomp(scale(m1),center=TRUE) 
res.m2<-prcomp(scale(m2),center=TRUE) 
res.m3<-prcomp(scale(m3),center=TRUE) 
res.m4<-prcomp(scale(m4),center=TRUE) 
res.m5<-prcomp(scale(m5),center=TRUE) 
res.m6<-prcomp(scale(m6),center=TRUE) 
res.m7<-prcomp(scale(m7),center=TRUE) 
res.m8<-prcomp(scale(m8),center=TRUE) 
res.m9<-prcomp(scale(m9),center=TRUE) 

#show pca values
print(res.m1)
print(res.m2)
print(res.m3)
print(res.m4)
print(res.m5)
print(res.m6)
print(res.m7)
print(res.m8)
print(res.m9)

#  Results for Variables of m1
res.var <- get_pca_var(res.m1)#  13 variables 
res.var$coord          #  Coordinates
res.var$contrib        #  Contributions to the PCs
res.var$cos2           #  Quality of representation 

eig.m1 <- get_eigenvalue(res.m1) 
eig.m1
eig.m2 <- get_eigenvalue(res.m2) 
eig.m2
eig.m3 <- get_eigenvalue(res.m3) 
eig.m3
eig.m4 <- get_eigenvalue(res.m4) 
eig.m4
eig.m5 <- get_eigenvalue(res.m5) 
eig.m5
eig.m6 <- get_eigenvalue(res.m6) 
eig.m6
eig.m7 <- get_eigenvalue(res.m7) 
eig.m7
eig.m8 <- get_eigenvalue(res.m8) 
eig.m8
eig.m9 <- get_eigenvalue(res.m9) 
eig.m9

#-----------------------------------------------------PCA plots start-------------------------------------------------#

#  screeplot
p1<-fviz_eig(res.m1,addlabels = TRUE)
p2<-fviz_eig(res.m2,addlabels = TRUE)
p3<-fviz_eig(res.m3,addlabels = TRUE)
p4<-fviz_eig(res.m4,addlabels = TRUE)
p5<-fviz_eig(res.m5,addlabels = TRUE)
p6<-fviz_eig(res.m6,addlabels = TRUE)
p7<-fviz_eig(res.m7,addlabels = TRUE)
p8<-fviz_eig(res.m8,addlabels = TRUE)
p9<-fviz_eig(res.m9,addlabels = TRUE)

#  plot screeplot
bmp("screeplot.bmp", width = 1920, height = 1080)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3, ncol=3)
dev.off()

#  biplot
p1<-fviz_pca_biplot(res.m1, label ="var")
p2<-fviz_pca_biplot(res.m2, label ="var")
p3<-fviz_pca_biplot(res.m3, label ="var")
p4<-fviz_pca_biplot(res.m4, label ="var")
p5<-fviz_pca_biplot(res.m5, label ="var")
p6<-fviz_pca_biplot(res.m6, label ="var")
p7<-fviz_pca_biplot(res.m7, label ="var")
p8<-fviz_pca_biplot(res.m8, label ="var")
p9<-fviz_pca_biplot(res.m9, label ="var")

#  plot biplot
bmp("biplot.bmp", width = 1920, height = 1080)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3, ncol=3)
dev.off()

#  plot contribution to dimensions

p1d1<-fviz_contrib(res.m1, choice = "var", axes = 1)
p1d2<-fviz_contrib(res.m1, choice = "var", axes = 2)
p2d1<-fviz_contrib(res.m2, choice = "var", axes = 1)
p2d2<-fviz_contrib(res.m2, choice = "var", axes = 2)
p3d1<-fviz_contrib(res.m3, choice = "var", axes = 1)
p3d2<-fviz_contrib(res.m3, choice = "var", axes = 2)
p4d1<-fviz_contrib(res.m4, choice = "var", axes = 1)
p4d2<-fviz_contrib(res.m4, choice = "var", axes = 2)
p5d1<-fviz_contrib(res.m5, choice = "var", axes = 1)
p5d2<-fviz_contrib(res.m5, choice = "var", axes = 2)
p6d1<-fviz_contrib(res.m6, choice = "var", axes = 1)
p6d2<-fviz_contrib(res.m6, choice = "var", axes = 2)
p7d1<-fviz_contrib(res.m7, choice = "var", axes = 1)
p7d2<-fviz_contrib(res.m7, choice = "var", axes = 2)
p8d1<-fviz_contrib(res.m8, choice = "var", axes = 1)
p8d2<-fviz_contrib(res.m8, choice = "var", axes = 2)
p9d1<-fviz_contrib(res.m9, choice = "var", axes = 1)
p9d2<-fviz_contrib(res.m9, choice = "var", axes = 2)


bmp("dimension contribution.bmp", width = 1920, height = 1080)
grid.arrange(p1d1,p1d2,p2d1,p2d2,p3d1,p3d2,p4d1,p4d2,p5d1,p5d2,p6d1,p6d2,p7d1,p7d2,p8d1,p8d2,p9d1,p9d2, nrow = 5, ncol=4)
dev.off()

rm(p1d1,p1d2,p2d1,p2d2,p3d1,p3d2,p4d1,p4d2,p5d1,p5d2,p6d1,p6d2,p7d1,p7d2,p8d1,p8d2,p9d1,p9d2) #  remove contribution plot variables


# Color by contributions to the PC uses cos^2 value

p1<-fviz_pca_var(res.m1,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

p2<-fviz_pca_var(res.m2,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

p3<-fviz_pca_var(res.m3,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

p4<-fviz_pca_var(res.m4,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

p5<-fviz_pca_var(res.m5,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

p6<-fviz_pca_var(res.m6,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

p7<-fviz_pca_var(res.m7,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

p8<-fviz_pca_var(res.m8,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

p9<-fviz_pca_var(res.m9,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

#plot contribution wise pca
bmp("cos2 contribution.bmp", width = 1920, height = 1080)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3, ncol=3)
dev.off()

rm(eig.m1,eig.m2,eig.m3,eig.m4,eig.m5,eig.m6,eig.m7,eig.m8,eig.m9)

rm(p1,p2,p3,p4,p5,p6,p7,p8,p9)
rm(res.m1,res.m2,res.m3,res.m4,res.m5,res.m6,res.m7,res.m8,res.m9)

#-----------------------------------------------------PCA plots ends-------------------------------------------------#

#------------------------------------------PCA end-------------------------------------#

#------------------------------------------------Tsne---------------------------------------------------#

ncen_tr<-transpose(ncentrality) #  transpose ncentrality for tsne colors
ncen_tr<-data.frame(names = c("degree", "eigenvector","pagerank","authorities","hubscore","betweenness","closeness",
                              "eccentricity","dmnc","lobby","leverage","localbridge"),ncen_tr) #y label

#-------------------------------tsne parameters----------------------------#

prx1<-50
prx2<-100
prx3<-30
prx4<-43
th1<-0.20
th2<-0.10
th3<-0.50
th4<-0.10
mit1<-2000
mit2<-1500
mit3<-1000
mit4<-1500
nthr<-11

#-------------------------------------------------------tsne model 1 start--------------------------------------------#

#  m1
set.seed(32)  
tsne_model_1_m1 = Rtsne(m1, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m2
set.seed(32)  
tsne_model_1_m2 = Rtsne(m2, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m3
set.seed(32)  
tsne_model_1_m3 = Rtsne(m3, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m3
set.seed(32)  
tsne_model_1_m4 = Rtsne(m4, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m5
set.seed(32)  
tsne_model_1_m5 = Rtsne(m5, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m6
set.seed(32)  
tsne_model_1_m6 = Rtsne(m6, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m7
set.seed(32)  
tsne_model_1_m7 = Rtsne(m7, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m8
set.seed(32)  
tsne_model_1_m8 = Rtsne(m8, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m9
set.seed(32)  
tsne_model_1_m9 = Rtsne(m9, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)

#------------------------------------------------TSNE model 1 plots------------------------------------------------#

bmp("tsne_model_1_m1.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m1$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_1_m2.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m2$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_1_m3.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m3$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_1_m4.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m4$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_1_m5.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m5$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_1_m6.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m6$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_1_m7.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m7$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_1_m8.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m8$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_1_m9.bmp", width = 1920, height = 1080)
plot(tsne_model_1_m9$Y,col=factor(ncen_tr$names), asp=1)
dev.off()

#-------------------------------------------------------tsne model 1 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 2 start--------------------------------------------#

#  m1
set.seed(31)  
tsne_model_2_m1 = Rtsne(m1, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m2
set.seed(31)  
tsne_model_2_m2 = Rtsne(m2, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m3
set.seed(31)  
tsne_model_2_m3 = Rtsne(m3, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m4
set.seed(31)  
tsne_model_2_m4 = Rtsne(m4, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m5
set.seed(31)  
tsne_model_2_m5 = Rtsne(m5, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m6
set.seed(31)  
tsne_model_2_m6 = Rtsne(m6, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m7
set.seed(31)  
tsne_model_2_m7 = Rtsne(m7, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)

#  m8
set.seed(31)  
tsne_model_2_m8 = Rtsne(m8, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)

#  m9
set.seed(31)  
tsne_model_2_m9 = Rtsne(m9, check_duplicates=FALSE, pca=TRUE, perplexity=prx2, theta=th2, dims=2, max_iter = mit2,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)

#------------------------------------------------TSNE model 2 plots------------------------------------------------#

bmp("tsne_model_2_m1.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m1$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_2_m2.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m2$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_2_m3.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m3$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_2_m4.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m4$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_2_m5.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m5$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_2_m6.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m6$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_2_m7.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m7$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_2_m8.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m8$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_2_m9.bmp", width = 1920, height = 1080)
plot(tsne_model_2_m9$Y,col=factor(ncen_tr$names), asp=1)
dev.off()

#-------------------------------------------------------tsne model 2 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 3 start--------------------------------------------#

#  m1
set.seed(30)  
tsne_model_3_m1 = Rtsne(m1, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m2
set.seed(30)  
tsne_model_3_m2 = Rtsne(m2, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m3
set.seed(30)  
tsne_model_3_m3 = Rtsne(m3, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m4
set.seed(30)  
tsne_model_3_m4 = Rtsne(m4, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m5
set.seed(30)  
tsne_model_3_m5 = Rtsne(m5, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m6
set.seed(30)  
tsne_model_3_m6 = Rtsne(m6, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m7
set.seed(30)  
tsne_model_3_m7 = Rtsne(m7, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m8
set.seed(30)  
tsne_model_3_m8 = Rtsne(m8, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m9
set.seed(30)  
tsne_model_3_m9 = Rtsne(m9, check_duplicates=FALSE, pca=TRUE, perplexity=prx3, theta=th3, dims=2, max_iter = mit3,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)

#------------------------------------------------TSNE model 3 plots------------------------------------------------#

bmp("tsne_model_3_m1.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m1$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_3_m2.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m2$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_3_m3.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m3$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_3_m4.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m4$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_3_m5.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m5$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_3_m6.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m6$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_3_m7.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m7$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_3_m8.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m8$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_3_m9.bmp", width = 1920, height = 1080)
plot(tsne_model_3_m9$Y,col=factor(ncen_tr$names), asp=1)
dev.off()

#-------------------------------------------------------tsne model 3 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 4 start--------------------------------------------#

#  m1
set.seed(29)  
tsne_model_4_m1 = Rtsne(m1, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m2
set.seed(29)  
tsne_model_4_m2 = Rtsne(m2, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m3
set.seed(29)  
tsne_model_4_m3 = Rtsne(m3, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m4
set.seed(29)  
tsne_model_4_m4 = Rtsne(m4, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m5
set.seed(29)  
tsne_model_4_m5 = Rtsne(m5, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m6
set.seed(29)  
tsne_model_4_m6 = Rtsne(m6, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m7
set.seed(29)  
tsne_model_4_m7 = Rtsne(m7, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m8
set.seed(29)  
tsne_model_4_m8 = Rtsne(m8, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
#  m9
set.seed(29)  
tsne_model_4_m9 = Rtsne(m9, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=2, max_iter = mit4,
                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)

#------------------------------------------------TSNE model 4 plots------------------------------------------------#

bmp("tsne_model_4_m1.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m1$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_4_m2.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m2$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_4_m3.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m3$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_4_m4.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m4$Y,col=factor(ncen_tr$names), asp=1)
dev.off()


bmp("tsne_model_4_m5.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m5$Y,col=factor(ncen_tr$names), asp=1)
dev.off()

bmp("tsne_model_4_m6.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m6$Y,col=factor(ncen_tr$names), asp=1)
dev.off()

bmp("tsne_model_4_m7.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m7$Y,col=factor(ncen_tr$names), asp=1)
dev.off()

bmp("tsne_model_4_m8.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m8$Y,col=factor(ncen_tr$names), asp=1)
dev.off()

bmp("tsne_model_4_m9.bmp", width = 1920, height = 1080)
plot(tsne_model_4_m9$Y,col=factor(ncen_tr$names), asp=1)
dev.off()

#-------------------------------------------------------tsne model 4 end--------------------------------------------#