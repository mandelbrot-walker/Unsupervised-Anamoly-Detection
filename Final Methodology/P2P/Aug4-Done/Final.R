#--------------index-------------------------#

#Lines       Code segment
# 64-198     Centrality calculation
# 171-250    dataset preparation
# 255-280    Correlation matrix and boxplot
# 282        PCA start 
# 284-330    PCA calculation
# 334-364    screeplot and biplot
# 368-461    contribution plot
# 463        PCA end
# 465        TSNE start 
# 473-485    tsne parameter macros
# 487-525    tsne model 1
# 526-572    tsne model 1 plots
# 576-616    tsne model 2
# 617-663    tsne model 2 plots
# 667-705    tsne model 3
# 708-750    tsne model 3 plots
# 756-794    tsne model 4
# 795-835    tsne model 4 plots
# 837        TSNE end
# 839        Spherical kmeans start
# 842-1037   cluster calculation
# 1040-1215  tsne m1 and kmeans
# 1217-1315  tsne m2 and kmeans
# 1317-1415  tsne m3 and kmeans
# 1417-1515  tsne m4 and kmeans
# 1517       tsne and kmeans end 
# 1519       Sammon's map
# 1521-1611  Sammon's map on models 
# 1613       Sammon's map End
# 1614       End of code
#--------------index-------------------------#

#--------------libraries--------------------#
library(igraph) #  centralities
library(centiserve) #  centralities
library(factoextra)  # fviz_cluster()
library(Rtsne) #  TSNE
library(Plasmidprofiler) #  normalize()
library(data.table) #  transpose()
library(corrplot) #  correlation plot
library(gridExtra) #  multiplot
library(CINNA) #  centralities
#library(mclust) #  gmm
library(Rdimtools) #  sammon's map
library(fpc) #  calinhara()
library(clusternor) # skmeans  
# library(scatterplot3d) #  3D plot
# library(rgl) #  3D plot
# library(kohonen) #  SOM
# library(kernlab) # kpca
# library(ClusterR) # GMM check later
# library(tidyverse)
# library(MASS)
# library(tibble)
# library(caret)
# library(plyr)
# library(uwot) #  UMAP
#--------------libraries--------------------#

#------------------------------------------Data loader and centrality calculation start-------------------------------------# 
edges<-read.delim("p2p-Gnutella04.txt",header = TRUE, sep = "\t")

g<-graph.data.frame(edges) #graph data frame for igraph

transitivity(g) #  Check for cross clique calculation

dg<-degree(g) # Calculation of Degree centrality
btn<-betweenness(g) # Calculation of Betweenness centrality
eig<-evcent(g)$vector # Calculation of Eigenvector centrality
clsn<-closeness.latora(g) # Calculation of Closeness centrality
pgr<-page_rank(g)$vector #  Calculation of Page Rank centrality
auth<-authority.score(g)$vector 
hubs<-hub.score(g)$vector              #  Hub centrality
denmnc<-dmnc(g)
lby<-lobby(g)
lvg<-leverage(g)
ecc<-eccentricity(g)
#infc<-calculate_centralities(g, include = "Information Centrality") #  takes a long time
lbc<-local_bridging_centrality(g)

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
#V(gr)$informationcent<-unlist(infc)
V(gr)$eccentricity<-ecc                          #  Eccentricity centrality  
V(gr)$dmnc<-denmnc                               #  Density of Maximum Neighborhood Component Centrality  
V(gr)$lobby<-lby                                 #  Lobby Centrality 
V(gr)$leverage<-lvg                              #  Leverage Centrality  
V(gr)$localbrigdecent<-lbc                       #  Local Bridging Centrality


centrality <- data.frame(row.names   = V(gr)$name,
                         degree      = V(gr)$degree,
                         eigenvector = V(gr)$eig,
                         pagerank    = V(gr)$pagerank,
                         authorities = V(gr)$authorities,
                         hubscore    = V(gr)$hubs,
                         betweenness = V(gr)$betweenness,
                         closeness   = V(gr)$closeness,
                         #informationcent = V(gr)$informationcent,
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
#ninformationcent = normalize(unlist(infc))
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
                           #informationcent = ninformationcent, #  8
                           eccentricity = neccentricity,       #  9
                           dmnc        = ndmnc,                #  10
                           lobby       = nlobby,               #  11
                           leverage    = nleverage,            #  12
                           localbridge = nlocalbridge          #  13
) #  normalized values 13 variables

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
prx4<-71
th1<-0.20
th2<-0.10
th3<-0.50
th4<-0.10
mit1<-2000
mit2<-1500
mit3<-1000
mit4<-1500
nthr<-7

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

#---------------------------------------------Spherical kmeans start------------------------------------------------------#

#---------------------------------------------------------clusters--------------------------------------#
#--------------clusters using different k values m1
ck2<-Skmeans(data=as.matrix(m1),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m1),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m1),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m1),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m1),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()  # garbage collection is used for stack imbalance warning. run gc() more than once if the warning persists
ck7<-Skmeans(data=as.matrix(m1),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck8<-Skmeans(data=as.matrix(m1),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck9<-Skmeans(data=as.matrix(m1),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck10<-Skmeans(data=as.matrix(m1),centers=10,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m1,ck2$cluster),digits=2) #  15207.12  Highest
round(calinhara(m1,ck3$cluster),digits=3) #  12288.01
round(calinhara(m1,ck4$cluster),digits=4) #  10007.39
round(calinhara(m1,ck5$cluster),digits=5) #  9231.524
round(calinhara(m1,ck6$cluster),digits=6) #  7951.821
round(calinhara(m1,ck7$cluster),digits=7) #  8351.8
round(calinhara(m1,ck8$cluster),digits=8) #  8019.982
round(calinhara(m1,ck9$cluster),digits=9) #  7610.551
round(calinhara(m1,ck10$cluster),digits=10) #  7284.249

#--------------clusters using different k values m2

ck2<-Skmeans(data=as.matrix(m2),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m2),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m2),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m2),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m2),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck7<-Skmeans(data=as.matrix(m2),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m2,ck2$cluster),digits=2) #  15275.51 Highest
round(calinhara(m2,ck3$cluster),digits=3) #  12261.01
round(calinhara(m2,ck4$cluster),digits=4) #  9913.505
round(calinhara(m2,ck5$cluster),digits=5) #  9177.806
round(calinhara(m2,ck6$cluster),digits=6) #  7913.423
round(calinhara(m2,ck7$cluster),digits=7) #  8281.966

#--------------clusters using different k values m3

ck2<-Skmeans(data=as.matrix(m3),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m3),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m3),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m3),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m3),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck7<-Skmeans(data=as.matrix(m3),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m3,ck2$cluster),digits=2) #  4264.96
round(calinhara(m3,ck3$cluster),digits=3) #  2896.255
round(calinhara(m3,ck4$cluster),digits=4) #  3756.562
round(calinhara(m3,ck5$cluster),digits=5) #  5066.883  Highest
round(calinhara(m3,ck6$cluster),digits=6) #  4665.315
round(calinhara(m3,ck7$cluster),digits=7) #  4139.602


#--------------clusters using different k values m4

ck2<-Skmeans(data=as.matrix(m4),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m4),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m4),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m4),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m4),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck7<-Skmeans(data=as.matrix(m4),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck8<-Skmeans(data=as.matrix(m4),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m4,ck2$cluster),digits=2) #  21664.47  Highest
round(calinhara(m4,ck3$cluster),digits=3) #  16172.21
round(calinhara(m4,ck4$cluster),digits=4) #  11952.2  
round(calinhara(m4,ck5$cluster),digits=5) #  10157.21
round(calinhara(m4,ck6$cluster),digits=6) #  9828.79
round(calinhara(m4,ck7$cluster),digits=7) #  8690.346
round(calinhara(m4,ck8$cluster),digits=8) #  7934.881

#--------------clusters using different k values m5

ck2<-Skmeans(data=as.matrix(m5),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m5),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m5),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m5),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m5),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck7<-Skmeans(data=as.matrix(m5),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck8<-Skmeans(data=as.matrix(m5),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m5,ck2$cluster),digits=2) #  18489.77  Highest
round(calinhara(m5,ck3$cluster),digits=3) #  16133.64
round(calinhara(m5,ck4$cluster),digits=4) #  13112.88
round(calinhara(m5,ck5$cluster),digits=5) #  11411.81
round(calinhara(m5,ck6$cluster),digits=6) #  10703.69
round(calinhara(m5,ck7$cluster),digits=7) #  10250.08
round(calinhara(m5,ck8$cluster),digits=8) #  9295.133

#--------------clusters using different k values m6

ck2<-Skmeans(data=as.matrix(m6),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m6),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m6),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m6),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m6),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m6,ck2$cluster),digits=2) #  5696.55  Highest
round(calinhara(m6,ck3$cluster),digits=3) #  5206.248
round(calinhara(m6,ck4$cluster),digits=4) #  4234.737
round(calinhara(m6,ck5$cluster),digits=5) #  3772.685
round(calinhara(m6,ck6$cluster),digits=6) #  3349.911 

#--------------clusters using different k values m7

ck2<-Skmeans(data=as.matrix(m7),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m7),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m7),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m7),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m7),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m7,ck2$cluster),digits=2) #  4227.35  Highest
round(calinhara(m7,ck3$cluster),digits=3) #  4279.856
round(calinhara(m7,ck4$cluster),digits=4) #  3399.06  
round(calinhara(m7,ck5$cluster),digits=5) #  3747.733
round(calinhara(m7,ck6$cluster),digits=6) #  3519.926

#--------------clusters using different k values m8

ck2<-Skmeans(data=as.matrix(m8),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m8),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m8),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m8),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m8),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m8,ck2$cluster),digits=2) #  18373.12  Highest
round(calinhara(m8,ck3$cluster),digits=3) #  15258.98
round(calinhara(m8,ck4$cluster),digits=4) #  15253.42
round(calinhara(m8,ck5$cluster),digits=5) #  13077.41
round(calinhara(m8,ck6$cluster),digits=6) #  10217.41

#--------------clusters using different k values m9

ck2<-Skmeans(data=as.matrix(m9),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck3<-Skmeans(data=as.matrix(m9),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck4<-Skmeans(data=as.matrix(m9),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck5<-Skmeans(data=as.matrix(m9),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck6<-Skmeans(data=as.matrix(m9),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m9,ck2$cluster),digits=2) #  17202.42  Highest
round(calinhara(m9,ck3$cluster),digits=3) #  12415.04
round(calinhara(m9,ck4$cluster),digits=4) #  9323.659
round(calinhara(m9,ck5$cluster),digits=5) #  8433.265
round(calinhara(m9,ck6$cluster),digits=6) #  7677.063

rm(ck2,ck3,ck4,ck5,ck6,ck7,ck8,ck9,ck10)
#-------------kmeans on dataset and cluster onto TSNE start-------------------------------------------------#

#-----------------------------------------Kmeans for tsne model 1------------------------------------------#

ckm1<-Skmeans(data=as.matrix(m1),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm2<-Skmeans(data=as.matrix(m2),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm3<-Skmeans(data=as.matrix(m3),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm4<-Skmeans(data=as.matrix(m4),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm5<-Skmeans(data=as.matrix(m5),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm6<-Skmeans(data=as.matrix(m6),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm7<-Skmeans(data=as.matrix(m7),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm8<-Skmeans(data=as.matrix(m8),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm9<-Skmeans(data=as.matrix(m9),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()  # garbage collection is used for stack imbalance warning. run gc() more than once if the warning persists

ck1<-kmeans(m1, 3, iter.max = 20, nstart = 25,
            algorithm = c("Hartigan-Wong"), trace=FALSE) #  dummy kmeans
ck2<-ck1
ck3<-ck1
ck4<-ck1
ck5<-ck1
ck6<-ck1
ck7<-ck1
ck8<-ck1
ck9<-ck1

#  m1
bmp("tsne_model1_m1_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m1$Y), col = ckm1$cluster)
dev.off()


ck1$cluster<-ckm1$cluster
ck1$centers<-ckm1$centers
ck1$size<-ckm1$size
ck1$iter<-ckm1$iters

p1 <- fviz_cluster(ck1, geom = "point",  data = as.data.frame(tsne_model_1_m1$Y)) 

bmp("tsne_model1_m1_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m2
bmp("tsne_model1_m2_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m2$Y), col = ckm2$cluster)
dev.off()

ck2$cluster<-ckm2$cluster
ck2$centers<-ckm2$centers
ck2$size<-ckm2$size
ck2$iter<-ckm2$iters

p1 <- fviz_cluster(ck2, geom = "point",  data = as.data.frame(tsne_model_1_m2$Y)) 

bmp("tsne_model1_m2_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m3
bmp("tsne_model1_m3_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m3$Y), col = ckm3$cluster)
dev.off()

ck3$cluster<-ckm3$cluster
ck3$centers<-ckm3$centers
ck3$size<-ckm3$size
ck3$iter<-ckm3$iters

p1 <- fviz_cluster(ck3, geom = "point",  data = as.data.frame(tsne_model_1_m3$Y)) 

bmp("tsne_model1_m3_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m4
bmp("tsne_model1_m4_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m4$Y), col = ckm4$cluster)
dev.off()

ck4$cluster<-ckm4$cluster
ck4$centers<-ckm4$centers
ck4$size<-ckm4$size
ck4$iter<-ckm4$iters

p1 <- fviz_cluster(ck4, geom = "point",  data = as.data.frame(tsne_model_1_m4$Y)) 

bmp("tsne_model1_m4_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m5
bmp("tsne_model1_m5_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m5$Y), col = ckm5$cluster)
dev.off()

ck5$cluster<-ckm5$cluster
ck5$centers<-ckm5$centers
ck5$size<-ckm5$size
ck5$iter<-ckm5$iters

p1 <- fviz_cluster(ck5, geom = "point",  data = as.data.frame(tsne_model_1_m5$Y)) 

bmp("tsne_model1_m5_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m6
bmp("tsne_model1_m6_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m6$Y), col = ckm6$cluster)
dev.off()

ck6$cluster<-ckm6$cluster
ck6$centers<-ckm6$centers
ck6$size<-ckm6$size
ck6$iter<-ckm6$iters

p1 <- fviz_cluster(ck6, geom = "point",  data = as.data.frame(tsne_model_1_m6$Y)) 

bmp("tsne_model1_m6_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m7
bmp("tsne_model1_m7_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m7$Y), col = ckm7$cluster)
dev.off()

ck7$cluster<-ckm7$cluster
ck7$centers<-ckm7$centers
ck7$size<-ckm7$size
ck7$iter<-ckm7$iters

p1 <- fviz_cluster(ck7, geom = "point",  data = as.data.frame(tsne_model_1_m7$Y)) 

bmp("tsne_model1_m7_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m8
bmp("tsne_model1_m8_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m8$Y), col = ckm8$cluster)
dev.off()

ck8$cluster<-ckm8$cluster
ck8$centers<-ckm8$centers
ck8$size<-ckm8$size
ck8$iter<-ckm8$iters

p1 <- fviz_cluster(ck8, geom = "point",  data = as.data.frame(tsne_model_1_m8$Y)) 

bmp("tsne_model1_m8_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m9
bmp("tsne_model1_m9_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_1_m9$Y), col = ckm9$cluster)
dev.off()

ck9$cluster<-ckm9$cluster
ck9$centers<-ckm9$centers
ck9$size<-ckm9$size
ck9$iter<-ckm9$iters

p1 <- fviz_cluster(ck9, geom = "point",  data = as.data.frame(tsne_model_1_m9$Y)) 

bmp("tsne_model1_m9_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#-------------------------------------------------Kmeans for tsne model 2--------------------------------#

bmp("tsne_model2_m1_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m1$Y), col = ckm1$cluster)
dev.off()

p1 <- fviz_cluster(ck1, geom = "point",  data = as.data.frame(tsne_model_2_m1$Y)) 

bmp("tsne_model2_m1_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m2
bmp("tsne_model2_m2_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m2$Y), col = ckm2$cluster)
dev.off()

p1 <- fviz_cluster(ck2, geom = "point",  data = as.data.frame(tsne_model_2_m2$Y)) 

bmp("tsne_model2_m2_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m3
bmp("tsne_model2_m3_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m3$Y), col = ckm3$cluster)
dev.off()

p1 <- fviz_cluster(ck3, geom = "point",  data = as.data.frame(tsne_model_2_m3$Y)) 

bmp("tsne_model2_m3_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m4
bmp("tsne_model2_m4_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m4$Y), col = ckm4$cluster)
dev.off()

p1 <- fviz_cluster(ck4, geom = "point",  data = as.data.frame(tsne_model_2_m4$Y)) 

bmp("tsne_model2_m4_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m5
bmp("tsne_model2_m5_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m5$Y), col = ckm5$cluster)
dev.off()

p1 <- fviz_cluster(ck5, geom = "point",  data = as.data.frame(tsne_model_2_m5$Y)) 

bmp("tsne_model2_m5_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m6
bmp("tsne_model2_m6_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m6$Y), col = ckm6$cluster)
dev.off()

p1 <- fviz_cluster(ck6, geom = "point",  data = as.data.frame(tsne_model_2_m6$Y)) 

bmp("tsne_model2_m6_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m7
bmp("tsne_model2_m7_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m7$Y), col = ckm7$cluster)
dev.off()

p1 <- fviz_cluster(ck7, geom = "point",  data = as.data.frame(tsne_model_2_m7$Y)) 

bmp("tsne_model2_m7_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m8
bmp("tsne_model2_m8_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m8$Y), col = ckm8$cluster)
dev.off()

p1 <- fviz_cluster(ck8, geom = "point",  data = as.data.frame(tsne_model_2_m8$Y)) 

bmp("tsne_model2_m8_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m9
bmp("tsne_model2_m9_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_2_m9$Y), col = ckm9$cluster)
dev.off()

p1 <- fviz_cluster(ck9, geom = "point",  data = as.data.frame(tsne_model_2_m9$Y)) 

bmp("tsne_model2_m9_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#--------------------------------------------------Kmeans for tsne model 3------------------------------#

bmp("tsne_model3_m1_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m1$Y), col = ckm1$cluster)
dev.off()

p1 <- fviz_cluster(ck1, geom = "point",  data = as.data.frame(tsne_model_3_m1$Y)) 

bmp("tsne_model3_m1_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m2
bmp("tsne_model3_m2_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m2$Y), col = ckm2$cluster)
dev.off()

p1 <- fviz_cluster(ck2, geom = "point",  data = as.data.frame(tsne_model_3_m2$Y)) 

bmp("tsne_model3_m2_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m3
bmp("tsne_model3_m3_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m3$Y), col = ckm3$cluster)
dev.off()

p1 <- fviz_cluster(ck3, geom = "point",  data = as.data.frame(tsne_model_3_m3$Y)) 

bmp("tsne_model3_m3_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m4
bmp("tsne_model3_m4_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m4$Y), col = ckm4$cluster)
dev.off()

p1 <- fviz_cluster(ck4, geom = "point",  data = as.data.frame(tsne_model_3_m4$Y)) 

bmp("tsne_model3_m4_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m5
bmp("tsne_model3_m5_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m5$Y), col = ckm5$cluster)
dev.off()

p1 <- fviz_cluster(ck5, geom = "point",  data = as.data.frame(tsne_model_3_m5$Y)) 

bmp("tsne_model3_m5_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m6
bmp("tsne_model3_m6_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m6$Y), col = ckm6$cluster)
dev.off()

p1 <- fviz_cluster(ck6, geom = "point",  data = as.data.frame(tsne_model_3_m6$Y)) 

bmp("tsne_model3_m6_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m7
bmp("tsne_model3_m7_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m7$Y), col = ckm7$cluster)
dev.off()

p1 <- fviz_cluster(ck7, geom = "point",  data = as.data.frame(tsne_model_3_m7$Y)) 

bmp("tsne_model3_m7_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m8
bmp("tsne_model3_m8_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m8$Y), col = ckm8$cluster)
dev.off()

p1 <- fviz_cluster(ck8, geom = "point",  data = as.data.frame(tsne_model_3_m8$Y)) 

bmp("tsne_model3_m8_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m9
bmp("tsne_model3_m9_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_3_m9$Y), col = ckm9$cluster)
dev.off()

p1 <- fviz_cluster(ck9, geom = "point",  data = as.data.frame(tsne_model_3_m9$Y)) 

bmp("tsne_model3_m9_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#-------------------------------------------------Kmeans for tsne model 4----------------------------------#

bmp("tsne_model4_m1_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m1$Y), col = ckm1$cluster)
dev.off()

p1 <- fviz_cluster(ck1, geom = "point",  data = as.data.frame(tsne_model_4_m1$Y)) 

bmp("tsne_model4_m1_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m2
bmp("tsne_model4_m2_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m2$Y), col = ckm2$cluster)
dev.off()

p1 <- fviz_cluster(ck2, geom = "point",  data = as.data.frame(tsne_model_4_m2$Y)) 

bmp("tsne_model4_m2_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m3
bmp("tsne_model4_m3_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m3$Y), col = ckm3$cluster)
dev.off()

p1 <- fviz_cluster(ck3, geom = "point",  data = as.data.frame(tsne_model_4_m3$Y)) 

bmp("tsne_model4_m3_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m4
bmp("tsne_model4_m4_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m4$Y), col = ckm4$cluster)
dev.off()

p1 <- fviz_cluster(ck4, geom = "point",  data = as.data.frame(tsne_model_4_m4$Y)) 

bmp("tsne_model4_m4_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m5
bmp("tsne_model4_m5_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m5$Y), col = ckm5$cluster)
dev.off()

p1 <- fviz_cluster(ck5, geom = "point",  data = as.data.frame(tsne_model_4_m5$Y)) 

bmp("tsne_model4_m5_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m6
bmp("tsne_model4_m6_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m6$Y), col = ckm6$cluster)
dev.off()

p1 <- fviz_cluster(ck6, geom = "point",  data = as.data.frame(tsne_model_4_m6$Y)) 

bmp("tsne_model4_m6_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m7
bmp("tsne_model4_m7_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m7$Y), col = ckm7$cluster)
dev.off()

p1 <- fviz_cluster(ck7, geom = "point",  data = as.data.frame(tsne_model_4_m7$Y)) 

bmp("tsne_model4_m7_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m8
bmp("tsne_model4_m8_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m8$Y), col = ckm8$cluster)
dev.off()

p1 <- fviz_cluster(ck8, geom = "point",  data = as.data.frame(tsne_model_4_m8$Y)) 

bmp("tsne_model4_m8_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#  m9
bmp("tsne_model4_m9_kmeans.bmp", width = 1980, height = 1280)
plot(as.data.frame(tsne_model_4_m9$Y), col = ckm9$cluster)
dev.off()

p1 <- fviz_cluster(ck9, geom = "point",  data = as.data.frame(tsne_model_4_m9$Y)) 

bmp("tsne_model4_m9_kmeans_ch.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#-----------------------kmeans on dataset and cluster onto TSNE end------------------------------------------------------#

#------------Sammon's map----------------#

x<-as.matrix(m1) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m1.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=factor(ncen_tr$names), main="m1")
plot(sammon$Y, pch=19, col=ckm1$cluster, main="m1 kmeans")
par(opar)
dev.off()

x<-as.matrix(m2) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m2.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:12, main="m2")
plot(sammon$Y, pch=19, col=ckm2$cluster, main="m2 kmeans")
par(opar)
dev.off()

x<-as.matrix(m3) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m3.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:7, main="m3")
plot(sammon$Y, pch=19, col=ckm3$cluster, main="m3 kmeans")
par(opar)
dev.off()

x<-as.matrix(m4) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m4.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:6, main="m4")
plot(sammon$Y, pch=19, col=ckm4$cluster, main="m4 kmeans")
par(opar)
dev.off()

x<-as.matrix(m5) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m5.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:9, main="m5")
plot(sammon$Y, pch=19, col=ckm5$cluster, main="m5 kmeans")
par(opar)
dev.off()

x<-as.matrix(m6) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m6.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:9, main="m6")
plot(sammon$Y, pch=19, col=ckm6$cluster, main="m6 kmeans")
par(opar)
dev.off()

x<-as.matrix(m7) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m7.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:6, main="m7")
plot(sammon$Y, pch=19, col=ckm7$cluster, main="m7 kmeans")
par(opar)
dev.off()

x<-as.matrix(m8) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m8.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:6, main="m8")
plot(sammon$Y, pch=19, col=ckm8$cluster, main="m8 kmeans")
par(opar)
dev.off()

x<-as.matrix(m9) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m9.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:6, main="m9")
plot(sammon$Y, pch=19, col=ckm9$cluster, main="m9 kmeans")
par(opar)
dev.off()

rm(x,opar)

#------------Sammon's map----------------#

cluster1_degree_m1<-as.data.frame(m1$degree[ck1$cluster ==1]) #  Better than subset
colnames(cluster1_degree_m1)[1] <- "degree"

denseplot <- function(data, p_title) {
        ggplot(data, aes(x=degree)) + 
                geom_histogram(aes(y=..density..), colour="black", fill="white",binwidth = .1,size=.1)  + 
                geom_density(alpha=.2, fill="#FF6666") + 
                geom_vline(aes(xintercept=mean(degree)), color="blue", linetype="dashed", size=1) +
                xlim(range(m1$degree)) +
                ggtitle(p_title) + theme_classic() +
                theme(plot.title = element_text(hjust = 0.5)) #  histogram of SR with density and mean
} 
p1<-denseplot(cluster1_degree_m1, "degree of cluster 1") 

cluster2_degree_m1<-as.data.frame(m1$degree[ck1$cluster ==2]) #  Better than subset
colnames(cluster2_degree_m1)[1] <- "degree"
p2<-denseplot(cluster2_degree_m1, "degree of cluster 2")

grid.arrange(p1,p2, nrow = 2, ncol=1)
rm(p1,p2)
rm(cluster1_degree_m1,cluster2_degree_m1)

#----------------------------------------Tests----------------------------------------------------------------------#
# library(MASS)
# library(plotly)
# 
# x <- tsne_model_4_m1$Y[,1]
# y <- tsne_model_4_m1$Y[,2]
# den3d <- kde2d(x, y)
# 
# plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()



# fit<-glm(m1$degree~.,family=gaussian,data=m1)
# summary(fit)
# par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
# plot(fit)
# 
# fit<-glm(m2$degree~.,family=gaussian,data=m1)
# summary(fit)
# par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
# plot(fit)
# 
# fviz_cluster(cc, geom = "point",  data = as.data.frame(m1)) 
# 
# km2<-Skmeans(data=as.matrix(m1),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# km3<-Skmeans(data=as.matrix(m1),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# km<-Skmeans(data=as.matrix(m1),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# km4<-Skmeans(data=as.matrix(m1),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# km6<-Skmeans(data=as.matrix(m1),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# 
# round(calinhara(m1,km2$cluster),digits=2)
# round(calinhara(m1,km3$cluster),digits=3)
# round(calinhara(m1,km$cluster),digits=4)
# round(calinhara(m1,km4$cluster),digits=5)
# round(calinhara(m1,km6$cluster),digits=6)
# 
# plot(tsne_model_1_m1$Y,col=km2$cluster,asp=1)
# 
# cc$cluster<-as.data.frame(spc@.Data)
# cc$centers<-as.data.frame(spc@centers)
# 
# fviz_cluster(cc, geom = "point",  data = as.data.frame(tsne_model_1_m1$Y)) 
# 
# plot(tsne_model_1_m1$Y, col = spc@.Data)
# 
# x1<-as.matrix(m1)
# 
# kmat<-kernelMatrix(rbfdot(sigma = 0.75), x1)
# 
# spc<-specc(x1, centers=3,
#            kernel = "rbfdot", kpar = "automatic", 
#            nystrom.red = TRUE, nystrom.sample = dim(x1)[1]/6,
#            iterations = 50, mod.sample = 0.75, na.action = na.omit)
# 
# library(HDclassif)
# 
# hdgmm<-hddc(m1, K = 1:10, model = c("ALL"), threshold = 0.45,
#             criterion = "bic", com_dim = 2, itermax = 50, eps = 0.001,
#             algo = "EM", d_select = "Cattell", init = "param", show = getHDclassif.show(),scaling = TRUE,
#             min.individuals = 10, noise.ctrl = 1e-08, mc.cores = 5,
#             nb.rep = 2, keepAllRes = FALSE, d_max = 50, subset = Inf)
# 
# hdgmmc<-hdgmm$class
# 
# plot(tsne_model_1_m1$Y, col = hdgmm$class)
# 
# table(hdgmm$class)
# 
# plot(hdgmm, method = "BIC")
# plot(tsne_model_1_m1$Y, col = xyMclust1$classification)
# 
# m1tr<-transpose(m1)
# label<-names(m1)
# 
# tsne_model_1_m1_3d = Rtsne(m1tr, check_duplicates=FALSE, pca=TRUE, perplexity=4, theta=th1, dims=2, max_iter = mit1,
#                             verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
# 
# plot(tsne_model_1_m1_3d$Y,col=factor(label),asp=1)
# 
# tsne_model_chk = Rtsne(m1[,-7], check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=2, max_iter = mit1,
#                        verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
# plot(tsne_model_chk$Y,col=factor(ncen_tr$names),asp=1)
# 
# library(ClusterR)
# 
# opt_gmm = Optimal_Clusters_GMM(m1, max_clusters = 13, criterion = "BIC", 
#                                
#                                dist_mode = "eucl_dist", seed_mode = "random_subset",
#                                
#                                km_iter = 10, em_iter = 10, var_floor = 1e-10, 
#                                
#                                plot_data = T)

#----------------------------------------Tests----------------------------------------------------------------------#

#------------------------------------------Old code-------------------------------------------------------------#

#------------------GMM-----------------------------#

##  m1

# set.seed(10)
# subset <- sample(1:nrow(m1), 500)
# 
# summary(mclustBIC(as.matrix(m1), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                   initialization=list(subset = subset),
#                   control = emControl(), 
#                   warn = mclust.options("warn"),
#                   verbose = TRUE))
# summary(mclustICL(as.matrix(m1), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                   initialization=list(subset = subset),
#                   control = emControl(), 
#                   warn = mclust.options("warn"),
#                   verbose = TRUE))
# 
# xyMclust1 <- Mclust(as.matrix(m1), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust1, parameters = F)

# plot(xyMclust1, what=c("classification"))
# plot(xyMclust1, "density")
# plot(xyMclust1, what=c("uncertainty"))
# plot(xyMclust1, what=c("BIC"))

# p1<-fviz_mclust(xyMclust1,what = c("classification"), geom = "point",ellipse.type = "norm",palette = "jco" )
# p2<-fviz_mclust(xyMclust1,what = c("uncertainty"),ellipse.type = "norm",palette = "jco" )
# 
# bmp("GMM_fit.bmp", width = 1440, height = 620)
# grid.arrange(p1,p2, nrow = 1)
# dev.off()

##  m2
# xyMclust2 <- Mclust(as.matrix(m2), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust2, parameters = F)
# 
# ##  m3
# xyMclust3 <- Mclust(as.matrix(m3), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust3, parameters = F)
# 
# ##  m4
# xyMclust4 <- Mclust(as.matrix(m4), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust4, parameters = F)
# 
# ##  m5
# xyMclust5 <- Mclust(as.matrix(m5), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust5, parameters = F)
# 
# ##  m6
# xyMclust6 <- Mclust(as.matrix(m6), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust6, parameters = F)
# 
# ##  m7
# xyMclust7 <- Mclust(as.matrix(m7), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust7, parameters = F)
# 
# ##  m8
# xyMclust8 <- Mclust(as.matrix(m8), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust8, parameters = F)
# 
# ##  m9
# xyMclust9 <- Mclust(as.matrix(m9), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust9, parameters = F)
# 
# 
# xyMclust9 <- Mclust(as.matrix(m9), prior = priorControl(functionName="defaultPrior", shrinkage=0.1),
#                     initialization=list(subset = subset),
#                     control = emControl(), 
#                     warn = mclust.options("warn"),
#                     verbose = TRUE)
# summary(xyMclust9, parameters = F)

#------------------GMM-----------------------------#


#-------------------------3d plotting-------------------------------------#
# set.seed(32)  
# tsne_model_1_m1_3d = Rtsne(m1, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=3, max_iter = mit1,
#                            verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
# 
# set.seed(32)  
# tsne_model_1_m9_3d = Rtsne(m9, check_duplicates=FALSE, pca=TRUE, perplexity=prx1, theta=th1, dims=3, max_iter = mit1,
#                            verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = nthr)
# 
# plot(tsne_model_1_m1_3d$Y,col=ckm1$cluster, asp=1)
# 
# plot(tsne_model_1_m1$Y,col=cc$cluster,asp=1)
# 
# tempo<-as.data.frame(res.m9$x)
# 
# plot(tempo$PC1,tempo$PC2,col=ckm1$cluster, asp=1)
# 
# plot3d(x=tsne_model_1_m1_3d$Y[,1],y=tsne_model_1_m1_3d$Y[,2],z=tsne_model_1_m1_3d$Y[,3],
#        col=cc$cluster,
#        type="s",radius=0.5)
# 
# plot3d(x=tempo$PC1,y=tempo$PC2,z=tempo$PC3,
#        col=1:10876,
#        type="s",radius=0.1)
# 
# plot3d(x=tsne_model_1_m9_3d$Y[,1],y=tsne_model_1_m9_3d$Y[,2],z=tsne_model_1_m9_3d$Y[,3],
#        col=ckm1$cluster,
#        type="s",radius=0.5)
#-------------------------------------------------------------------------------#


#------------------------different kmeans---------------------#

# library(sparcl)
# 
# kout <- KMeansSparseCluster(as.matrix(m1),K=4,wbounds=seq(1.3,4,len=8))
# print(kout)
# ks<-KMeansSparseCluster(as.matrix(m1), K=4, wbounds = kout$bestw, nstart = 20, silent =
#                       TRUE, maxiter=25, centers=NULL)
# 
# print(ks)
# plot(ks)
# 
# library(RSKC)
# rkm<-RSKC(as.matrix(m1), ncl=4, alpha=0.5, L1 = 1.5, nstart = 25, 
#      silent=FALSE, scaling = FALSE, correlation = FALSE)
# 
# plot(tsne_model_1_m1$Y,col=rkm$labels,asp=1)
# 
# cc$cluster<-rkm$labels
# 
# fviz_cluster(cc, geom = "point",  data = as.data.frame(tsne_model_1_m1$Y), asp=1.5) 


#------------------------------------------------------------#



#---------------------------------------------kneeplot----------------------------------------------------------#

##  m1
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m1, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m1.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()
# 
# ##  m2
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m2, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m2.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()
# 
# ##  m3
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m3, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m3.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()
# 
# ##  m4
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m4, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m4.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()
# 
# ##  m5
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m5, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m5.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()
# 
# ##  m6
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m6, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m6.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()
# 
# ##  m7
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m7, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m7.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()
# 
# ##  m8
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m8, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m8.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()
# 
# ##  m9
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(m9, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# bmp("kmeans_kneeplot_m9.bmp", width = 1280, height = 720)
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares") #kneeplot
# dev.off()

#---------------------------------------------kneeplot----------------------------------------------------------#
#---------------------------------------------kmeans------------------------------------------------------------#
#  m1
# ck1<-kmeans(m1, 3, iter.max = 20, nstart = 25,
#            algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# ck2<-kmeans(m1, 4, iter.max = 20, nstart = 25,
#            algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m1, 5, iter.max = 20, nstart = 25,
#            algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m1, 6, iter.max = 20, nstart = 25,
#            algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m1,ck1$cluster),digits=3) #  7733.985
# round(calinhara(m1,ck2$cluster),digits=4) #  7985.291 highest
# round(calinhara(m1,ck3$cluster),digits=5) #  6959.628
# round(calinhara(m1,ck4$cluster),digits=6) #  6348.964
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m1) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m1) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m1) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m1) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m1.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  m2
# ck1<-kmeans(m2, 3, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck2<-kmeans(m2, 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m2, 5, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m2, 6, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m2,ck1$cluster),digits=3) #  7703.294
# round(calinhara(m2,ck2$cluster),digits=4) #  8037.741 highest
# round(calinhara(m2,ck3$cluster),digits=5) #  7001.544
# round(calinhara(m2,ck4$cluster),digits=6) #  6370.787
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m2) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m2) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m2) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m2) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m2.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  m3
# ck1<-kmeans(m3, 3, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck2<-kmeans(m3, 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m3, 5, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m3, 6, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m3,ck1$cluster),digits=3) #  7649.803
# round(calinhara(m3,ck2$cluster),digits=4) #  7698.748 
# round(calinhara(m3,ck3$cluster),digits=5) #  7738.866
# round(calinhara(m3,ck4$cluster),digits=6) #  7787.33 highest
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m3) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m3) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m3) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m3) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m3.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  m4
# ck1<-kmeans(m4, 3, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck2<-kmeans(m4, 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m4, 5, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m4, 6, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m4,ck1$cluster),digits=3) #  9030.761
# round(calinhara(m4,ck2$cluster),digits=4) #  10646.54 highest
# round(calinhara(m4,ck3$cluster),digits=5) #  9468.393
# round(calinhara(m4,ck4$cluster),digits=6) #  8854.677
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m4) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m4) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m4) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m4) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m4.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  m5
# ck1<-kmeans(m5, 3, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck2<-kmeans(m5, 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m5, 5, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m5, 6, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m5,ck1$cluster),digits=3) #  8841.458
# round(calinhara(m5,ck2$cluster),digits=4) #  9615.578 highest
# round(calinhara(m5,ck3$cluster),digits=5) #  8563.109
# round(calinhara(m5,ck4$cluster),digits=6) #  7967.682
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m5) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m5) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m5) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m5) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m5.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  m6
# ck1<-kmeans(m6, 3, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck2<-kmeans(m6, 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m6, 5, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m6, 6, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m6,ck1$cluster),digits=3) #  16450.43 highest
# round(calinhara(m6,ck2$cluster),digits=4) #  14770.27 
# round(calinhara(m6,ck3$cluster),digits=5) #  14848.32
# round(calinhara(m6,ck4$cluster),digits=6) #  14466.01
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m6) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m6) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m6) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m6) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m6.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  m7
# ck1<-kmeans(m7, 3, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck2<-kmeans(m7, 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m7, 5, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m7, 6, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m7,ck1$cluster),digits=3) #  8412.904 
# round(calinhara(m7,ck2$cluster),digits=4) #  8478.609 highest
# round(calinhara(m7,ck3$cluster),digits=5) #  7935.222
# round(calinhara(m7,ck4$cluster),digits=6) #  8114.244
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m7) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m7) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m7) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m7) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m7.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  m8
# ck1<-kmeans(m8, 3, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck2<-kmeans(m8, 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m8, 5, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m8, 6, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m8,ck1$cluster),digits=3) #  18717.68 highest
# round(calinhara(m8,ck2$cluster),digits=4) #  17385.23 
# round(calinhara(m8,ck3$cluster),digits=5) #  16981.55
# round(calinhara(m8,ck4$cluster),digits=6) #  16352.07
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m8) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m8) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m8) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m8) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m8.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  m9
# ck1<-kmeans(m9, 3, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck2<-kmeans(m9, 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck3<-kmeans(m9, 5, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# ck4<-kmeans(m9, 6, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m9,ck1$cluster),digits=3) #  14642.38 highest
# round(calinhara(m9,ck2$cluster),digits=4) #  12906.56 
# round(calinhara(m9,ck3$cluster),digits=5) #  11824.88
# round(calinhara(m9,ck4$cluster),digits=6) #  11635.55
# 
# 
# #plot of clusters
# p1 <- fviz_cluster(ck1, geom = "point", data = m9) + ggtitle("k = 3")
# p2 <- fviz_cluster(ck2, geom = "point", data = m9) + ggtitle("k = 4")
# p3 <- fviz_cluster(ck3, geom = "point", data = m9) + ggtitle("k = 5")
# p4 <- fviz_cluster(ck4, geom = "point", data = m9) + ggtitle("k = 6")
# 
# bmp("kmeans_pca_m9.bmp", width = 1920, height = 1280)
# grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()

#  Full dataset
# c_ncentrality<-kmeans(ncentrality, 4, iter.max = 20, nstart = 25,
#                       algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# bmp("kmeans_ncentrality_k4.bmp", width = 1920, height = 1280)
# plot(ncentrality, col = c_ncentrality$cluster)
# points(c_ncentrality$centers, col = 1:8, pch = 8)
# dev.off()

# tc1<-kmeans(as.data.frame(tsne_model_1_m1$Y), 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# fviz_cluster(tc1, geom = "point", data = as.data.frame(tsne_model_1_m1$Y))
# 
# tc7<-kmeans(as.data.frame(tsne_model_1_m7$Y), 4, iter.max = 20, nstart = 25,
#             algorithm = c("Hartigan-Wong"), trace=FALSE)
# fviz_cluster(tc7, geom = "point", data = as.data.frame(tsne_model_1_m7$Y))

# #-------------------------------------------dbscan start-----------------------------------------------------------------#
# 
# #---------------------------------------calculating h start------------------------------------------#
# 
# bmp("dbscan_kneeplot_m1.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m1, k =  2)
# abline(h = 0.165, lty = 2)
# dev.off()
# 
# bmp("dbscan_kneeplot_m2.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m2, k =  2)
# abline(h = 0.165, lty = 2)
# dev.off()
# 
# bmp("dbscan_kneeplot_m2.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m3, k =  2)
# abline(h = 0.07, lty = 2)
# dev.off()
# 
# bmp("dbscan_kneeplot_m4.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m4, k =  2)
# abline(h = 0.075, lty = 2)
# dev.off()
# 
# bmp("dbscan_kneeplot_m5.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m5, k =  2)
# abline(h = 0.12, lty = 2)
# dev.off()
# 
# bmp("dbscan_kneeplot_m6.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m6, k =  2)
# abline(h = 0.08, lty = 2)
# dev.off()
# 
# bmp("dbscan_kneeplot_m7.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m7, k =  2)
# abline(h = 0.06, lty = 2)
# dev.off()
# 
# bmp("dbscan_kneeplot_m8.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m8, k =  2)
# abline(h = 0.05, lty = 2)
# dev.off()
# 
# bmp("dbscan_kneeplot_m9.bmp", width = 841, height = 477)
# dbscan::kNNdistplot(m9, k =  2)
# abline(h = 0.045, lty = 2)
# dev.off()
# #---------------------------------------calculating h end------------------------------------------#
# 
# #---------------------------------------plotting dbscan start--------------------------------------------------#
# 
# res.db <- dbscan::dbscan(m1, 0.165, 2)
# gc()
# bmp("dbscan_m1.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m1, geom = "point")
# dev.off()
# 
# res.db <- dbscan::dbscan(m2, 0.165, 2)
# gc()
# bmp("dbscan_m2.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m2, geom = "point")
# dev.off()
# 
# res.db <- dbscan::dbscan(m3, 0.07, 2)
# gc()
# bmp("dbscan_m3.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m3, geom = "point")
# dev.off()
# 
# res.db <- dbscan::dbscan(m4, 0.075, 2)
# gc()
# bmp("dbscan_m4.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m4, geom = "point")
# dev.off()
# 
# res.db <- dbscan::dbscan(m5, 0.12, 2)
# gc()
# bmp("dbscan_m5.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m5, geom = "point")
# dev.off()
# 
# res.db <- dbscan::dbscan(m6, 0.08, 2)
# gc()
# bmp("dbscan_m6.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m6, geom = "point")
# dev.off()
# 
# res.db <- dbscan::dbscan(m7, 0.06, 2)
# gc()
# bmp("dbscan_m7.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m7, geom = "point")
# dev.off()
# 
# res.db <- dbscan::dbscan(m8, 0.05, 2)
# gc()
# bmp("dbscan_m8.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m8, geom = "point")
# dev.off()
# 
# res.db <- dbscan::dbscan(m9, 0.045, 2)
# gc()
# bmp("dbscan_m9.bmp", width = 1980, height = 1280)
# fviz_cluster(res.db, m9, geom = "point")
# dev.off()
# 
# rm(res.db)
# 
# #---------------------------------------plotting dbscan end--------------------------------------------------#
# 
# #-----------------------------------------------dbscan end-----------------------------------------------------------------#


#----------------UMAP----------------------------#

#  parameters

# ngb<-90
# lr<-1.1
# 
# #  Umap m1
# umap_calc <- umap(m1, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m1 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# #  Umap plot
# 
# ggplot(umap_m1, aes(
#   x = UMAP1, y = UMAP2,
#   col = UMAP1
# )) +
#   geom_point()
# 
# bmp("umap_m1.bmp", width = 1920, height = 1280)
# plot(umap_m1$UMAP1,umap_m1$UMAP2, col=factor(umap_m1$UMAP1))
# dev.off()
# 
# bmp("umap_m1_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm1, geom = "point",  data = umap_m1) 
# dev.off()
# 
# #  Umap m2
# umap_calc <- umap(m2, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m2 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# bmp("umap_m2.bmp", width = 1920, height = 1280)
# plot(umap_m2$UMAP1,umap_m2$UMAP2, col=factor(umap_m2$UMAP1))
# dev.off()
# 
# bmp("umap_m2_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm2, geom = "point",  data = umap_m2) 
# dev.off()
# 
# #  Umap m3
# umap_calc <- umap(m3, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m3 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# bmp("umap_m3.bmp", width = 1920, height = 1280)
# plot(umap_m3$UMAP1,umap_m3$UMAP2, col=factor(umap_m3$UMAP1))
# dev.off()
# 
# bmp("umap_m3_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm3, geom = "point",  data = umap_m3) 
# dev.off()
# 
# #  Umap m4
# umap_calc <- umap(m4, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m4 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# bmp("umap_m4.bmp", width = 1920, height = 1280)
# plot(umap_m4$UMAP1,umap_m4$UMAP2, col=factor(umap_m4$UMAP1))
# dev.off()
# 
# bmp("umap_m4_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm4, geom = "point",  data = umap_m4) 
# dev.off()
# 
# #  Umap m5
# umap_calc <- umap(m5, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m5 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# bmp("umap_m5.bmp", width = 1920, height = 1280)
# plot(umap_m5$UMAP1,umap_m5$UMAP2, col=factor(umap_m5$UMAP1))
# dev.off()
# 
# bmp("umap_m5_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm5, geom = "point",  data = umap_m5) 
# dev.off()
# 
# #  Umap m6
# umap_calc <- umap(m6, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m6 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# bmp("umap_m6.bmp", width = 1920, height = 1280)
# plot(umap_m6$UMAP1,umap_m6$UMAP2, col=factor(umap_m6$UMAP1))
# dev.off()
# 
# bmp("umap_m6_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm6, geom = "point",  data = umap_m6) 
# dev.off()
# 
# #  Umap m7
# umap_calc <- umap(m7, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m7 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# bmp("umap_m7.bmp", width = 1920, height = 1280)
# plot(umap_m7$UMAP1,umap_m7$UMAP2, col=factor(umap_m7$UMAP1))
# dev.off()
# 
# bmp("umap_m7_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm7, geom = "point",  data = umap_m7) 
# dev.off()
# 
# #  Umap m8
# umap_calc <- umap(m8, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m8 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# bmp("umap_m8.bmp", width = 1920, height = 1280)
# plot(umap_m8$UMAP1,umap_m8$UMAP2, col=factor(umap_m8$UMAP1))
# dev.off()
# 
# bmp("umap_m8_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm8, geom = "point",  data = umap_m8) 
# dev.off()
# 
# #  Umap m9
# umap_calc <- umap(m9, n_neighbors = ngb, n_components = 2,metric = "cosine",
#                   learning_rate = lr, scale = T, init = "spectral",
#                   bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
#                   search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
#                   ret_nn = T, n_threads = nthr, verbose = getOption("verbose",TRUE),
#                   grain_size = 5,  min_dist = 10, spread = 50 )
# 
# umap_m9 <- data.frame(
#   UMAP1 = umap_calc$embedding[, 1],
#   UMAP2 = umap_calc$embedding[, 2]
# )
# 
# bmp("umap_m9.bmp", width = 1920, height = 1280)
# plot(umap_m9$UMAP1,umap_m9$UMAP2, col=factor(umap_m9$UMAP1))
# dev.off()
# 
# bmp("umap_m9_km.bmp", width = 1920, height = 1280)
# fviz_cluster(ckm9, geom = "point",  data = umap_m9) 
# dev.off()
# 
# rm(umap_calc)
#----------------Umap---------------------------#