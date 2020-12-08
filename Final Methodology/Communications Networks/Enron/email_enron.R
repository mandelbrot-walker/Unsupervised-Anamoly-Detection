library(igraph) #  centralities
library(centiserve) #  centralities
library(factoextra)  # fviz_cluster()
library(Rtsne) #  TSNE
library(Plasmidprofiler) #  normalize()
library(data.table) #  transpose()
library(corrplot) #  correlation plot
library(gridExtra) #  multiplot
library(CINNA) #  centralities
library(Rdimtools) #  sammon's map
library(fpc) #  calinhara()
library(clusternor) # skmeans 

#------------------------------------------Data loader and centrality calculation start-------------------------------------# 
edges<-read.delim("Email-Enron.txt",header = TRUE, sep = "\t")


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
round(calinhara(m1,ck2$cluster),digits=2) #  12936.01  
round(calinhara(m1,ck3$cluster),digits=3) #  34570.29
round(calinhara(m1,ck4$cluster),digits=4) #  44222.15  
round(calinhara(m1,ck5$cluster),digits=5) #  43585.78
round(calinhara(m1,ck6$cluster),digits=6) #  42654.56
round(calinhara(m1,ck7$cluster),digits=7) #  42110.6
round(calinhara(m1,ck8$cluster),digits=8) #  43480.34
round(calinhara(m1,ck9$cluster),digits=9) #  44695.7  Highest
round(calinhara(m1,ck10$cluster),digits=10) #  43537.76

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
ck8<-Skmeans(data=as.matrix(m2),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck9<-Skmeans(data=as.matrix(m2),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck10<-Skmeans(data=as.matrix(m2),centers=10,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m2,ck2$cluster),digits=2) #  12997.19 
round(calinhara(m2,ck3$cluster),digits=3) #  34960.27
round(calinhara(m2,ck4$cluster),digits=4) #  44978.74  
round(calinhara(m2,ck5$cluster),digits=5) #  44174.32
round(calinhara(m2,ck6$cluster),digits=6) #  43284.42
round(calinhara(m2,ck7$cluster),digits=7) #  42648.32
round(calinhara(m2,ck8$cluster),digits=8) #  44910.6
round(calinhara(m2,ck9$cluster),digits=9) #  45455.04  Highest
round(calinhara(m2,ck10$cluster),digits=10) #  44112.46

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
ck8<-Skmeans(data=as.matrix(m3),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck9<-Skmeans(data=as.matrix(m3),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck10<-Skmeans(data=as.matrix(m3),centers=10,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck11<-Skmeans(data=as.matrix(m3),centers=11,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck12<-Skmeans(data=as.matrix(m3),centers=12,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck13<-Skmeans(data=as.matrix(m3),centers=13,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck14<-Skmeans(data=as.matrix(m3),centers=14,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck15<-Skmeans(data=as.matrix(m3),centers=15,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()


#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m3,ck2$cluster),digits=2) #  31899.08
round(calinhara(m3,ck3$cluster),digits=3) #  36938.21
round(calinhara(m3,ck4$cluster),digits=4) #  37607.78
round(calinhara(m3,ck5$cluster),digits=5) #  59122.76  Highest  
round(calinhara(m3,ck6$cluster),digits=6) #  52820.22
round(calinhara(m3,ck7$cluster),digits=7) #  48335.24
round(calinhara(m3,ck8$cluster),digits=8) #  41940.8
round(calinhara(m3,ck9$cluster),digits=9) #  37155.67
round(calinhara(m3,ck10$cluster),digits=10) #  33434.91
round(calinhara(m3,ck11$cluster),digits=11) #  31363.69
round(calinhara(m3,ck12$cluster),digits=12) #  31706.2
round(calinhara(m3,ck13$cluster),digits=13) #  29048.65  
round(calinhara(m3,ck14$cluster),digits=15) #  27777.12
round(calinhara(m3,ck15$cluster),digits=15) #  30789.98

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
ck9<-Skmeans(data=as.matrix(m4),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck10<-Skmeans(data=as.matrix(m4),centers=10,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck11<-Skmeans(data=as.matrix(m4),centers=11,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m4,ck2$cluster),digits=2) #  18016.54
round(calinhara(m4,ck3$cluster),digits=3) #  30123.89
round(calinhara(m4,ck4$cluster),digits=4) #  46641.11  
round(calinhara(m4,ck5$cluster),digits=5) #  45708.24
round(calinhara(m4,ck6$cluster),digits=6) #  44403.59
round(calinhara(m4,ck7$cluster),digits=7) #  54236.1  Highest
round(calinhara(m4,ck8$cluster),digits=8) #  51224.19
round(calinhara(m4,ck9$cluster),digits=9) #  47669.8
round(calinhara(m4,ck10$cluster),digits=10) #  45405.62
round(calinhara(m4,ck11$cluster),digits=11) #  42463.55

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
ck7<-Skmeans(data=as.matrix(m5),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck8<-Skmeans(data=as.matrix(m5),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck9<-Skmeans(data=as.matrix(m5),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck10<-Skmeans(data=as.matrix(m5),centers=10,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m5,ck2$cluster),digits=2) #  13265.42  
round(calinhara(m5,ck3$cluster),digits=3) #  36626.3
round(calinhara(m5,ck4$cluster),digits=4) #  48404.81  Highest
round(calinhara(m5,ck5$cluster),digits=5) #  46690.67
round(calinhara(m5,ck6$cluster),digits=6) #  45826.62
round(calinhara(m5,ck7$cluster),digits=7) #  45826.62
round(calinhara(m5,ck7$cluster),digits=8) #  45826.62
round(calinhara(m5,ck7$cluster),digits=9) #  45826.62
round(calinhara(m5,ck7$cluster),digits=10) #  45826.62

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
ck7<-Skmeans(data=as.matrix(m6),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m6,ck2$cluster),digits=2) #  26209.5  Highest
round(calinhara(m6,ck3$cluster),digits=3) #  14362.43
round(calinhara(m6,ck4$cluster),digits=4) #  9886.359
round(calinhara(m6,ck5$cluster),digits=5) #  7736.287
round(calinhara(m6,ck6$cluster),digits=6) #  6307.078 
round(calinhara(m6,ck7$cluster),digits=7) #  6307.078 

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
ck7<-Skmeans(data=as.matrix(m7),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m7,ck2$cluster),digits=2) #  45689.63  Highest
round(calinhara(m7,ck3$cluster),digits=3) #  33279
round(calinhara(m7,ck4$cluster),digits=4) #  26569.88  
round(calinhara(m7,ck5$cluster),digits=5) #  23089.9
round(calinhara(m7,ck6$cluster),digits=6) #  18624.55
round(calinhara(m7,ck7$cluster),digits=7) #  18553.56

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
ck7<-Skmeans(data=as.matrix(m8),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck8<-Skmeans(data=as.matrix(m8),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck9<-Skmeans(data=as.matrix(m8),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck10<-Skmeans(data=as.matrix(m8),centers=10,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck11<-Skmeans(data=as.matrix(m8),centers=11,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck12<-Skmeans(data=as.matrix(m8),centers=12,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck13<-Skmeans(data=as.matrix(m8),centers=13,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck14<-Skmeans(data=as.matrix(m8),centers=14,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck15<-Skmeans(data=as.matrix(m8),centers=15,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m8,ck2$cluster),digits=2) #  17915.95  Highest  
round(calinhara(m8,ck3$cluster),digits=3) #  10564.62
round(calinhara(m8,ck4$cluster),digits=4) #  9003.982
round(calinhara(m8,ck5$cluster),digits=5) #  7950.469
round(calinhara(m8,ck6$cluster),digits=6) #  7501.928
round(calinhara(m8,ck7$cluster),digits=7) #  8069.734
round(calinhara(m8,ck8$cluster),digits=8) #  7145.015
round(calinhara(m8,ck9$cluster),digits=9) #  6575.735  
round(calinhara(m8,ck10$cluster),digits=10) #  5896.039
round(calinhara(m8,ck11$cluster),digits=11) #  5326.471
round(calinhara(m8,ck12$cluster),digits=12) #  4858.539
round(calinhara(m8,ck13$cluster),digits=13) #  4453.02
round(calinhara(m8,ck14$cluster),digits=14) #  4117.116
round(calinhara(m8,ck15$cluster),digits=15) #  3832.258

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
ck7<-Skmeans(data=as.matrix(m8),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck8<-Skmeans(data=as.matrix(m4),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck9<-Skmeans(data=as.matrix(m4),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck10<-Skmeans(data=as.matrix(m4),centers=10,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck11<-Skmeans(data=as.matrix(m4),centers=11,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m9,ck2$cluster),digits=2) #  2627.33  
round(calinhara(m9,ck3$cluster),digits=3) #  13716.85
round(calinhara(m9,ck4$cluster),digits=4) #  10191.85
round(calinhara(m9,ck5$cluster),digits=5) #  9063.424
round(calinhara(m9,ck6$cluster),digits=6) #  7812.304
round(calinhara(m9,ck7$cluster),digits=7) #  8118.445  
round(calinhara(m9,ck8$cluster),digits=8) #  50158.82  Highest
round(calinhara(m9,ck9$cluster),digits=9) #  46688.42
round(calinhara(m9,ck10$cluster),digits=10) #  42286.55
round(calinhara(m9,ck11$cluster),digits=11) #  40026.56



rm(ck2,ck3,ck4,ck5,ck6,ck7,ck8,ck9,ck10,ck11,ck12,ck13,ck14,ck15)
#-------------kmeans on dataset and cluster onto TSNE start-------------------------------------------------#

#-----------------------------------------Kmeans for tsne model 1------------------------------------------#

ckm1<-Skmeans(data=as.matrix(m1),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm2<-Skmeans(data=as.matrix(m2),centers=9,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm3<-Skmeans(data=as.matrix(m3),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm4<-Skmeans(data=as.matrix(m4),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm5<-Skmeans(data=as.matrix(m5),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm6<-Skmeans(data=as.matrix(m6),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm7<-Skmeans(data=as.matrix(m7),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm8<-Skmeans(data=as.matrix(m8),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm9<-Skmeans(data=as.matrix(m9),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
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

