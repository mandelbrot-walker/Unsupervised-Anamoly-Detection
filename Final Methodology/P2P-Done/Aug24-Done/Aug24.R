#--------------index-------------------------#

#Lines       Code segment
# 54-135     Centrality calculation
# 137-167    dataset preparation
# 169-195    Correlation matrix and boxplot
# 197        PCA start 
# 199-245    PCA calculation
# 250-279    screeplot and biplot
# 283-374    contribution plot
# 378        PCA end
# 380        TSNE start 
# 388-400    tsne parameter macros
# 402-441    tsne model 1
# 443-485    tsne model 1 plots
# 491-531    tsne model 2
# 534-576    tsne model 2 plots
# 582-620    tsne model 3
# 621-665    tsne model 3 plots
# 671-709    tsne model 4
# 712-750    tsne model 4 plots
# 752        TSNE end
# 754        Spherical kmeans start
# 756-952    cluster calculation
# 955-1130   tsne m1 and kmeans
# 1132-1230  tsne m2 and kmeans
# 1232-1330  tsne m3 and kmeans
# 1332-1430  tsne m4 and kmeans
# 1432       tsne and kmeans end 
# 1434       Sammon's map
# 1436-1526  Sammon's map on models 
# 1528       Sammon's map End
# 1529       End of code
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
library(Rdimtools) #  sammon's map
library(fpc) #  calinhara()
library(clusternor) # skmeans  
#--------------libraries--------------------#



#------------------------------------------Data loader and centrality calculation start-------------------------------------# 
edges<-read.delim("p2p-Gnutella24.txt",header = TRUE, sep = "\t")

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
lbc<-local_bridging_centrality(g)

gr<-g # temporary variable gr

V(gr)$degree <- dg                               #  Degree centrality
V(gr)$eig <- eig                                 #  Eigenvector centrality
V(gr)$pagerank <- pgr                            #  Pagerank centrality
V(gr)$authorities <- auth                        #  Authority centrality
V(gr)$hubs <- hubs                               #  Hub centrality
V(gr)$betweenness <- btn                         #  Vertex betweenness centrality
V(gr)$closeness <- clsn                          #  Closeness centrality
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
round(calinhara(m1,ck2$cluster),digits=2) #  53610.28  Highest
round(calinhara(m1,ck3$cluster),digits=3) #  32490.01
round(calinhara(m1,ck4$cluster),digits=4) #  48439.63
round(calinhara(m1,ck5$cluster),digits=5) #  45391.63
round(calinhara(m1,ck6$cluster),digits=6) #  47293.93
round(calinhara(m1,ck7$cluster),digits=7) #  44038.83
round(calinhara(m1,ck8$cluster),digits=8) #  40455.15
round(calinhara(m1,ck9$cluster),digits=9) #  37139.96
round(calinhara(m1,ck10$cluster),digits=10)#  35067.73

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
round(calinhara(m2,ck2$cluster),digits=2) #  53589.66 Highest
round(calinhara(m2,ck3$cluster),digits=3) #  32538.24
round(calinhara(m2,ck4$cluster),digits=4) #  48534.09
round(calinhara(m2,ck5$cluster),digits=5) #  45516.34
round(calinhara(m2,ck6$cluster),digits=6) #  47431.1
round(calinhara(m2,ck7$cluster),digits=7) #  44164.47

#--------------clusters using different k values m3 spec

# set.seed(294)  
# m3_ldim = Rtsne(m3, check_duplicates=FALSE, pca=TRUE, perplexity=prx4, theta=th4, dims=3, max_iter = mit4,
#                         verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 5)
# 
# library(kernlab)
# 
# ck2<-specc(as.matrix(m3_ldim$Y), centers=2,kernel = "rbfdot", kpar = "automatic", 
#       nystrom.red = TRUE, nystrom.sample = dim(m3_ldim$Y)[1]/6,
#       iterations = 50, mod.sample = 0.75, na.action = na.omit)
# 
# ck2<-Skmeans(data=as.matrix(m3_ldim),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck3<-Skmeans(data=as.matrix(m3_ldim),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck4<-Skmeans(data=as.matrix(m3),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck5<-Skmeans(data=as.matrix(m3),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck6<-Skmeans(data=as.matrix(m3),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck7<-Skmeans(data=as.matrix(m3),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# 
# #  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m3,ck2$assignments),digits=4) #  4264.96
# round(calinhara(m3,ck3$cluster),digits=3) #  
# round(calinhara(m3,ck4$cluster),digits=4) #  
# round(calinhara(m3,ck5$cluster),digits=5) #  
# round(calinhara(m3,ck6$cluster),digits=6) #  
# round(calinhara(m3,ck7$cluster),digits=7) #  


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

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m4,ck2$cluster),digits=2) #  83465.41  Highest
round(calinhara(m4,ck3$cluster),digits=3) #  74451.03
round(calinhara(m4,ck4$cluster),digits=4) #  63745.7  
round(calinhara(m4,ck5$cluster),digits=5) #  65263.82
round(calinhara(m4,ck6$cluster),digits=6) #  56688.88

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
round(calinhara(m5,ck2$cluster),digits=2) #  75452.45  Highest
round(calinhara(m5,ck3$cluster),digits=3) #  63766.47
round(calinhara(m5,ck4$cluster),digits=4) #  54328.42
round(calinhara(m5,ck5$cluster),digits=5) #  52831.27
round(calinhara(m5,ck6$cluster),digits=6) #  47484.62
round(calinhara(m5,ck7$cluster),digits=7) #  42416.53
round(calinhara(m5,ck8$cluster),digits=8) #  38759.66

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
ck7<-Skmeans(data=as.matrix(m6),centers=7,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ck8<-Skmeans(data=as.matrix(m6),centers=8,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m6,ck2$cluster),digits=2) #  19080.74  Highest
round(calinhara(m6,ck3$cluster),digits=3) #  15462.16
round(calinhara(m6,ck4$cluster),digits=4) #  14105.89
round(calinhara(m6,ck5$cluster),digits=5) #  12613.1
round(calinhara(m6,ck6$cluster),digits=6) #  15736.07 
round(calinhara(m6,ck7$cluster),digits=7) #  13133.19 
round(calinhara(m6,ck8$cluster),digits=8) #  11437.27 

#--------------clusters using different k values m7 spec

# ck2<-Skmeans(data=as.matrix(m7),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck3<-Skmeans(data=as.matrix(m7),centers=3,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck4<-Skmeans(data=as.matrix(m7),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck5<-Skmeans(data=as.matrix(m7),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()
# ck6<-Skmeans(data=as.matrix(m7),centers=6,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
# gc()

#  Checking for correct no of clusters. Higher the index value better the cluster
# round(calinhara(m7,ck2$cluster),digits=2) #  
# round(calinhara(m7,ck3$cluster),digits=3) #  
# round(calinhara(m7,ck4$cluster),digits=4) #  
# round(calinhara(m7,ck5$cluster),digits=5) #  
# round(calinhara(m7,ck6$cluster),digits=6) #  

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

#  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(m8,ck2$cluster),digits=2) #  57394.79  
round(calinhara(m8,ck3$cluster),digits=3) #  63905.46
round(calinhara(m8,ck4$cluster),digits=4) #  70264.44 Highest
round(calinhara(m8,ck5$cluster),digits=5) #  64950.7
round(calinhara(m8,ck6$cluster),digits=6) #  57677.19
round(calinhara(m8,ck7$cluster),digits=7) #  13419.63
round(calinhara(m8,ck8$cluster),digits=8) #  11770.35

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
round(calinhara(m9,ck2$cluster),digits=2) #  53365.65  Highest
round(calinhara(m9,ck3$cluster),digits=3) #  42551.87
round(calinhara(m9,ck4$cluster),digits=4) #  50662.36
round(calinhara(m9,ck5$cluster),digits=5) #  40831.18
round(calinhara(m9,ck6$cluster),digits=6) #  35619.3

rm(ck2,ck3,ck4,ck5,ck6,ck7,ck8,ck9)
#-------------kmeans on dataset and cluster onto TSNE start-------------------------------------------------#

#-----------------------------------------Kmeans for tsne model 1------------------------------------------#

ckm1<-Skmeans(data=as.matrix(m1),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm2<-Skmeans(data=as.matrix(m2),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
#ckm3<-Skmeans(data=as.matrix(m3),centers=5,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm4<-Skmeans(data=as.matrix(m4),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm5<-Skmeans(data=as.matrix(m5),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm6<-Skmeans(data=as.matrix(m6),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
#ckm7<-Skmeans(data=as.matrix(m7),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm8<-Skmeans(data=as.matrix(m8),centers=4,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()
ckm9<-Skmeans(data=as.matrix(m9),centers=2,iter.max = 25,nthread = 5,init = c("kmeanspp"),tolerance = 0.0005)
gc()  # garbage collection is used for stack imbalance warning. run gc() more than once if the warning persists

ck1<-kmeans(m1, 3, iter.max = 20, nstart = 25,
            algorithm = c("Hartigan-Wong"), trace=FALSE) #  dummy kmeans
ck2<-ck1
#ck3<-ck1
ck4<-ck1
ck5<-ck1
ck6<-ck1
#ck7<-ck1
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
# bmp("tsne_model1_m3_kmeans.bmp", width = 1980, height = 1280)
# plot(as.data.frame(tsne_model_1_m3$Y), col = ckm3$cluster)
# dev.off()
# 
# ck3$cluster<-ckm3$cluster
# ck3$centers<-ckm3$centers
# ck3$size<-ckm3$size
# ck3$iter<-ckm3$iters
# 
# p1 <- fviz_cluster(ck3, geom = "point",  data = as.data.frame(tsne_model_1_m3$Y)) 
# 
# bmp("tsne_model1_m3_kmeans_ch.bmp", width = 1980, height = 1280)
# grid.arrange(p1, ncol = 1, nrow = 1)
# dev.off()

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
# bmp("tsne_model1_m7_kmeans.bmp", width = 1980, height = 1280)
# plot(as.data.frame(tsne_model_1_m7$Y), col = ckm7$cluster)
# dev.off()
# 
# ck7$cluster<-ckm7$cluster
# ck7$centers<-ckm7$centers
# ck7$size<-ckm7$size
# ck7$iter<-ckm7$iters
# 
# p1 <- fviz_cluster(ck7, geom = "point",  data = as.data.frame(tsne_model_1_m7$Y)) 
# 
# bmp("tsne_model1_m7_kmeans_ch.bmp", width = 1980, height = 1280)
# grid.arrange(p1, ncol = 1, nrow = 1)
# dev.off()

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
# bmp("tsne_model2_m3_kmeans.bmp", width = 1980, height = 1280)
# plot(as.data.frame(tsne_model_2_m3$Y), col = ckm3$cluster)
# dev.off()
# 
# p1 <- fviz_cluster(ck3, geom = "point",  data = as.data.frame(tsne_model_2_m3$Y)) 
# 
# bmp("tsne_model2_m3_kmeans_ch.bmp", width = 1980, height = 1280)
# grid.arrange(p1, ncol = 1, nrow = 1)
# dev.off()

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
# bmp("tsne_model2_m7_kmeans.bmp", width = 1980, height = 1280)
# plot(as.data.frame(tsne_model_2_m7$Y), col = ckm7$cluster)
# dev.off()
# 
# p1 <- fviz_cluster(ck7, geom = "point",  data = as.data.frame(tsne_model_2_m7$Y)) 
# 
# bmp("tsne_model2_m7_kmeans_ch.bmp", width = 1980, height = 1280)
# grid.arrange(p1, ncol = 1, nrow = 1)
# dev.off()

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
# bmp("tsne_model3_m3_kmeans.bmp", width = 1980, height = 1280)
# plot(as.data.frame(tsne_model_3_m3$Y), col = ckm3$cluster)
# dev.off()
# 
# p1 <- fviz_cluster(ck3, geom = "point",  data = as.data.frame(tsne_model_3_m3$Y)) 
# 
# bmp("tsne_model3_m3_kmeans_ch.bmp", width = 1980, height = 1280)
# grid.arrange(p1, ncol = 1, nrow = 1)
# dev.off()

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
# bmp("tsne_model3_m7_kmeans.bmp", width = 1980, height = 1280)
# plot(as.data.frame(tsne_model_3_m7$Y), col = ckm7$cluster)
# dev.off()
# 
# p1 <- fviz_cluster(ck7, geom = "point",  data = as.data.frame(tsne_model_3_m7$Y)) 
# 
# bmp("tsne_model3_m7_kmeans_ch.bmp", width = 1980, height = 1280)
# grid.arrange(p1, ncol = 1, nrow = 1)
# dev.off()

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
# bmp("tsne_model4_m3_kmeans.bmp", width = 1980, height = 1280)
# plot(as.data.frame(tsne_model_4_m3$Y), col = ckm3$cluster)
# dev.off()
# 
# p1 <- fviz_cluster(ck3, geom = "point",  data = as.data.frame(tsne_model_4_m3$Y)) 
# 
# bmp("tsne_model4_m3_kmeans_ch.bmp", width = 1980, height = 1280)
# grid.arrange(p1, ncol = 1, nrow = 1)
# dev.off()

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
# bmp("tsne_model4_m7_kmeans.bmp", width = 1980, height = 1280)
# plot(as.data.frame(tsne_model_4_m7$Y), col = ckm7$cluster)
# dev.off()
# 
# p1 <- fviz_cluster(ck7, geom = "point",  data = as.data.frame(tsne_model_4_m7$Y)) 
# 
# bmp("tsne_model4_m7_kmeans_ch.bmp", width = 1980, height = 1280)
# grid.arrange(p1, ncol = 1, nrow = 1)
# dev.off()

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

set.seed(729)
x<-as.matrix(m1) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m1.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=factor(ncen_tr$names), main="m1")
plot(sammon$Y, pch=19, col=ckm1$cluster, main="m1 kmeans")
par(opar)
dev.off()

set.seed(728)
x<-as.matrix(m2) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m2.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:12, main="m2")
plot(sammon$Y, pch=19, col=ckm2$cluster, main="m2 kmeans")
par(opar)
dev.off()

set.seed(727)
x<-as.matrix(m3) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m3.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:7, main="m3")
plot(sammon$Y, pch=19, col=ckm3$cluster, main="m3 kmeans")
par(opar)
dev.off()

set.seed(726)
x<-as.matrix(m4) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m4.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:6, main="m4")
plot(sammon$Y, pch=19, col=ckm4$cluster, main="m4 kmeans")
par(opar)
dev.off()

set.seed(725)
x<-as.matrix(m5) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m5.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:9, main="m5")
plot(sammon$Y, pch=19, col=ckm5$cluster, main="m5 kmeans")
par(opar)
dev.off()

set.seed(724)
x<-as.matrix(m6) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m6.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:9, main="m6")
plot(sammon$Y, pch=19, col=ckm6$cluster, main="m6 kmeans")
par(opar)
dev.off()

set.seed(723)
x<-as.matrix(m7) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m7.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:6, main="m7")
plot(sammon$Y, pch=19, col=ckm7$cluster, main="m7 kmeans")
par(opar)
dev.off()

set.seed(722)
x<-as.matrix(m8) 
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "random")
bmp("sammon_m8.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=1:6, main="m8")
plot(sammon$Y, pch=19, col=ckm8$cluster, main="m8 kmeans")
par(opar)
dev.off()

set.seed(721)
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

tsne_m1_m1<-as.data.frame(tsne_model_1_m1$Y)

tsne_m1_m1$V1<-(tsne_m1_m1$V1-min(tsne_m1_m1$V1))/(max(tsne_m1_m1$V1)-min(tsne_m1_m1$V1))
tsne_m1_m1$V2<-(tsne_m1_m1$V2-min(tsne_m1_m1$V2))/(max(tsne_m1_m1$V2)-min(tsne_m1_m1$V2))

library(MASS)
den3d <- kde2d(tsne_m1_m1$V1, tsne_m1_m1$V2)

library(plotly)
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()

hist(tsne_m1_m1$V1, col=rgb(1,0,0,0.5), main = "axis 1 and 2")
hist(tsne_m1_m1$V2, col=rgb(0,0,1,0.5), add=T)
box()

# library(Spectrum)  # spectral clustering
# 
# ck2 <- Spectrum(m3, method = 2, silent = FALSE, showres = TRUE,
#                 diffusion = TRUE, kerneltype = c("density", "stsc"), maxk = 5,
#                 NN = 25, NN2 = 35, showpca = FALSE, frac = 2, thresh = 7,
#                 fontsize = 18, dotsize = 3, tunekernel = T,
#                 clusteralg = "GMM", FASP = FALSE, FASPk = NULL, fixk = NULL,
#                 krangemax = 10, runrange = FALSE, diffusion_iters = 4,
#                 KNNs_p = 50, missing = FALSE)

# library(proxy)
# 
# smat<-dist(as.matrix(m3), y = NULL, method = "eJaccard",diag = FALSE, upper = FALSE,
#            pairwise = F, by_rows = F, convert_similarities = T,
#            auto_convert_data_frames = F)
# 
# smat2<-as.matrix(smat)

