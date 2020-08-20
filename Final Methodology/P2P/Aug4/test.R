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
infc<-calculate_centralities(g, include = "Information Centrality") #  takes a long time
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

#------------------------------------------------Custom tsne---------------------------------------------------#

ncen_tr<-transpose(ncentrality) #  transpose ncentrality for tsne colors
ncen_tr<-data.frame(names = c("degree", "eigenvector", "closeness", "pagerank","betweenness","hubscore","authorities",
                              "radiality","clusterrank","dmnc","lobby","leverage","subgraph","topologicalcoeff","eccentricity",
                              "geodkpath","stress","informationcent","localbridge"),ncen_tr) #y label

#-------------------------------------------------------tsne model 1 start--------------------------------------------#

set.seed(32)  
tsne_model_1_var8 = Rtsne(var8, check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.20, dims=2, max_iter = 2000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 8)
bmp("tsne_model_1_var8.bmp", width = 1920, height = 1080)
plot(tsne_model_1_var8$Y,col=1:19, asp=1)
dev.off()

#-------------------------------------------------------tsne model 1 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 2 start--------------------------------------------#

set.seed(323)  
tsne_model_2_var8 = Rtsne(var8, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 8) #lowest error
bmp("tsne_model_2_var8.bmp", width = 1920, height = 1080)
plot(tsne_model_2_var8$Y,col=1:19, asp=1)
dev.off()
#-------------------------------------------------------tsne model 2 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 3 start--------------------------------------------#

set.seed(333)  
tsne_model_3_var8 = Rtsne(var8, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 8) #lowest error
bmp("tsne_model_3_var8.bmp", width = 1920, height = 1080)
plot(tsne_model_3_var8$Y,col=1:19, asp=1)
dev.off()
#-------------------------------------------------------tsne model 3 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 4 start--------------------------------------------#

set.seed(358)  
tsne_model_4_var8 = Rtsne(var8, check_duplicates=FALSE, pca=TRUE, perplexity=43, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 8)
bmp("tsne_model_4_var8.bmp", width = 1920, height = 1080)
plot(tsne_model_4_var8$Y,col=1:19, asp=1)
dev.off()

#-------------------------------------------------------tsne model 4 end--------------------------------------------#

#---------------------------------------------kmeans start------------------------------------------------------#

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(var8, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

bmp("kmeans_kneeplot_var8.bmp", width = 1280, height = 720)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares") #kneeplot
dev.off()

#clusters using different k values var8
c1<-kmeans(var8, 2, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(var8, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(var8, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(var8, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c1<-kmeans(var8, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
#plot of clusters
p1 <- fviz_cluster(c1, geom = "point",  data = var8) + ggtitle("k = 2")
p2 <- fviz_cluster(c2, geom = "point", data = var8) + ggtitle("k = 3")
p3 <- fviz_cluster(c3, geom = "point",  data = var8) + ggtitle("k = 4")
p4 <- fviz_cluster(c4, geom = "point",  data = var8) + ggtitle("k = 5")

bmp("kmeans_pca_var8.bmp", width = 1920, height = 1280)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

c_ncentrality<-kmeans(ncentrality, 2, iter.max = 20, nstart = 25,
                      algorithm = c("Hartigan-Wong"), trace=FALSE)

bmp("kmeans_ncentrality_k2.bmp", width = 1920, height = 1280)
plot(ncentrality, col = c_ncentrality$cluster)
points(c_ncentrality$centers, col = 1:8, pch = 8)
dev.off()

c_ncentrality<-kmeans(ncentrality, 3, iter.max = 20, nstart = 25,
                      algorithm = c("Hartigan-Wong"), trace=FALSE)

bmp("kmeans_ncentrality_k3.bmp", width = 1920, height = 1280)
plot(ncentrality, col = c_ncentrality$cluster)
points(c_ncentrality$centers, col = 1:8, pch = 8)
dev.off()

rm(p1,p2,p3,p4,p5)

round(calinhara(var8,c1$cluster),digits=2) #  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(var8,c2$cluster),digits=3)
round(calinhara(var8,c3$cluster),digits=4)
round(calinhara(var8,c4$cluster),digits=5)
round(calinhara(var8,c4$cluster),digits=6)

#-------------TSNE data frames for kmeans start--------------------------------------#
#~~~~~~~~~~~~~tsne model 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d_tsne_1_var8 = as.data.frame(tsne_model_1_var8$Y) #list to dataframe var 8

#~~~~~~~~~~~~~tsne model 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d_tsne_2_var8 = as.data.frame(tsne_model_2_var8$Y) #list to dataframe var 8

#~~~~~~~~~~~~~tsne model 3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d_tsne_3_var8 = as.data.frame(tsne_model_3_var8$Y) #list to dataframe var 8

#~~~~~~~~~~~~~tsne model 4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d_tsne_4_var8 = as.data.frame(tsne_model_4_var8$Y) #list to dataframe var 8

#-------------TSNE data frames for kmeans end--------------------------------------#

#-------------kmeans on dataset and cluster onto TSNE start-------------------------------------------------#

#----------------------kmenas k=2------------------------------------#

c1<-kmeans(var8, 2, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)

plot(d_tsne_1_var8, col = c1$cluster)

#----------------------kmenas k=3------------------------------------#

c2<-kmeans(var8, 3, iter.max = 20, nstart = 25,
            algorithm = c("Hartigan-Wong"), trace=FALSE)

#Kmeans for tsne model 1

p1 <- fviz_cluster(c1, geom = "point",  data = d_tsne_1_var8) + ggtitle("k = 2 var8")


bmp("tsne_model1_kmeans_k2.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

p1 <- fviz_cluster(c2, geom = "point",  data = d_tsne_1_var8) + ggtitle("k = 3 var8")

bmp("tsne_model1_kmeans_k3.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#Kmeans for tsne model 2

p1 <- fviz_cluster(c1, geom = "point",  data = d_tsne_2_var8) + ggtitle("k = 2 var8")

bmp("tsne_model2_kmeans_k2.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

p1 <- fviz_cluster(c2, geom = "point",  data = d_tsne_2_var8) + ggtitle("k = 3 var8")

bmp("tsne_model2_kmeans_k3.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#Kmeans for tsne model 3

p1 <- fviz_cluster(c1, geom = "point",  data = d_tsne_3_var8) + ggtitle("k = 2 var8")

bmp("tsne_model3_kmeans_k2.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

p1 <- fviz_cluster(c2, geom = "point",  data = d_tsne_3_var8) + ggtitle("k = 3 var8")

bmp("tsne_model3_kmeans_k3.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

#Kmeans for tsne model 4

p1 <- fviz_cluster(c1, geom = "point",ellipse.type = "convex", data = d_tsne_4_var8) + ggtitle("k = 2 var8")

bmp("tsne_model4_kmeans_k2.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

p1 <- fviz_cluster(c2, geom = "point",  data = d_tsne_4_var8) + ggtitle("k = 3 var8")

bmp("tsne_model4_kmeans_k3.bmp", width = 1980, height = 1280)
grid.arrange(p1, ncol = 1, nrow = 1)
dev.off()

round(calinhara(tsne_model_1_var8$Y,c1$cluster),digits=2)
round(calinhara(tsne_model_1_var8$Y,c2$cluster),digits=3)

rm(c1,c2,c3,c4,c1,c6,c7,c8,c9,c2,p1,p2,p3,p4,p5)
rm(d_tsne_1_var6_degree,d_tsne_1_var5,d_tsne_1_var6, d_tsne_1_var7, d_tsne_1_var8)
rm(d_tsne_2_var6_degree,d_tsne_2_var5,d_tsne_2_var6, d_tsne_2_var7, d_tsne_2_var8)
rm(d_tsne_3_var6_degree,d_tsne_3_var5,d_tsne_3_var6, d_tsne_3_var7, d_tsne_3_var8)
rm(d_tsne_4_var6_degree,d_tsne_4_var5,d_tsne_4_var6, d_tsne_4_var7, d_tsne_4_var8)

#-----------------------kmeans on dataset and cluster onto TSNE end------------------------------------------------------#

#-------------------------------------------dbscan start-----------------------------------------------------------------#

#---------------------------------------calculating h start------------------------------------------#

bmp("dbscan_kneeplot_var8.bmp", width = 841, height = 477)
dbscan::kNNdistplot(var8, k =  2)
abline(h = 0.21, lty = 2)
dev.off()
#---------------------------------------calculating h end------------------------------------------#

#---------------------------------------plotting dbscan start--------------------------------------------------#
#ranvar<-sample_n(var8, 100000, replace = TRUE)
#set.seed(123)
res.db <- dbscan::dbscan(var8, 0.068, 2)
gc()
bmp("dbscan_var8.bmp", width = 1980, height = 1280)
fviz_cluster(res.db, var8, geom = "point")
dev.off()

#---------------------------------------plotting dbscan end--------------------------------------------------#

#-----------------------------------------------dbscan end-----------------------------------------------------------------#

#------------------GMM-----------------------------#

# dat = as.matrix(var5)
# 
# dat = center_scale(dat)
# 
# gmm = GMM(dat, 2, "maha_dist", "random_subset", 10, 10)
# 
# plot(gmm$Log_likelihood,tsne_model_4$Y, col=c(1,2,3,4,5,6,7,8))
# 
# plot(as.data.frame(tsne_model_1_var8$Y), col=factor(gmm$centroids))

xyMclust <- Mclust(as.matrix(var8), prior = priorControl(), 
                   control = emControl(), 
                   warn = mclust.options("warn"),
                   verbose = TRUE)
plot(xyMclust)

summary(xyMclust, parameters = TRUE)

plot(mclustBIC(precip), legendArgs =  list(x = "bottomleft"))
plot(mclustBIC(faithful))
plot(mclustBIC(var8))

clustgmm <- Mclust(var8, G=2)

plot(clustgmm, what=c("classification"))
plot(clustgmm, "density")
plot(clustgmm, what=c("uncertainty"))

c1$cluster<-xyMclust$classification
c1$cluster<-clustgmm$classification

fviz_cluster(c1, geom = "point",  data = var8) + ggtitle("k = 9 var8")
#------------------GMM-----------------------------#

#----------------UMAP----------------------------#
library(uwot)

umap_calc <- umap(var8, n_neighbors = 90, n_components = 2,metric = "cosine",
                  learning_rate = 1.1, scale = T, init = "spectral",
                  bandwidth = 30, negative_sample_rate = 20, n_trees = 50,
                  search_k = 2*90*50, pca_center = T, pcg_rand = T, ret_model = T,
                  ret_nn = T, n_threads = 7, verbose = getOption("verbose",TRUE),
                  grain_size = 5,  min_dist = 10, spread = 50 )

umap_var8 <- data.frame(
  UMAP1 = umap_calc$embedding[, 1],
  UMAP2 = umap_calc$embedding[, 2]
)

ggplot(umap_ranvar, aes(
  x = UMAP1, y = UMAP2,
  col = UMAP1
)) +
  geom_point()# +
#ggrepel::geom_text_repel(cex = 2.5)

plot(umap_var8$UMAP1,umap_var8$UMAP2, col=factor(umap_var8$UMAP1))
fviz_cluster(c1, geom = "point",  data = umap_var8) + ggtitle("k = 9 var8")

n<-centrality[,-13]
n$rnames<-row.names(centrality)

df_melt <- reshape2::melt(n, id.var ='rnames')

ggplot(df_melt, aes(x = factor(rnames), y = value, colour = variable)) + 
  geom_point() + xlab('nodes')
fviz_mclust(xyMclust,what = c("classification"), geom = "point",ellipse.type = "norm",palette = "jco" )
fviz_mclust(xyMclust,what = c("uncertainty"),ellipse.type = "norm",palette = "jco" )
fviz_mclust(clustgmm,what = c("classification"), geom = "point",ellipse.type = "norm",palette = "jco" )
fviz_mclust(clustgmm,what = c("uncertainty"),ellipse.type = "norm",palette = "jco" )
#----------------Umap---------------------------#

#------------Sammon's map----------------#

x<-as.matrix(var8) #var8
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "pca")
bmp("sammon_var8.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=factor(c), main="var8")
plot(sammon$Y, pch=19, col=c1$cluster, main="var8 k=2")
par(opar)
dev.off()

rm(x,opar)

#------------Sammon's map----------------#
