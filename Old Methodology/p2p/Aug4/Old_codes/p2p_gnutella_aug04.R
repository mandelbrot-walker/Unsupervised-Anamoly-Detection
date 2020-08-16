library(igraph)
library(centiserve)
library(tidyverse)
library(factoextra)
library(fastnet)
library(Rtsne)
library(Plasmidprofiler)
library(MASS)
library(data.table)
library(corrplot)
library(tibble)
library(caret)
library(plyr)
library(gridExtra) 
library(resample)
library(earth)
library(CINNA)
library(ClusterR)
library(mclust)
library(kohonen)
library(scatterplot3d)
#------------------------------------------Data loader and centrality calculation start-------------------------------------# 
edges<-read.delim("p2p-Gnutella04.txt",header = TRUE, sep = "\t")

g<-graph.data.frame(edges) #graph data frame for igraph

cent<-proper_centralities(g)

calculate_centralities(g, include = cent[1:50])%>%
  pca_centralities(scale.unit = TRUE, ncp = 50) # takes indefinite time

dg<-degree(g) # Calculation of Degre centrality
write.csv(dg, "dg_p2p_Gnutella04.csv")

btn<-betweenness(g) # Calculation of Betweeness centrality
write.csv(btn, "btn_p2p_Gnutella04.csv")

eig<-evcent(g)$vector # Calculation of Eigenvector centrality
write.csv(eig, "eig_p2p_Gnutella04.csv")

clsn<-closeness(g) # Calculation of Closeness centrality
write.csv(clsn, "clsn_p2p_Gnutella04.csv")

pgr<-page_rank(g)$vector #Calculation of Page Rank centrality
write.csv(pgr, "pgr_p2p_Gnutella04.csv")

katz<-katzcent(g) # Error in alpha >= maxEigenvalue : invalid comparison with complex values

crsc<-crossclique(g) # Calculation of Cross-Clique centrality
write.csv(crsc, "crsc_p2p_Gnutella04.csv")

frm<-closeness.freeman(g) # Not calculatable as graphis not strongly connected

edge_connectivity(g) # Outputs 0

#--------Plotting takes too long--------------------------------------#
lay <- layout.fruchterman.reingold(g) # For plotting graph
lay2 <- layout_with_fr(g) # For plotting graph

plot.igraph(g, layout=lay, 
            vertex.size=btn*0.25,    # Rescaled by multiplying by 25
            main="Betweenness")
plot(g, layout = lay, vertex.label = NA)
plot(g, layout = lay2, vertex.label = NA)
#----------Plotting section end--------------------------------------#


gr<-g # temporary variable gr

V(gr)$degree <- dg                               # Degree centrality
V(gr)$eig <- eig                                 # Eigenvector centrality
V(gr)$closeness <- clsn                          # Closeness centrality
V(gr)$pagerank <- pgr                            # Pagerank centrality
V(gr)$betweenness <- btn                         # Vertex betweenness centrality
V(gr)$crossclique <- crsc                        # Crossclique centrality
V(gr)$hubs <- hub.score(g)$vector                # "Hub" centrality
V(gr)$authorities <- authority.score(g)$vector   # "Authority" centrality

centrality <- data.frame(row.names   = V(gr)$name,
                         degree      = V(gr)$degree,
                         eigenvector = V(gr)$eig,
                         closeness   = V(gr)$closeness,
                         pagerank    = V(gr)$pagerank,
                         crossclique = V(gr)$crossclique,
                         betweenness = V(gr)$betweenness,
                         hubscore    = V(gr)$hubs,
                         authorities = V(gr)$authorities
)

centrality <- centrality[order(row.names(centrality)),] # These values are not normalized


head(centrality) #check centrality variables

ndegree      = normalize(dg) 
neigenvector = normalize(eig) 
ncloseness   = normalize(clsn)
npagerank    = normalize(pgr)
ncrossclique = normalize(crsc)
nbetweenness = normalize(btn)
nhubscore    = normalize(V(gr)$hubs)
nauthorities = normalize(V(gr)$authorities)

ncentrality <- data.frame(degree      = ndegree,
                          eigenvector = neigenvector,
                          closeness   = ncloseness,
                          pagerank    = npagerank,
                          crossclique = ncrossclique,
                          betweenness = nbetweenness
) #normalized values 6 variables

ncentrality2 <- data.frame(degree      = ndegree,
                           eigenvector = neigenvector,
                           closeness   = ncloseness,
                           pagerank    = npagerank,
                           crossclique = ncrossclique,
                           betweenness = nbetweenness,
                           hubscore    = nhubscore,
                           authorities = nauthorities 
)#normalized values 8 variables

#------------------------------------------Data loader and centrality calculation End-------------------------------------#

#------------------------------------------Write data start-------------------------------------#
write.csv(centrality,"Centrality_p2p_Gnutella04.csv")
write.csv(ncentrality,"NCentrality_p2p_Gnutella04.csv")
write.csv(ncentrality2,"NCentrality_p2p_Gnutella04.csv")

#------------------------------------------Write data end-------------------------------------#

#------------------------------------------Boxplot start-------------------------------------#
M <- cor(ncentrality2)
corrplot(M, method = "circle") # correlation matrix 

boxplot(ndegree, neigenvector, ncloseness, npagerank,ncrossclique,nbetweenness,nhubscore,nauthorities,
        main = "Multiple boxplots for comparision",
        at = c(1,2,3,4,5,6,7,8),
        names = c("ndegree", "neigenvector", "ncloseness", "npagerank","ncrossclique","nbetweenness","nhubscore","nauthorities"),
        las = 2,
        col = c("orange","red"),
        border = "brown",
        horizontal = TRUE,
        notch = FALSE
) #multiple boxplot 
#------------------------------------------Boxplot end-------------------------------------#

#------------------------------------------PCA start-------------------------------------#

#principal component analysis

pcacen<-ncentrality#keeping ncentrality as it is
pcacen<-within(pcacen, rm(degree)) #without degree total 5 variables

pcacen2<-ncentrality2#keeping ncentrality2 as it is 
pcacen2<-within(pcacen2, rm(degree)) #without degree total 7 variables

res.pca<-prcomp(scale(ncentrality),center=TRUE) #with degree 6 variables
res.pca2<-prcomp(scale(pcacen),center=TRUE) #without degree 5 variables 
res.pca3<-prcomp(scale(ncentrality2),center=TRUE) #with degree 8 variables
res.pca4<-prcomp(scale(pcacen2),center=TRUE) #without degree 7 variables 

print(res.pca) #show pca values
print(res.pca2)#show pca values

fviz_eig(res.pca) #Scree plot 6 var
fviz_eig(res.pca2)#Scree plot 5 var
fviz_eig(res.pca3)#Scree plot 8 var
fviz_eig(res.pca4)#Scree plot 7 var

biplot(res.pca,scale=0, cex=.7)#biplot 6 var
biplot(res.pca2,scale=0, cex=.7)#biplot 5 vat
biplot(res.pca3,scale=0, cex=.7)#biplot 8 var
biplot(res.pca4,scale=0, cex=.7)#biplot 7 var

#hist(dg, breaks=5, main="With breaks=30") #histogram did not work

# Color by the quality of representation !--takes too long--!
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)#with degree 6 variables

fviz_pca_ind(res.pca2,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)#without degree 5 variables 

# Color by contributions to the PC

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)#with degree 6 variables

fviz_pca_var(res.pca2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)#without degree 5 variables 

fviz_pca_var(res.pca3,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)#with degree 8 variables

fviz_pca_var(res.pca4,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)#without degree 7 variables 

# Biplot with colored variables and no overlap

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)#with degree 6 variables


fviz_pca_biplot(res.pca2, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)#without degree 5 variables 

res.pca22<-prcomp(ncen_tr[,-1],  scale = TRUE,center=TRUE) #extra pca to plot with ellipse

fviz_pca_biplot(res.pca22, 
                # Individuals
                geom.ind = "point",
                fill.ind = ncen_tr$names, col.ind = "white",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdBu"
)+  labs(fill = "Variables", color = "Contrib", alpha = "Contrib") # ellipse does work 


eig.val <- get_eigenvalue(res.pca)
eig.val #shows eigenvalues 

# Results for Variables
res.var <- get_pca_var(res.pca4)#without degree 7 variables 
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca4)#without degree 7 variables 
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation
#------------------------------------------PCA end-------------------------------------#

#------------------------------------------------Custom tsne---------------------------------------------------#

ncen_tr<-transpose(ncentrality2) #transpose ncentrality for tsne
ncen_tr<-data.frame(names = c('degree','eigenvector','closeness','pagerank','crossclique','betweeness','hubscore','authorities'),ncen_tr) #y label
colors<-c("red","purple","blue","green","black","yellow","orange","magenta")
#tsne model 4
set.seed(32)  
tsne_model_4 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.20, dims=2, max_iter = 2000,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
plot(tsne_model_4$Y,col=colors, asp=1)

#tsne model 5
set.seed(323)  
tsne_model_5 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
plot(tsne_model_5$Y,col=ncen_tr$names, asp=1)

#tsne model 6
set.seed(333)  
tsne_model_6 = Rtsne(pcacen2, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
plot(tsne_model_6$Y,col=ncen_tr$names, asp=1)


#tsne model 7
set.seed(358)  
tsne_model_7 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=43, theta=0.10, dims=2, max_iter = 5000,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
plot(tsne_model_7$Y,col=ncen_tr$names, asp=1)

d_tsne_1 = as.data.frame(tsne_model_5$Y) #list to dataframe
d_tsne_1_original=d_tsne_1 #keeping the original

set.seed(358)  
tsne_model_t = Rtsne(as.matrix(ncentrality2), check_duplicates=FALSE, pca=TRUE, perplexity=43, theta=0.10, dims=3, max_iter = 1000,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
plot(tsne_model_t$Y,col=ncen_tr$names, asp=1)

scatterplot3d(x=tsne_model_t$Y[,1],y=tsne_model_t$Y[,2],z=tsne_model_t$Y[,3],
              color = c("blue"), asp=1)

library(rgl)
plot3d(x=tsne_model_t$Y[,1],y=tsne_model_t$Y[,2],z=tsne_model_t$Y[,3],
       col=c("red","purple","blue","green","black","yellow","orange","magenta"),
       type="s",radius=0.5)

set.seed(358)  
tsne_model_t5 = Rtsne(as.matrix(ncentrality2), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=3, max_iter = 1000,
                      verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6)

plot3d(x=tsne_model_t5$Y[,1],y=tsne_model_t5$Y[,2],z=tsne_model_t5$Y[,3],
       col=c("red","purple","blue","green","black","yellow","orange","magenta"),
       type="s",radius=0.5)

set.seed(358)  
tsne_model_t5_1 = Rtsne(ncentrality2, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                      verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6)
plot(tsne_model_t5_1$Y,col=ncen_tr$names, asp=1)
#---------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------kmeans starts from here------------------------------------------------------#

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(d_tsne_1, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares") #kneeplot

#clusters using different k values
c1<-kmeans(pcacen2, 2, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(pcacen2, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(pcacen2, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(pcacen2, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)

c4<-kmeans(pcacen, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)

plot(pcacen, col = c4$cluster)
points(c4$centers, col = 1:6, pch = 8)

plot(tsne_model_6$Y, col = c4$cluster)
points(c4$centers, col = 1:6, pch = 8)

#plot of clusters
p1 <- fviz_cluster(c1, geom = "point",  data = pcacen2) + ggtitle("k = 3")
p2 <- fviz_cluster(c2, geom = "point", data = pcacen2) + ggtitle("k = 4")
p3 <- fviz_cluster(c3, geom = "point",  data = pcacen2) + ggtitle("k = 5")
p4 <- fviz_cluster(c4, geom = "point",  data = pcacen2) + ggtitle("k = 6")

#grid arrangement
grid.arrange(p1, p2, p3, p4, nrow = 2)

#clustering for tsne model
c5<-kmeans(d_tsne_1, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c6<-kmeans(d_tsne_1, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)

p5 <- fviz_cluster(c5, geom = "point",  data = d_tsne_1) + ggtitle("k = 4")
p6 <- fviz_cluster(c6, geom = "point",  data = d_tsne_1) + ggtitle("k = 5")

grid.arrange(p5, p6, nrow = 1)

#---------------------------------------------kmeans ends from here------------------------------------------------------#

#-----------------------------------------------dbscan-----------------------------------------------------------------#

dbscan::kNNdistplot(pcacen2, k =  2)
abline(h = 0.04, lty = 2)

set.seed(123)
res.db <- dbscan::dbscan(pcacen2, 0.04, 2)
fviz_cluster(res.db, pcacen2, geom = "point")

#-------------------test---------------------------#
## keeping original data
d_tsne_1_original=as.data.frame(tsne_model_7$Y)
d_tsne_1 =tsne_model_7
## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1), 3)
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

## setting 3 clusters as output
d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=3))

plot_cluster=function(data, var_cluster, palette)
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.25) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
}


plot_k=plot_cluster(d_tsne_1_original, "cl_kmeans", "Accent")
plot_h=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")

## and finally: putting the plots side by side with gridExtra lib...
library(gridExtra)
grid.arrange(plot_k, plot_h,  ncol=2)
#--------------------test--------------------------#

#-------------------test2---------------------------#
d <- dist(tsne_model_7$Y, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

c1<-kmeans(ncentrality2[,-1], 2, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(ncentrality2, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(ncentrality2, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(ncentrality2, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)

c4<-kmeans(pcacen, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)

plot(ncentrality2[-1], col = c4$cluster)
points(c4$centers, col = 1:6, pch = 8)

plot(tsne_model_6$Y, col = c4$cluster)
points(c4$centers, col = 1:6, pch = 8)

#plot of clusters
p1 <- fviz_cluster(c1, geom = "point",  data = ncentrality2) + ggtitle("k = 3")
p2 <- fviz_cluster(c2, geom = "point", data = ncentrality2) + ggtitle("k = 4")
p3 <- fviz_cluster(c3, geom = "point",  data = ncentrality2) + ggtitle("k = 5")
p4 <- fviz_cluster(c4, geom = "point",  data = ncentrality2) + ggtitle("k = 6")

#grid arrangement
grid.arrange(p1, p2, p3, p4, nrow = 2)
#------------------test2-----------------------------#

#------------------test3-----------------------------#

dat = as.matrix(ncentrality2)

dat = center_scale(dat)

gmm = GMM(dat, 2, "maha_dist", "random_subset", 10, 10)

plot(gmm$Log_likelihood,tsne_model_7$Y, col=c(1,2,3,4,5,6,7,8))

xyMclust <- Mclust(ncentrality2)
plot(xyMclust)
summary(xyMclust, parameters = TRUE)

clustgmm <- Mclust(ncentrality2,5)

plot(clustgmm, what=c("classification"))
plot(clustgmm, "density")
plot(clustgmm, what=c("BIC"))
#------------------test3-----------------------------#

#------------------SOM-----------------------------#

data_train_matrix <- as.matrix(scale(ncentrality2))

# Create the SOM Grid - you generally have to specify the size of the 
# training grid prior to training the SOM. Hexagonal and Circular 
# topologies are possible
som_grid <- somgrid(xdim = 20, ydim=20, topo="hexagonal")

# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(data_train_matrix, 
                       grid = som.grid, 
                       rlen = 100,
                       alpha = c(0.05,0.01),
                       keep.data = FALSE,
                       n.hood = "circular",
                       toroidal = T)
plot(som_model, type="changes")
plot(som_model, type="count")
plot(som_model, type="dist.neighbours")
plot(som_model, type="codes")

plot(som_model, type = "property", property = som_model$codes[,4], main=names(som_model$data)[4], palette.name=coolBlueHotRed)

var <- 2 #define the variable to plot 
var_unscaled <- aggregate(as.numeric(data_train[,var]), by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,2] 
plot(som_model, type = "property", property=var_unscaled, main=names(data_train)[var], palette.name=coolBlueHotRed)

mydata <- som_model$codes 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
for (i in 2:15) {
  wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
}
plot(wss)

## use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(som_model$codes)), 6)
# plot these results:
plot(som_model, type="mapping", bgcol = pretty_palette[som_cluster], main = "Clusters") 
add.cluster.boundaries(som_model, som_cluster)

#------------------SOM-----------------------------#
