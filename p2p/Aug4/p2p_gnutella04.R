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

edges<-read.delim("p2p-Gnutella04.txt",header = TRUE, sep = "\t")

g<-graph.data.frame(edges) #graph data frame for igraph

cent<-proper_centralities(g)

calculate_centralities(g, include = cent[1:50])%>%
  pca_centralities(scale.unit = TRUE, ncp = 50)

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

lay <- layout.fruchterman.reingold(g) # For plotting graph
lay2 <- layout_with_fr(g) # For plotting graph

plot.igraph(g, layout=lay, 
            vertex.size=btn*0.25,    # Rescaled by multiplying by 25
            main="Betweenness")
plot(g, layout = lay, vertex.label = NA)
plot(g, layout = lay2, vertex.label = NA)

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
write.csv(centrality,"Centrality_p2p_Gnutella04.csv")

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

write.csv(ncentrality,"NCentrality_p2p_Gnutella04.csv")
write.csv(ncentrality2,"NCentrality_p2p_Gnutella04.csv")

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

fviz_eig(res.pca) #Scree plot
fviz_eig(res.pca2)#Scree plot 
fviz_eig(res.pca3)#Scree plot
fviz_eig(res.pca4)#Scree plot

biplot(res.pca,scale=0, cex=.7)#biplot
biplot(res.pca2,scale=0, cex=.7)#biplot
biplot(res.pca3,scale=0, cex=.7)#biplot
biplot(res.pca4,scale=0, cex=.7)#biplot

#hist(dg, breaks=5, main="With breaks=30") #histogram did not work

# Color by the quality of representation
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

ncen_tr<-transpose(ncentrality2) #transpose ncentrality for lda
ncen_tr<-data.frame(names = c('degree','eigenvector','closeness','pagerank','crossclique','betweeness','hubscore','authorities'),ncen_tr) #y label

#tsne model 4
set.seed(32)  
tsne_model_4 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.20, dims=2, max_iter = 2000,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
plot(tsne_model_4$Y,col=ncen_tr$names, asp=1)

#tsne model 5
set.seed(323)  
tsne_model_5 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
plot(tsne_model_5$Y,col=ncen_tr$names, asp=1)

set.seed(333)  
tsne_model_6 = Rtsne(pcacen2, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error

plot(tsne_model_6$Y,col=ncen_tr$names, asp=1)

d_tsne_1 = as.data.frame(tsne_model_5$Y) #list to dataframe
d_tsne_1_original=d_tsne_1 #keeping the original

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

dbscan::kNNdistplot(pcacen2, k =  2)
abline(h = 0.04, lty = 2)

set.seed(123)
res.db <- dbscan::dbscan(pcacen2, 0.04, 2)
fviz_cluster(res.db, pcacen2, geom = "point")

