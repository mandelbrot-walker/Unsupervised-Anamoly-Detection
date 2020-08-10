library(igraph)
library(centiserve)
library(tidyverse)
library(factoextra)
library(fastnet)
library(Rtsne)
library(Plasmidprofiler) #normalize function 
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


edges<-read_csv("Amazon0302.csv") #Reading in amazon0302 data

#might need to sample main dataset to address np hard issues

g<-graph.data.frame(edges) #graph data frame for igraph

cent<-proper_centralities(g) #returns 50 

#centrality feature decesion pending

transitivity(g) #for cross-clique centrality 

btn<-betweenness(g) # Calculation of Betweeness centrality
write.csv(btn, "btn_amazon0302.csv")

dg<-degree(g) # Calculation of Degre centrality
write.csv(dg, "dg_amazon0302.csv")

eig<-evcent(g)$vector # Calculation of Eigenvector centrality
write.csv(eig, "eig_amazon0302.csv")

clsn<-closeness(g) # Calculation of Closeness centrality
write.csv(clsn, "clsn_amazon0302.csv")

pgr<-page_rank(g)$vector #Calculation of Page Rank centrality
write.csv(pgr, "pgr_amazon0302.csv")

#katz<-katzcent(g) # Could not calculate this due to insufficient RAM

crsc<-crossclique(g) # Calculation of Cross-Clique centrality
write.csv(crsc, "crsc_amazon0302.csv")

frm<-closeness.freeman(g) # Not calculatable as graphis not strongly connected

edge_connectivity(g) # Outputs 0

cent<-proper_centralities(edges)

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
write.csv(centrality,"Centrality_amazon0302.csv")

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

write.csv(ncentrality,"NCentrality_amazon0302.csv")
write.csv(ncentrality2,"NCentrality_amazon0302.csv")

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
        notch = TRUE
) #multiple boxplot

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

#attempt at sampling

#tmp1<-resample(pcacen, samp.bootstrap(262111,5000), R = 5000, seed = 321,
#               block.size = 1000,trace = TRUE, observedIndices = 1:n)
#tmp<-bootstrap(pcacen, mean(pcacen), R = 50000,
#               seed = 321, sampler = samp.bootstrap,
#               block.size = 1000,
#               trace = TRUE)
#do.call(rbind.data.frame, tmp)
#coud not resample due to memory 

#attempt at lda

ncen_tr<-transpose(ncentrality) #transpose ncentrality for lda
ncen_tr<-data.frame(names = c('degree','eigenvector','closeness','pagerank','crossclique','betweeness'),ncen_tr) #y label

tmp<-ddply(ncen_tr, .(names), function(d,parallel=TRUE) { d[sample(nrow(d), pmin(nrow(d), 10000)),] }) #again tried resampling did not work

trranpcacen2<-transpose(ranpcacen2)
trranpcacen2<-data.frame(name = c('eigenvector','closeness','pagerank','crossclique','betweeness','hubscore','authorities'),trranpcacen2)

trsampcen<-transpose(sampcen)
trsampcen<-data.frame(name = c('degree','eigenvector','closeness','pagerank','crossclique','betweeness','hubscore','authorities'),trsampcen)

fit <- lda(trsampcen[,-1],grouping=trsampcen[,1], CV=TRUE) #10k data Fails cause: https://stat.ethz.ch/pipermail/r-help/2003-December/043870.html
fit # show results


## First attempt at tsne
set.seed(9)  
tsne_model_1 = Rtsne(as.matrix(ncentrality), check_duplicates=FALSE, pca=TRUE, perplexity=0.6, theta=0.5, dims=2)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)

## keeping original data
d_tsne_1_original=d_tsne_1

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1), 6)  
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1))) #doesn't work due to lack of memory

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


plot_k=plot_cluster(d_tsne_1_original, "cl_kmeans", "Accent") #plot of tsne using kmeans 
#plot_h=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")

#grid.arrange(plot_k, plot_k,  ncol=2)  
grid.arrange(plot_k,ncol=1) #plotting tnse model 1 with kmeans

#Custom tsne
#tsne model 1
set.seed(42)  
tsne_model_2 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=63, theta=0.5, dims=2,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 3)
plot(tsne_model_2$Y,col=ncen_tr$names, asp=1) #did not save the final error need to run again

plot(tsne_model_1$Y,col=ncen_tr$names, asp=1) #plot of tsne

d_tsne_2 = as.data.frame(tsne_model_2$Y) #list to dataframe
d_tsne_2_original=d_tsne_2 #keeping the original

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans_1=kmeans(scale(d_tsne_2), 6)  
d_tsne_2_original$cl_kmeans = factor(fit_cluster_kmeans_1$cluster)

plot_k1=plot_cluster(d_tsne_2_original, "cl_kmeans", "Accent")  

grid.arrange(plot_k1,ncol=1)  

#tsne model 3
set.seed(321)  
tsne_model_3 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=35, theta=0.5, dims=2,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 4)
plot(tsne_model_3$Y,col=ncen_tr$names, asp=1) #did not save the final error need to run again

#tsne model 4
set.seed(32)  
tsne_model_4 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.20, dims=2, max_iter = 1000,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
plot(tsne_model_4$Y,col=ncen_tr$names, asp=1)

#tsne model 5
set.seed(323)  
tsne_model_5 = Rtsne(ncentrality, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
plot(tsne_model_5$Y,col=ncen_tr$names, asp=1)

#tsne model 6 without degree
set.seed(333)  
tsne_model_6 = Rtsne(pcacen, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error

plot(tsne_model_6$Y,col=fornames$name, asp=1)

fornames<-transpose(pcacen)
fornames<-data.frame(name = c('eigenvector','closeness','pagerank','crossclique','betweeness'),fornames) #y label

set.seed(133)  
tsne_model_7 = Rtsne(pcacen, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 5) #lowest error

plot(tsne_model_7$Y,col=fornames$name, asp=1)

set.seed(13)  
tsne_model_8 = Rtsne(pcacen, check_duplicates=FALSE, pca=TRUE, perplexity=37, theta=0.10, dims=2, max_iter = 1500,
                     verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 5) #lowest error

plot(tsne_model_8$Y,col=fornames$name, asp=1)

d_tsne_3 = as.data.frame(tsne_model_5$Y) #list to dataframe
d_tsne_3_original=d_tsne_3 #keeping the original

#kmeans starts from here

fviz_nbclust(d_tsne_2,FUNcluster = kmeans,method = c("silhouette", "wss", "gap_stat"),diss = NULL,k.max = 10, nboot = 100,
             verbose = TRUE,barfill = "steelblue",barcolor = "steelblue",linecolor = 000000,
             print.summary = TRUE) #tried to find kneeplot for k

set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(ranpcacen2, k, nstart = 10 )$tot.withinss
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
c1<-kmeans(pcacen2, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(pcacen2, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(pcacen2, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(pcacen2, 6, iter.max = 20, nstart = 25,
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
c5<-kmeans(d_tsne_3, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c6<-kmeans(d_tsne_3, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)

p5 <- fviz_cluster(c5, geom = "point",  data = d_tsne_3) + ggtitle("k = 4")
p6 <- fviz_cluster(c6, geom = "point",  data = d_tsne_3) + ggtitle("k = 5")

grid.arrange(p5, p6, nrow = 1)

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans_1=kmeans(scale(d_tsne_3), 5)  
d_tsne_3_original$cl_kmeans = factor(fit_cluster_kmeans_1$cluster)

plot_t3=plot_cluster(d_tsne_3_original, "cl_kmeans", "Accent")  

grid.arrange(plot_t3,ncol=1)  

#tried dbscan here

library(fpc)
library(dbscan)

#kneeplot calculation 
dbscan::kNNdistplot(d_tsne_1, k =  7)
abline(h = 0.18, lty = 2)

dbscan::kNNdistplot(d_tsne_1, k =  7)
abline(h = 0.18, lty = 2)

#dbscan
set.seed(123)
# fpc package
res.fpc <- fpc::dbscan(d_tsne_1, eps = 0.18, MinPts = 2)
# dbscan package
res.db <- dbscan::dbscan(d_tsne_1, 0.18, 2)

#all(res.fpc$cluster == res.db)


fviz_cluster(res.fpc, d_tsne_1, geom = "point") #does not show expected result
#fviz_cluster(res.db, d_tsne_1, geom = "point")

#Random sampled dbscan
ranpcacen<-sample_n(pcacen, 50000, replace = TRUE)

dbscan::kNNdistplot(ranpcacen, k =  5)
abline(h = 0.017, lty = 2)

set.seed(123)
res.db <- dbscan::dbscan(ranpcacen, 0.017, 5)
fviz_cluster(res.db, ranpcacen, geom = "point")

res.db1 <- dbscan::dbscan(ranpcacen, 0.017, 5)
fviz_cluster(res.db1, ranpcacen, geom = "point")

ranpcacen2<-sample_n(pcacen2, 10000, replace = FALSE)#used 150k for dbscan
sampcen<-sample_n(centrality, 10000, replace = FALSE)

dbscan::kNNdistplot(ranpcacen2, k =  3)
abline(h = 0.016, lty = 2)
set.seed(123)
res.db2 <- dbscan::dbscan(ranpcacen2, 0.016, 3)
fviz_cluster(res.db2, ranpcacen2, geom = "point")

pairs(ranpcacen, col = res.db1$cluster + 1L)

#Mars model

# fit model
fit <- earth(trsampcen$name~., trsampcen,Use.beta.cache = TRUE)
# summarize the fit
summary(fit)
# summarize the importance of input variables
evimp(fit)
# make predictions
predictions <- predict(fit, trranpcacen2)
# summarize accuracy
mse <- mean((trranpcacen2$names - predictions)^2)
print(mse)

