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
edges<-read.delim("CA-GrQc.txt",header = TRUE, sep = "\t")

g<-graph.data.frame(edges) #graph data frame for igraph

transitivity(g) #  Check for cross clique calculation
transitivity(g, type = "localaverage")

# cent<-proper_centralities(g)
# 
# calculate_centralities(g, include = cent[1:50])%>%
#   pca_centralities(scale.unit = TRUE, ncp = 50) # takes indefinite time

dg<-degree(g) # Calculation of Degree centrality
write.csv(dg, "dg_gqrc.csv")

btn<-betweenness(g) # Calculation of Betweenness centrality
write.csv(btn, "btn_gqrc.csv")

eig<-evcent(g)$vector # Calculation of Eigenvector centrality
write.csv(eig, "eig_gqrc.csv")

clsn<-closeness(g) # Calculation of Closeness centrality
write.csv(clsn, "clsn_gqrc.csv")

pgr<-page_rank(g)$vector #Calculation of Page Rank centrality
write.csv(pgr, "pgr_gqrc.csv")

#katz<-katzcent(g) # Error in alpha >= maxEigenvalue : invalid comparison with complex values

crsc<-crossclique(g) # Calculation of Cross-Clique centrality
write.csv(crsc, "crsc_gqrc.csv")

frm<-closeness.freeman(g) # Not calculatable as graphis not strongly connected

edge_connectivity(g) # Outputs 0

gr<-g # temporary variable gr

V(gr)$degree <- dg                               #  Degree centrality
V(gr)$eig <- eig                                 #  Eigenvector centrality
V(gr)$closeness <- clsn                          #  Closeness centrality
V(gr)$pagerank <- pgr                            #  Pagerank centrality
V(gr)$betweenness <- btn                         #  Vertex betweenness centrality
V(gr)$crossclique <- crsc                        #  Crossclique centrality
V(gr)$hubs <- hub.score(g)$vector                #  Hub centrality
V(gr)$authorities <- authority.score(g)$vector   #  Authority centrality

centrality <- data.frame(row.names   = V(gr)$name,
                         degree      = V(gr)$degree,
                         eigenvector = V(gr)$eig,
                         closeness   = V(gr)$closeness,
                         pagerank    = V(gr)$pagerank,
                         crossclique = V(gr)$crossclique,
                         betweenness = V(gr)$betweenness,
                         hubscore    = V(gr)$hubs,
                         authorities = V(gr)$authorities
) #  Non-normalized centrality values

centrality <- centrality[order(row.names(centrality)),] #  These values are not normalized


head(centrality) #  check centrality variables

ndegree      = normalize(dg) 
neigenvector = normalize(eig) 
ncloseness   = normalize(clsn)
npagerank    = normalize(pgr)
ncrossclique = normalize(crsc)
nbetweenness = normalize(btn)
nhubscore    = normalize(V(gr)$hubs)
nauthorities = normalize(V(gr)$authorities)

ncentrality  <- data.frame(degree      = ndegree,
                           eigenvector = neigenvector,
                           closeness   = ncloseness,
                           pagerank    = npagerank,
                           crossclique = ncrossclique,
                           betweenness = nbetweenness,
                           hubscore    = nhubscore,
                           authorities = nauthorities 
) #  normalized values 8 variables


#------------------------------------------Data loader and centrality calculation End-------------------------------------#

var8<-ncentrality
var7<-within(ncentrality, rm(degree)) #  without degree total 7 variables
var6<-within(var7, rm(hubscore)) #  without hubscore total 6 variables
var5<-within(var6, rm(authorities)) #  without degree, authorities, hubscore total 5 variables
var6_degree<-within(ncentrality, rm(authorities,hubscore)) #  without authorities, hubscore total 6 variables

#------------------------------------------Data loader and centrality calculation End-------------------------------------#

#------------------------------------------Write data start-------------------------------------#
write.csv(centrality,"Centrality_gqrc.csv")
write.csv(ncentrality,"NCentrality_gqrc.csv")

#------------------------------------------Write data end-------------------------------------#

#------------------------------------------Boxplot and corelation matrix start-------------------------------------#
M <- cor(ncentrality)

#plot correlation matrix
bmp("corrplot.bmp", width = 1280, height = 720)
c<-corrplot(M, method = "circle") # correlation matrix 
dev.off()

#plot boxplot
bmp("boxplot.bmp", width = 1280, height = 720)
boxplot(ndegree, neigenvector, ncloseness, npagerank,ncrossclique,nbetweenness,nhubscore,nauthorities,
        main = "Multiple boxplots for comparision",
        at = c(1,2,3,4,5,6,7,8),
        names = c("degree", "eigenvector", "closeness", "pagerank","crossclique","betweenness","hubscore","authorities"),
        las = 2,
        col = c("orange","red"),
        border = "brown",
        horizontal = TRUE,
        notch = FALSE
) #multiple boxplot 
dev.off()
#------------------------------------------Boxplot and corelation matrix end-------------------------------------#

#------------------------------------------PCA start-------------------------------------#

#principal component analysis

res.pca6_degree<-prcomp(scale(var6_degree),center=TRUE) #  6 variables with degree 
res.pca5<-prcomp(scale(var5),center=TRUE) #  5 variables 
res.pca6<-prcomp(scale(var6),center=TRUE) #  6 variables
res.pca7<-prcomp(scale(var7),center=TRUE) #  7 variables 
res.pca8<-prcomp(scale(var8),center=TRUE) #  8 variables

#show pca values
print(res.pca6_degree)
print(res.pca5)
print(res.pca6)
print(res.pca7)
print(res.pca8)

p1<-fviz_eig(res.pca6_degree,addlabels = TRUE)#Scree plot 5 var
p2<-fviz_eig(res.pca5,addlabels = TRUE)#Scree plot 5 var
p3<-fviz_eig(res.pca6,addlabels = TRUE) #Scree plot 6 var
p4<-fviz_eig(res.pca7,addlabels = TRUE)#Scree plot 7 var
p5<-fviz_eig(res.pca8,addlabels = TRUE)#Scree plot 8 var

#plot screeplot
bmp("screeplot.bmp", width = 1920, height = 1080)
grid.arrange(p1, p2, p3, p4,p5, nrow = 2)
dev.off()

p1<-fviz_pca_biplot(res.pca6_degree, label ="var")#biplot 5 var
p2<-fviz_pca_biplot(res.pca5, label ="var")#biplot 5 var
p3<-fviz_pca_biplot(res.pca6, label ="var")#biplot 6 var
p4<-fviz_pca_biplot(res.pca7, label ="var")#biplot 7 var
p5<-fviz_pca_biplot(res.pca8, label ="var")#biplot 8 var

#plot biplot
bmp("biplot.bmp", width = 1920, height = 1080)
grid.arrange(p1, p2, p3, p4,p5, nrow = 2)
dev.off()


#plot contribution to dimensions
p1d1<-fviz_contrib(res.pca6_degree, choice = "var", axes = 1)
p1d2<-fviz_contrib(res.pca6_degree, choice = "var", axes = 2)
p2d1<-fviz_contrib(res.pca5, choice = "var", axes = 1)
p2d2<-fviz_contrib(res.pca5, choice = "var", axes = 2)
p3d1<-fviz_contrib(res.pca6, choice = "var", axes = 1)
p3d2<-fviz_contrib(res.pca6, choice = "var", axes = 2)
p4d1<-fviz_contrib(res.pca7, choice = "var", axes = 1)
p4d2<-fviz_contrib(res.pca7, choice = "var", axes = 2)
p5d1<-fviz_contrib(res.pca8, choice = "var", axes = 1)
p5d2<-fviz_contrib(res.pca8, choice = "var", axes = 2)

bmp("dimension contribution.bmp", width = 1920, height = 1080)
grid.arrange(p1d1,p1d2,p2d1,p2d2,p3d1,p3d2,p4d1,p4d2,p5d1,p5d2, nrow = 3, ncol = 4)
dev.off()

rm(p1d1,p1d2,p2d1,p2d2,p3d1,p3d2,p4d1,p4d2,p5d1,p5d2) #  remove contribution plot variables


# Color by contributions to the PC uses cos^2 value

p1<-fviz_pca_var(res.pca6_degree,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
) #  6 variables 


p2<-fviz_pca_var(res.pca5,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
) #  5 variables 

p3<-fviz_pca_var(res.pca6,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
) #  6 variables

p4<-fviz_pca_var(res.pca7,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
) #  7 variables 

p5<-fviz_pca_var(res.pca8,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
) #  8 variables

#plot contribution wise pca
#plot biplot
bmp("cos2 contribution.bmp", width = 1920, height = 1080)
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)
dev.off()

eig.val6_degree<-get_eigenvalue(res.pca6_degree)  #  gets eigenvalues 
eig.val5 <- get_eigenvalue(res.pca5) 
eig.val6 <- get_eigenvalue(res.pca6) 
eig.val7 <- get_eigenvalue(res.pca7) 
eig.val8 <- get_eigenvalue(res.pca8) 
eig.val6_degree #  shows eigenvalues 
eig.val5        
eig.val6
eig.val7
eig.val8


#  Results for Variables
res.var <- get_pca_var(res.pca8)#  8 variables 
res.var$coord          #  Coordinates
res.var$contrib        #  Contributions to the PCs
res.var$cos2           #  Quality of representation 

#------------------------------------------PCA end-------------------------------------#

save.image(".Rdata",safe = TRUE)

#------------------------------------------------Custom tsne---------------------------------------------------#

ncen_tr<-transpose(ncentrality) #  transpose ncentrality for tsne colors
ncen_tr<-data.frame(names = c('degree','eigenvector','closeness','pagerank','crossclique','betweeness','hubscore','authorities'),ncen_tr) #y label

#  color legend
colors_v8<-c("red","purple","blue","green","black","yellow","orange","magenta")
#  red = degree , purple = eigenvector, blue = closeness, green = pagerank, black = crossclique 
#  yellow = betweeness, orange = hubscore, magenta = authorities
colors_v7<-c("purple","blue","green","black","yellow","orange","magenta")
colors_v6<-c("purple","blue","green","black","yellow","magenta")
colors_v5<-c("purple","blue","green","black","yellow")
colors_v6_degree<-c("red","purple","blue","green","black","yellow","orange","magenta")
#  datasets: var6_degree var5 var6 var7 var8    

#-------------------------------------------------------tsne model 1 start--------------------------------------------#
set.seed(32)  
tsne_model_1_var6_degree = Rtsne(var6_degree, check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.20, dims=2, max_iter = 2000,
                                 verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_1_var6_degree.bmp", width = 1920, height = 1080)
plot(tsne_model_1_var6_degree$Y,col=colors_v6_degree, asp=1)
dev.off()

set.seed(32)  
tsne_model_1_var5 = Rtsne(var5, check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.20, dims=2, max_iter = 2000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_1_var5.bmp", width = 1920, height = 1080)
plot(tsne_model_1_var5$Y,col=colors_v5, asp=1)
dev.off()

set.seed(32)  
tsne_model_1_var6 = Rtsne(var6, check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.20, dims=2, max_iter = 2000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_1_var6.bmp", width = 1920, height = 1080)
plot(tsne_model_1_var6$Y,col=colors_v6, asp=1)
dev.off()

set.seed(32)  
tsne_model_1_var7 = Rtsne(var7, check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.20, dims=2, max_iter = 2000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_1_var7.bmp", width = 1920, height = 1080)
plot(tsne_model_1_var7$Y,col=colors_v7, asp=1)
dev.off()

set.seed(32)  
tsne_model_1_var8 = Rtsne(var8, check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.20, dims=2, max_iter = 2000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_1_var8.bmp", width = 1920, height = 1080)
plot(tsne_model_1_var8$Y,col=colors_v8, asp=1)
dev.off()

#-------------------------------------------------------tsne model 1 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 2 start--------------------------------------------#
set.seed(323)  
tsne_model_2_var6_degree = Rtsne(var6_degree, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                                 verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
bmp("tsne_model_2_var6_degree.bmp", width = 1920, height = 1080)
plot(tsne_model_2_var6_degree$Y,col=colors_v6_degree, asp=1)
dev.off()

set.seed(323)  
tsne_model_2_var5 = Rtsne(var5, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
bmp("tsne_model_2_var5.bmp", width = 1920, height = 1080)
plot(tsne_model_2_var5$Y,col=colors_v5, asp=1)
dev.off()

set.seed(323)  
tsne_model_2_var6 = Rtsne(var6, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
bmp("tsne_model_2_var6.bmp", width = 1920, height = 1080)
plot(tsne_model_2_var6$Y,col=colors_v6, asp=1)
dev.off()

set.seed(323)  
tsne_model_2_var7 = Rtsne(var7, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
bmp("tsne_model_2_var7.bmp", width = 1920, height = 1080)
plot(tsne_model_2_var7$Y,col=colors_v7, asp=1)
dev.off()

set.seed(323)  
tsne_model_2_var8 = Rtsne(var8, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 6) #lowest error
bmp("tsne_model_2_var8.bmp", width = 1920, height = 1080)
plot(tsne_model_2_var8$Y,col=colors_v8, asp=1)
dev.off()
#-------------------------------------------------------tsne model 2 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 3 start--------------------------------------------#
set.seed(333)  
tsne_model_3_var6_degree = Rtsne(var6_degree, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                                 verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7) #lowest error
bmp("tsne_model_3_var6_degree.bmp", width = 1920, height = 1080)
plot(tsne_model_3_var6_degree$Y,col=colors_v6_degree, asp=1)
dev.off()

set.seed(333)  
tsne_model_3_var5 = Rtsne(var5, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7) #lowest error
bmp("tsne_model_3_var5.bmp", width = 1920, height = 1080)
plot(tsne_model_3_var5$Y,col=colors_v5, asp=1)
dev.off()

set.seed(333)  
tsne_model_3_var6 = Rtsne(var6, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7) #lowest error
bmp("tsne_model_3_var6.bmp", width = 1920, height = 1080)
plot(tsne_model_3_var6$Y,col=colors_v6, asp=1)
dev.off()

set.seed(333)  
tsne_model_3_var7 = Rtsne(var7, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7) #lowest error
bmp("tsne_model_3_var7.bmp", width = 1920, height = 1080)
plot(tsne_model_3_var7$Y,col=colors_v7, asp=1)
dev.off()

set.seed(333)  
tsne_model_3_var8 = Rtsne(var8, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.50, dims=2, max_iter = 1000,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7) #lowest error
bmp("tsne_model_3_var8.bmp", width = 1920, height = 1080)
plot(tsne_model_3_var8$Y,col=colors_v8, asp=1)
dev.off()
#-------------------------------------------------------tsne model 3 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------------------tsne model 4 start--------------------------------------------#
set.seed(358)  
tsne_model_4_var6_degree = Rtsne(var6_degree, check_duplicates=FALSE, pca=TRUE, perplexity=43, theta=0.10, dims=2, max_iter = 1500,
                                 verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_4_var6_degree.bmp", width = 1920, height = 1080)
plot(tsne_model_4_var6_degree$Y,col=colors_v6_degree, asp=1)
dev.off()

set.seed(358)  
tsne_model_4_var5 = Rtsne(var5, check_duplicates=FALSE, pca=TRUE, perplexity=43, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_4_var5.bmp", width = 1920, height = 1080)
plot(tsne_model_4_var5$Y,col=colors_v5, asp=1)
dev.off()

set.seed(358)  
tsne_model_4_var6 = Rtsne(var6, check_duplicates=FALSE, pca=TRUE, perplexity=43, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_4_var6.bmp", width = 1920, height = 1080)
plot(tsne_model_4_var6$Y,col=colors_v6, asp=1)
dev.off()

set.seed(358)  
tsne_model_4_var7 = Rtsne(var7, check_duplicates=FALSE, pca=TRUE, perplexity=43, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_4_var7.bmp", width = 1920, height = 1080)
plot(tsne_model_4_var7$Y,col=colors_v7, asp=1)
dev.off()

set.seed(358)  
tsne_model_4_var8 = Rtsne(var5, check_duplicates=FALSE, pca=TRUE, perplexity=43, theta=0.10, dims=2, max_iter = 1500,
                          verbose = TRUE, is_distance = FALSE, pca_center = TRUE, pca_scale = TRUE, num_threads = 7)
bmp("tsne_model_4_var8.bmp", width = 1920, height = 1080)
plot(tsne_model_4_var8$Y,col=colors_v8, asp=1)
dev.off()

#-------------------------------------------------------tsne model 4 end--------------------------------------------#

save.image(".Rdata",safe = TRUE)
#d_tsne_1 = as.data.frame(tsne_model_2_8var$Y) #list to dataframe
#d_tsne_1_original=d_tsne_1 #keeping the original

#---------------------------------------------------------------------------------------------------------------------------#

save.image("backup.Rdata",safe = TRUE)

scatterplot3d(x=tsne_model_t$Y[,1],y=tsne_model_t$Y[,2],z=tsne_model_t$Y[,3],
              color = c("blue"), asp=1)

plot3d(x=tsne_model_t$Y[,1],y=tsne_model_t$Y[,2],z=tsne_model_t$Y[,3],
       col=c("red","purple","blue","green","black","yellow","orange","magenta"),
       type="s",radius=0.5)


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

wss <- function(k) {
  kmeans(var6, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

bmp("kmeans_kneeplot_var6.bmp", width = 1280, height = 720)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares") #kneeplot
dev.off()

wss <- function(k) {
  kmeans(var6_degree, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

bmp("kmeans_kneeplot_var6_degree.bmp", width = 1280, height = 720)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares") #kneeplot
dev.off()

#clusters using different k values var6_degree 
c1<-kmeans(var6_degree, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(var6_degree, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(var6_degree, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(var6_degree, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
#plot of clusters
p1 <- fviz_cluster(c1, geom = "point",  data = var6_degree) + ggtitle("k = 2")
p2 <- fviz_cluster(c2, geom = "point", data = var6_degree) + ggtitle("k = 3")
p3 <- fviz_cluster(c3, geom = "point",  data = var6_degree) + ggtitle("k = 4")
p4 <- fviz_cluster(c4, geom = "point",  data = var6_degree) + ggtitle("k = 5")

bmp("kmeans_pca_var6_degree.bmp", width = 1920, height = 1280)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

#clusters using different k values var5
c1<-kmeans(var5, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(var5, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(var5, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(var5, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
#plot of clusters
p1 <- fviz_cluster(c1, geom = "point",  data = var5) + ggtitle("k = 2")
p2 <- fviz_cluster(c2, geom = "point", data = var5) + ggtitle("k = 3")
p3 <- fviz_cluster(c3, geom = "point",  data = var5) + ggtitle("k = 4")
p4 <- fviz_cluster(c4, geom = "point",  data = var5) + ggtitle("k = 5")

bmp("kmeans_pca_var5.bmp", width = 1920, height = 1280)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

#clusters using different k values var6
c1<-kmeans(var6, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(var6, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(var6, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(var6, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
#plot of clusters
p1 <- fviz_cluster(c1, geom = "point",  data = var6) + ggtitle("k = 2")
p2 <- fviz_cluster(c2, geom = "point", data = var6) + ggtitle("k = 3")
p3 <- fviz_cluster(c3, geom = "point",  data = var6) + ggtitle("k = 4")
p4 <- fviz_cluster(c4, geom = "point",  data = var6) + ggtitle("k = 5")

bmp("kmeans_pca_var6.bmp", width = 1920, height = 1280)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

#clusters using different k values var7
c1<-kmeans(var7, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(var7, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(var7, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(var7, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
#plot of clusters
p1 <- fviz_cluster(c1, geom = "point",  data = var7) + ggtitle("k = 2")
p2 <- fviz_cluster(c2, geom = "point", data = var7) + ggtitle("k = 3")
p3 <- fviz_cluster(c3, geom = "point",  data = var7) + ggtitle("k = 4")
p4 <- fviz_cluster(c4, geom = "point",  data = var7) + ggtitle("k = 5")

bmp("kmeans_pca_var7.bmp", width = 1920, height = 1280)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

#clusters using different k values var8
c1<-kmeans(var8, 3, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(var8, 4, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(var8, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(var8, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
#plot of clusters
p1 <- fviz_cluster(c1, geom = "point",  data = var8) + ggtitle("k = 2")
p2 <- fviz_cluster(c2, geom = "point", data = var8) + ggtitle("k = 3")
p3 <- fviz_cluster(c3, geom = "point",  data = var8) + ggtitle("k = 4")
p4 <- fviz_cluster(c4, geom = "point",  data = var8) + ggtitle("k = 5")

bmp("kmeans_pca_var8.bmp", width = 1920, height = 1280)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

c_ncentrality<-kmeans(ncentrality, 5, iter.max = 20, nstart = 25,
                      algorithm = c("Hartigan-Wong"), trace=FALSE)

bmp("kmeans_ncentrality_k2.bmp", width = 1920, height = 1280)
plot(ncentrality, col = c_ncentrality$cluster)
points(c_ncentrality$centers, col = 1:8, pch = 8)
dev.off()

c_ncentrality<-kmeans(ncentrality, 6, iter.max = 20, nstart = 25,
                      algorithm = c("Hartigan-Wong"), trace=FALSE)

bmp("kmeans_ncentrality_k3.bmp", width = 1920, height = 1280)
plot(ncentrality, col = c_ncentrality$cluster)
points(c_ncentrality$centers, col = 1:8, pch = 8)
dev.off()

rm(p1,p2,p3,p4,p5)

round(calinhara(var8,c1$cluster),digits=3) #  Checking for correct no of clusters. Higher the index value better the cluster
round(calinhara(var8,c2$cluster),digits=4)
round(calinhara(var8,c3$cluster),digits=5)
round(calinhara(var8,c4$cluster),digits=6)

#-------------TSNE data frames for kmeans start--------------------------------------#
#~~~~~~~~~~~~~tsne model 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d_tsne_1_var6_degree = as.data.frame(tsne_model_1_var6_degree$Y) #list to dataframe var6 degree
d_tsne_1_var5 = as.data.frame(tsne_model_1_var5$Y) #list to dataframe var 5
d_tsne_1_var6 = as.data.frame(tsne_model_1_var6$Y) #list to dataframe var 6
d_tsne_1_var7 = as.data.frame(tsne_model_1_var7$Y) #list to dataframe var 7
d_tsne_1_var8 = as.data.frame(tsne_model_1_var8$Y) #list to dataframe var 8


#~~~~~~~~~~~~~tsne model 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d_tsne_2_var6_degree = as.data.frame(tsne_model_2_var6_degree$Y) #list to dataframe var6 degree
d_tsne_2_var5 = as.data.frame(tsne_model_2_var5$Y) #list to dataframe var 5
d_tsne_2_var6 = as.data.frame(tsne_model_2_var6$Y) #list to dataframe var 6
d_tsne_2_var7 = as.data.frame(tsne_model_2_var7$Y) #list to dataframe var 7
d_tsne_2_var8 = as.data.frame(tsne_model_2_var8$Y) #list to dataframe var 8

#~~~~~~~~~~~~~tsne model 3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d_tsne_3_var6_degree = as.data.frame(tsne_model_3_var6_degree$Y) #list to dataframe var6 degree
d_tsne_3_var5 = as.data.frame(tsne_model_3_var5$Y) #list to dataframe var 5
d_tsne_3_var6 = as.data.frame(tsne_model_3_var6$Y) #list to dataframe var 6
d_tsne_3_var7 = as.data.frame(tsne_model_3_var7$Y) #list to dataframe var 7
d_tsne_3_var8 = as.data.frame(tsne_model_3_var8$Y) #list to dataframe var 8

#~~~~~~~~~~~~~tsne model 4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

d_tsne_4_var6_degree = as.data.frame(tsne_model_4_var6_degree$Y) #list to dataframe var6 degree
d_tsne_4_var5 = as.data.frame(tsne_model_4_var5$Y) #list to dataframe var 5
d_tsne_4_var6 = as.data.frame(tsne_model_4_var6$Y) #list to dataframe var 6
d_tsne_4_var7 = as.data.frame(tsne_model_4_var7$Y) #list to dataframe var 7
d_tsne_4_var8 = as.data.frame(tsne_model_4_var8$Y) #list to dataframe var 8

#-------------TSNE data frames for kmeans end--------------------------------------#

#-------------kmeans on dataset and cluster onto TSNE start-------------------------------------------------#

#----------------------kmenas k=2------------------------------------#

c1<-kmeans(var6_degree, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c2<-kmeans(var5, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c3<-kmeans(var6, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c4<-kmeans(var7, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c5<-kmeans(var8, 5, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)

#----------------------kmenas k=3------------------------------------#

c6<-kmeans(var6_degree, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c7<-kmeans(var5, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c8<-kmeans(var6, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c9<-kmeans(var7, 6, iter.max = 20, nstart = 25,
           algorithm = c("Hartigan-Wong"), trace=FALSE)
c10<-kmeans(var8, 6, iter.max = 20, nstart = 25,
            algorithm = c("Hartigan-Wong"), trace=FALSE)

#Kmeans for tsne model 1

p1 <- fviz_cluster(c1, geom = "point",  data = d_tsne_1_var6_degree) + ggtitle("k = 2 var6 degree")
p2 <- fviz_cluster(c2, geom = "point",  data = d_tsne_1_var5) + ggtitle("k = 2 var5")
p3 <- fviz_cluster(c3, geom = "point",  data = d_tsne_1_var6) + ggtitle("k = 2 var6")
p4 <- fviz_cluster(c4, geom = "point",  data = d_tsne_1_var7) + ggtitle("k = 2 var7")
p5 <- fviz_cluster(c5, geom = "point",  data = d_tsne_1_var8) + ggtitle("k = 2 var8")


bmp("tsne_model1_kmeans_k2.bmp", width = 2560, height = 1280)
grid.arrange(p1, p2, p3,p4,p5, ncol = 5, nrow = 1)
dev.off()

p1 <- fviz_cluster(c6, geom = "point",  data = d_tsne_1_var6_degree) + ggtitle("k = 3 var6 degree")
p2 <- fviz_cluster(c7, geom = "point",  data = d_tsne_1_var5) + ggtitle("k = 3 var5")
p3 <- fviz_cluster(c8, geom = "point",  data = d_tsne_1_var6) + ggtitle("k = 3 var6")
p4 <- fviz_cluster(c9, geom = "point",  data = d_tsne_1_var7) + ggtitle("k = 3 var7")
p5 <- fviz_cluster(c10, geom = "point",  data = d_tsne_1_var8) + ggtitle("k = 3 var8")


bmp("tsne_model1_kmeans_k3.bmp", width = 2560, height = 1280)
grid.arrange(p1, p2, p3,p4,p5, ncol = 5, nrow = 1)
dev.off()

#Kmeans for tsne model 2

p1 <- fviz_cluster(c1, geom = "point",  data = d_tsne_2_var6_degree) + ggtitle("k = 2 var6 degree")
p2 <- fviz_cluster(c2, geom = "point",  data = d_tsne_2_var5) + ggtitle("k = 2 var5")
p3 <- fviz_cluster(c3, geom = "point",  data = d_tsne_2_var6) + ggtitle("k = 2 var6")
p4 <- fviz_cluster(c4, geom = "point",  data = d_tsne_2_var7) + ggtitle("k = 2 var7")
p5 <- fviz_cluster(c5, geom = "point",  data = d_tsne_2_var8) + ggtitle("k = 2 var8")


bmp("tsne_model2_kmeans_k2.bmp", width = 2560, height = 1280)
grid.arrange(p1, p2, p3,p4,p5, ncol = 5, nrow = 1)
dev.off()

p1 <- fviz_cluster(c6, geom = "point",  data = d_tsne_2_var6_degree) + ggtitle("k = 3 var6 degree")
p2 <- fviz_cluster(c7, geom = "point",  data = d_tsne_2_var5) + ggtitle("k = 3 var5")
p3 <- fviz_cluster(c8, geom = "point",  data = d_tsne_2_var6) + ggtitle("k = 3 var6")
p4 <- fviz_cluster(c9, geom = "point",  data = d_tsne_2_var7) + ggtitle("k = 3 var7")
p5 <- fviz_cluster(c10, geom = "point",  data = d_tsne_2_var8) + ggtitle("k = 3 var8")


bmp("tsne_model2_kmeans_k3.bmp", width = 2560, height = 1280)
grid.arrange(p1, p2, p3,p4,p5, ncol = 5, nrow = 1)
dev.off()

#Kmeans for tsne model 3

p1 <- fviz_cluster(c1, geom = "point",  data = d_tsne_3_var6_degree) + ggtitle("k = 2 var6 degree")
p2 <- fviz_cluster(c2, geom = "point",  data = d_tsne_3_var5) + ggtitle("k = 2 var5")
p3 <- fviz_cluster(c3, geom = "point",  data = d_tsne_3_var6) + ggtitle("k = 2 var6")
p4 <- fviz_cluster(c4, geom = "point",  data = d_tsne_3_var7) + ggtitle("k = 2 var7")
p5 <- fviz_cluster(c5, geom = "point",  data = d_tsne_3_var8) + ggtitle("k = 2 var8")


bmp("tsne_model3_kmeans_k2.bmp", width = 2560, height = 1280)
grid.arrange(p1, p2, p3,p4,p5, ncol = 5, nrow = 1)
dev.off()

p1 <- fviz_cluster(c6, geom = "point",  data = d_tsne_3_var6_degree) + ggtitle("k = 3 var6 degree")
p2 <- fviz_cluster(c7, geom = "point",  data = d_tsne_3_var5) + ggtitle("k = 3 var5")
p3 <- fviz_cluster(c8, geom = "point",  data = d_tsne_3_var6) + ggtitle("k = 3 var6")
p4 <- fviz_cluster(c9, geom = "point",  data = d_tsne_3_var7) + ggtitle("k = 3 var7")
p5 <- fviz_cluster(c10, geom = "point",  data = d_tsne_3_var8) + ggtitle("k = 3 var8")


bmp("tsne_model3_kmeans_k3.bmp", width = 2560, height = 1280)
grid.arrange(p1, p2, p3,p4,p5, ncol = 5, nrow = 1)
dev.off()

#Kmeans for tsne model 4

p1 <- fviz_cluster(c1, geom = "point",  data = d_tsne_4_var6_degree) + ggtitle("k = 2 var6 degree")
p2 <- fviz_cluster(c2, geom = "point",  data = d_tsne_4_var5) + ggtitle("k = 2 var5")
p3 <- fviz_cluster(c3, geom = "point",  data = d_tsne_4_var6) + ggtitle("k = 2 var6")
p4 <- fviz_cluster(c4, geom = "point",  data = d_tsne_4_var7) + ggtitle("k = 2 var7")
p5 <- fviz_cluster(c5, geom = "point",  data = d_tsne_4_var8) + ggtitle("k = 2 var8")


bmp("tsne_model4_kmeans_k2.bmp", width = 2560, height = 1280)
grid.arrange(p1, p2, p3,p4,p5, ncol = 5, nrow = 1)
dev.off()

p1 <- fviz_cluster(c6, geom = "point",  data = d_tsne_4_var6_degree) + ggtitle("k = 3 var6 degree")
p2 <- fviz_cluster(c7, geom = "point",  data = d_tsne_4_var5) + ggtitle("k = 3 var5")
p3 <- fviz_cluster(c8, geom = "point",  data = d_tsne_4_var6) + ggtitle("k = 3 var6")
p4 <- fviz_cluster(c9, geom = "point",  data = d_tsne_4_var7) + ggtitle("k = 3 var7")
p5 <- fviz_cluster(c10, geom = "point",  data = d_tsne_4_var8) + ggtitle("k = 3 var8")


bmp("tsne_model4_kmeans_k3.bmp", width = 2560, height = 1280)
grid.arrange(p1, p2, p3,p4,p5, ncol = 5, nrow = 1)
dev.off()

round(calinhara(tsne_model_1_var8$Y,c5$cluster),digits=5)
round(calinhara(tsne_model_1_var8$Y,c10$cluster),digits=6)

rm(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,p1,p2,p3,p4,p5)
rm(d_tsne_1_var6_degree,d_tsne_1_var5,d_tsne_1_var6, d_tsne_1_var7, d_tsne_1_var8)
rm(d_tsne_2_var6_degree,d_tsne_2_var5,d_tsne_2_var6, d_tsne_2_var7, d_tsne_2_var8)
rm(d_tsne_3_var6_degree,d_tsne_3_var5,d_tsne_3_var6, d_tsne_3_var7, d_tsne_3_var8)
rm(d_tsne_4_var6_degree,d_tsne_4_var5,d_tsne_4_var6, d_tsne_4_var7, d_tsne_4_var8)

#-----------------------kmeans on dataset and cluster onto TSNE end------------------------------------------------------#

save.image(".Rdata",safe = TRUE)

#-------------------------------------------dbscan start-----------------------------------------------------------------#

#---------------------------------------calculating h start------------------------------------------#
bmp("dbscan_kneeplot_var6_degree.bmp", width = 841, height = 477)
dbscan::kNNdistplot(var6_degree, k =  2)
abline(h = 0.001, lty = 2)
dev.off()

bmp("dbscan_kneeplot_var5.bmp", width = 841, height = 477)
dbscan::kNNdistplot(var5, k =  2)
abline(h = 0.001, lty = 2)
dev.off()

bmp("dbscan_kneeplot_var6.bmp", width = 841, height = 477)
dbscan::kNNdistplot(var6, k =  2)
abline(h = 0.001, lty = 2)
dev.off()

bmp("dbscan_kneeplot_var7.bmp", width = 841, height = 477)
dbscan::kNNdistplot(var7, k =  2)
abline(h = 0.001, lty = 2)
dev.off()

bmp("dbscan_kneeplot_var8.bmp", width = 841, height = 477)
dbscan::kNNdistplot(var8, k =  2)
abline(h = 0.001, lty = 2)
dev.off()
#---------------------------------------calculating h end------------------------------------------#

#---------------------------------------plotting dbscan start--------------------------------------------------#
#ranvar<-sample_n(var6_degree, 100000, replace = TRUE)
#set.seed(123)
res.db <- dbscan::dbscan(var6_degree, 0.067, 2)
gc()
bmp("dbscan_var6_degree.bmp", width = 1980, height = 1280)
fviz_cluster(res.db, var6_degree, geom = "point")
dev.off()

gc()
#ranvar<-sample_n(var5, 100000, replace = TRUE)
#set.seed(123)
res.db <- dbscan::dbscan(var5, 0.06, 2)
gc()
bmp("dbscan_var5.bmp", width = 1980, height = 1280)
fviz_cluster(res.db, var5, geom = "point")
dev.off()

gc()
#ranvar<-sample_n(var6, 100000, replace = TRUE)
#set.seed(123)
res.db <- dbscan::dbscan(var6, 0.06, 2)
gc()
bmp("dbscan_var6.bmp", width = 1980, height = 1280)
fviz_cluster(res.db, var6, geom = "point")
dev.off()

gc()
#ranvar<-sample_n(var7, 100000, replace = TRUE)
#set.seed(123)
res.db <- dbscan::dbscan(var7, 0.068, 2)
gc()
bmp("dbscan_var7.bmp", width = 1980, height = 1280)
fviz_cluster(res.db, var7, geom = "point")
dev.off()

gc()
#ranvar<-sample_n(var8, 100000, replace = TRUE)
#set.seed(123)
res.db <- dbscan::dbscan(var8, 0.068, 2)
gc()
bmp("dbscan_var8.bmp", width = 1980, height = 1280)
fviz_cluster(res.db, var8, geom = "point")
dev.off()

#---------------------------------------plotting dbscan end--------------------------------------------------#

#-----------------------------------------------dbscan end-----------------------------------------------------------------#

#-------------------TSNE_kmeans_without_convex_hull_start---------------------------#if necessary later

# plot_cluster=function(data, var_cluster, palette)
# {
#   ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
#     geom_point(size=0.25) +
#     guides(colour=guide_legend(override.aes=list(size=6))) +
#     xlab("") + ylab("") +
#     ggtitle("") +
#     theme_light(base_size=20) +
#     theme(axis.text.x=element_blank(),
#           axis.text.y=element_blank(),
#           legend.direction = "horizontal", 
#           legend.position = "bottom",
#           legend.box = "horizontal") + 
#     scale_colour_brewer(palette = palette) 
# }
# c5<-kmeans(var8, 2, iter.max = 20, nstart = 25,
#            algorithm = c("Hartigan-Wong"), trace=FALSE)
# 
# plot_k=plot_cluster(as.data.frame(tsne_model_1_var8$Y), factor(c5$cluster), "Accent")
# 
# grid.arrange(plot_k, nrow=1)

#--------------------TSNE_kmeans_without_convex_hull_end--------------------------#

#-------------------Dataset and TSNE dendogram---------------------------#
#d <- dist(var8, method = "euclidean")

# Hierarchical clustering using Complete Linkage
#hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
#plot(hc1, cex = 0.6, hang = -1)

#------------------Dataset and TSNE dendogram-----------------------------#

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

clustgmm <- Mclust(var8, G=3)

plot(clustgmm, what=c("classification"))
plot(clustgmm, "density")

#------------------GMM-----------------------------#

save.image(".Rdata",safe = TRUE)
save.image("backup.Rdata",safe = TRUE)

#----------------UMAP----------------------------#
library(uwot)
umap_calc <- umap(ranvar,
                  n_neighbors = 90,
                  min_dist = 10, spread = 50
)

umap_ranvar <- data.frame(
  UMAP1 = umap_calc[, 1],
  UMAP2 = umap_calc[, 2]
  #classification = c("eigenvector","closeness","pagerank","crossclique","betweenness")
)

ggplot(umap_ranvar, aes(
  x = UMAP1, y = UMAP2,
  col = classification
)) +
  geom_point()# +
#ggrepel::geom_text_repel(cex = 2.5)

plot(umap_ranvar$UMAP1,umap_ranvar$UMAP2, col=c22$cluster)
fviz_cluster(c22, geom = "point",  data = umap_ranvar) + ggtitle("k = 2 var5")
#----------------Umap---------------------------#

#--------------phate---------------------#
library(phateR)
library(ggplot2)
ranvar<-sample_n(var5, 80000, replace = FALSE)
x <- as.matrix(ranvar)
phate_calc <- phate(scale(x))

phate_ranvar5 <- data.frame(
  PHATE1 = phate_calc$embedding[, 1],
  PHATE2 = phate_calc$embedding[, 2],
  classification = c("eigenvector","closeness","pagerank","crossclique","betweenness")
)

ggplot(phate_ranvar5, aes(x = PHATE1, y = PHATE2, 
                          col = classification)) + 
  geom_point() 
plot(phate_ranvar5$PHATE1,phate_ranvar5$PHATE2, col=c2$cluster)

c22<-kmeans(ranvar, 2, iter.max = 20, nstart = 25,
            algorithm = c("Hartigan-Wong"), trace=FALSE)


fviz_cluster(c22, geom = "point",  data = phate_ranvar5) + ggtitle("k = 3 var6 degree")

#--------------phate---------------------#

#-------------autoencoders-------------=#
#-------------autoencoders-------------=#

#------------Sammon's map----------------#

x<-as.matrix(var6_degree) #var6_degree
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "pca")
bmp("sammon_var6_degree.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=colors_v6_degree, main="var6_degree")
plot(sammon$Y, pch=19, col=c6$cluster, main="var6_degree with k=2")
par(opar)
dev.off()

x<-as.matrix(var5) #var5
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "pca")
bmp("sammon_var5.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=colors_v5, main="var5")
plot(sammon$Y, pch=19, col=c7$cluster, main="var5 k=2")
par(opar)
dev.off()

x<-as.matrix(var6) #var6
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "pca")
bmp("sammon_var6.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=colors_v6, main="var6")
plot(sammon$Y, pch=19, col=c8$cluster, main="var6 k=2")
par(opar)
dev.off()

x<-as.matrix(var7) #var7
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "pca")
bmp("sammon_var7.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=colors_v7, main="var7")
plot(sammon$Y, pch=19, col=c9$cluster, main="var7 k=2")
par(opar)
dev.off()

x<-as.matrix(var8) #var8
sammon = do.sammon(x, ndim=2, preprocess = c("center"), initialize = "pca")
bmp("sammon_var8.bmp", width = 1280, height = 720)
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(sammon$Y, pch=19, col=colors_v8, main="var8")
plot(sammon$Y, pch=19, col=c10$cluster, main="var8 k=2")
par(opar)
dev.off()

rm(x,opar)

#------------Sammon's map----------------#

#------------NeRV----------------#
library(ProjectionBasedClustering)

#set.seed(100)
#ranvar<-sample_n(var5, 15000, replace = FALSE) #var5

projectionpoints=NeRV(as.matrix(var8), neighbors = 50, randominit = T, iterations = 5)
#Computation of Generalized Umatrix
library(GeneralizedUmatrix)
visualization=GeneralizedUmatrix(Data = as.matrix(var8),projectionpoints)
TopviewTopographicMap(visualization$Umatrix,visualization$Bestmatches)
#or in 3D if rgl package exists
library(rgl)
plotTopographicMap(visualization$Umatrix,visualization$Bestmatches)

##Interactive Island Generation 
## from a tiled Umatrix (toroidal assumption)
## Not run: 
Imx = ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(visualization$Umatrix,
                                                                     visualization$Bestmatches)
plotTopographicMap(visualization$Umatrix,visualization$Bestmatches, Imx = Imx)

## End(Not run)

# Automatic Clustering
LC=c(visualization$Lines,visualization$Columns)
# number of Cluster from dendrogram or visualization (PlotIt=TRUE)
Cls=ProjectionBasedClustering(k=2, as.matrix(var8),visualization$Bestmatches, LC,PlotIt=TRUE)
# Verification
TopviewTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls)
#or in 3D if rgl package exists
plotTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls)

## Sometimes you can improve a Clustering interactivly or mark additional Outliers manually
## Not run: 
Cls2 = interactiveClustering(visualization$Umatrix, visualization$Bestmatches, Cls)
#------------NeRV----------------#

#------------------SOM-----------------------------#

data_train <- var5
var5_n<-var5
data_train_matrix <- as.matrix(scale(var5))

coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}
pretty_palette <- c("#1f77b4","#ff7f0e","#2ca02c", "#d62728","#9467bd","#8c564b","#e377c2")
# Create the SOM Grid - you generally have to specify the size of the 
# training grid prior to training the SOM. Hexagonal and Circular 
# topologies are possible
som_grid <- somgrid(xdim = 50, ydim=50, topo="hexagonal")

# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(data_train_matrix, 
                 grid=som_grid, 
                 rlen=500, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE )

plot(som_model, type="changes")
plot(som_model, type="count", main="Node Counts")
plot(som_model, type="dist.neighbours", main = "SOM neighbour distances")
plot(som_model, type="codes")

plot(som_model, type = "property", property = getCodes(som_model)[,2], main=colnames(getCodes(som_model))[2], palette.name=coolBlueHotRed)

var <- 2 #define the variable to plot 
var_unscaled <- aggregate(as.numeric(data_train[,var]), by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,2] 
plot(som_model, type = "property", property=var_unscaled, main=names(data_train)[var], palette.name=coolBlueHotRed)

mydata <- as.data.frame(som_model$codes) 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
for (i in 2:15) {
  wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
}
plot(wss)

## use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(mydata)), 6)
# plot these results:
plot(som_model, type="mapping", bgcol = pretty_palette[som_cluster], main = "Clusters") 
add.cluster.boundaries(som_model, som_cluster)

# get vector with cluster value for each original data sample
cluster_assignment <- som_cluster[som_model$unit.classif]
# for each of analysis, add the assignment as a column in the original data:
var5_n$cluster <- cluster_assignment
#------------------SOM-----------------------------#


#  remove old variables
# rm(colors,xyMclust,tsne_model_4,tsne_model_5,tsne_model_6,tsne_model_7,
#   tsne_model_t,tsne_model_t5,tsne_model_t5_1,som_grid,res.db,res.pca,res.pca2,res.pca3,res.pca4,plot_h,plot_k,pcacen,pcacen2,
#   p1,p2,p3,p4,p5,p6,ncentrality,ncentrality2,ncen_tr,k14,hc1,gmm,fit_cluster_hierarchical,fit_cluster_kmeans,data_train_matrix,dat,
#   d_tsne_1_original,d_tsne_1,clustgmm,check,centrality,c1,c2,c3,c4,c5,c6)
