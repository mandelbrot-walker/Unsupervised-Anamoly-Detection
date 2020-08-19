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


edges<-read_csv("musae_PTBR_edges.csv") #Reading in musae_PTBR data

g<-graph.data.frame(edges) #graph data frame for igraph

transitivity(g)

ivs(g)

btn<-betweenness(g) # Calculation of Betweeness centrality
write.csv(btn, "btn_musae_PTBR.csv")

dg<-degree(g) # Calculation of Degre centrality
write.csv(dg, "dg_musae_PTBR.csv")

eig<-evcent(g)$vector # Calculation of Eigenvector centrality
write.csv(eig, "eig_musae_PTBR.csv")

clsn<-closeness(g) # Calculation of Closeness centrality
write.csv(clsn, "clsn_musae_PTBR.csv")

pgr<-page_rank(g)$vector #Calculation of Page Rank centrality
write.csv(pgr, "pgr_musae_PTBR.csv")

katz<-katzcent(g) # Could not calculate this due to insufficient RAM

crsc<-crossclique(g) # Calculation of Cross-Clique centrality
write.csv(crsc, "crsc_musae_PTBR.csv")

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
write.csv(centrality,"Centrality_musae_PTBR.csv")

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

write.csv(ncentrality,"NCentrality_musae_PTBR.csv")
write.csv(ncentrality2,"NCentrality_musae_PTBR.csv")

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