library("igraph")
library("centiserve")
library(tidyverse)
library(factoextra)
library(fastnet)
library(Rtsne)



edges<-read_csv("musae_PTBR_edges.csv")
g<-graph.data.frame(edges)

pr_cent<-proper_centralities(g)

plot(g)

dg<-degree(g)

lay <- layout.fruchterman.reingold(g) # For plotting
lay2 <- layout_with_fr(g) # For plotting

plot.igraph(g, layout=lay, 
            vertex.size=dg*0.25,    # Rescaled by multiplying by 25
            main="Degree")
plot(g, layout = lay, vertex.label = NA)
plot(g, layout = lay2, vertex.label = NA)

sample1=sample_n(g,100,replace = FALSE)