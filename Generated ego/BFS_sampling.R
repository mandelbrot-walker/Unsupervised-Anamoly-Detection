library(igraph)
library(centiserve)
library(tidyverse)
library(data.table)
library(tibble)
library(gridExtra) 


edges<-read.delim("Email-EuAll.txt",header = TRUE, sep = "\t")

g<-graph.data.frame(edges) #  graph data frame for igraph

ncentrality2<-read.csv("NCentrality_Email_EuAll.csv") # read dataset of centralities

#---for rowsum highest----#

ncentrality2_rowsum_withnode<-cbind(ncentrality2, row_sum = rowSums(ncentrality2[-1]))
ncentrality2_rowsum_withnode<-ncentrality2_rowsum_withnode[order(-ncentrality2_rowsum_withnode$row_sum),]
source_node<- first(ncentrality2_rowsum_withnode$X)

x=30000 #number of nodes
f1 <- function(graph, data, extra) {
  data['rank'] == x
}

bfs_highest_rowsum_source  <- bfs(g, root=source_node, "all", order=TRUE, callback = f1) #run bfs

sub_highest_rowsum_source <-bfs_highest_rowsum_sourcek$order[1:x] 

newgraph_highest_rowsum_source = induced.subgraph(g, sub_highest_rowsum_source) #generate igraph of subgraph

bmp("ego_Plot_highest_rowsum-source_.bmp", width = 1920, height = 1080) #insert x value
plot(newgraph_highest_rowsum_source, main=" ego network of highest row sum source node nearby ... nodes") #insert x value
dev.off()


#---for highest contributing variable----#

ncentrality2_rowsum_withnode<-ncentrality2_rowsum_withnode[order(-ncentrality2_rowsum_withnode$betweenness),] #insert which varaible/centrality needed
source_node<- first(ncentrality2_rowsum_withnode$X)


bfs_rowsum_highestcontributing_variable_source  <- bfs(g, root=source_node, "all", order=TRUE, callback = f1) #run bfs

sub_rowsum_highestcontributing_variable_source <-bfs_rowsum_highestcontributing_variable_source$order[1:x]

newgraph_rowsum_highestcontributing_variable_source = induced.subgraph(g, sub_rowsum_highestcontributing_variable_source) # generate igraph of subgraph

bmp("ego_Plot_highest_betweenness_source_.bmp", width = 1920, height = 1080) #insert x value
plot(newgraph_rowsum_highestcontributing_variable_source_30k, main=" ego network of highest betweenness source nearby .... nodes") #insert x value
dev.off()


#------ for multiple graphs together---------#

bfs1  <- bfs(g, root=1, "all", order=TRUE, callback = f)
bfs2  <- bfs(g, root=1000, "all", order=TRUE, callback = f)
bfs3  <- bfs(g, root=2000, "all", order=TRUE, callback = f)
bfs4  <- bfs(g, root=3000, "all", order=TRUE, callback = f)
bfs5  <- bfs(g, root=4000, "all", order=TRUE, callback = f)
bfs6  <- bfs(g, root=5000, "all", order=TRUE, callback = f)


sub1 <- bfs1$order[1:x]
sub2 <- bfs2$order[1:x]
sub3 <- bfs3$order[1:x]
sub4 <- bfs4$order[1:x]
sub5 <- bfs5$order[1:x]
sub6 <- bfs6$order[1:x]


newgraph1 = induced.subgraph(g, sub1)
newgraph2 = induced.subgraph(g, sub2)
newgraph3 = induced.subgraph(g, sub3)
newgraph4 = induced.subgraph(g, sub4)
newgraph5 = induced.subgraph(g, sub5)
newgraph6 = induced.subgraph(g, sub6)


bmp("ego_plots_8k_vertices.bmp", width = 1920, height = 1080)
par(mfrow=c(3,2), mar=c(1,1,1,1))
plot(newgraph1, main="Source vertex=1")
plot(newgraph2, main="Source vertex=1000")
plot(newgraph3, main="Source vertex=2000")
plot(newgraph4, main="Source vertex=3000")
plot(newgraph5, main="Source vertex=4000")
plot(newgraph6, main="Source vertex=5000")
dev.off()