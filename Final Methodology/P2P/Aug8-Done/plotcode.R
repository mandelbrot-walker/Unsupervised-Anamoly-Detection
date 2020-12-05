library(svglite)
library(ggplot2)
library(reshape2)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
melted_cormat <- melt(c)
head(melted_cormat)
ggplot(melted_cormat,aes(x=Var1, y=Var2, fill=value))+ 
  geom_tile()


# Get lower triangle of the correlation matrix
get_lower_tri<-function(c){
  c[upper.tri(c)] <- NA
  return(c)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(c){
  c[lower.tri(c)]<- NA
  return(c)
}

upper_tri <- get_upper_tri(c)
upper_tri

melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#999999", high = "#CC79A7", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggsave(file="test.svg", dpi= "print")
