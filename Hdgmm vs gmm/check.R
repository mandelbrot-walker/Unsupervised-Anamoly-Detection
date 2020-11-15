library(Rtsne)
library(kernlab)
library(rgl)

df<-data.frame(new1 = rep(sample.int(10, 10, replace = T), each=2, length.out=500))

for (i in 2:500) {
  new<- rep(sample.int(10, 10, replace = T), each=2, length.out=500)
  df[,ncol(df)+1] <-new
  colnames(df)[ncol(df)] <- paste0("new", i)
}

tsne_m1 = Rtsne(df, check_duplicates=FALSE, pca=T, perplexity=5, theta=0.20, dims=2, max_iter = 1000,
                        verbose = TRUE, is_distance = FALSE, pca_center = T, pca_scale = T, num_threads = 5)

plot(tsne_m1$Y,col=factor(1:10))


tsne_m1_3d = Rtsne(df, check_duplicates=FALSE, pca=T, perplexity=5, theta=0.20, dims=3, max_iter = 1000,
                verbose = TRUE, is_distance = FALSE, pca_center = T, pca_scale = T, num_threads = 5)

plot3d(x=tsne_m1_3d$Y[,1],y=tsne_m1_3d$Y[,2],z=tsne_m1_3d$Y[,3],
              col=factor(1:10),
              type="s",radius=1)

library(HDclassif)

hdgmm<-hddc(df, K = 1:10, model = c("ALL"), threshold = 0.45,
            criterion = "bic", com_dim = 2, itermax = 50, eps = 0.001,
            algo = "EM", d_select = "Cattell", init = "param", show = getHDclassif.show(),scaling = TRUE,
            min.individuals = 5, noise.ctrl = 1e-08, mc.cores = 5,
            nb.rep = 2, keepAllRes = FALSE, d_max = 50, subset = Inf)

plot(tsne_m1$Y, col = hdgmm$class)

plot3d(x=tsne_m1_3d$Y[,1],y=tsne_m1_3d$Y[,2],z=tsne_m1_3d$Y[,3],
       col=hdgmm$class,
       type="s",radius=1)

table(hdgmm$class)

library(mclust)

gmm <- Mclust(as.matrix(df), prior = priorControl(), 
                    control = emControl(), 
                    warn = mclust.options("warn"),
                    verbose = TRUE)

summary(gmm, parameters = F)

plot(tsne_m1$Y, col = gmm$classification)

plot3d(x=tsne_m1_3d$Y[,1],y=tsne_m1_3d$Y[,2],z=tsne_m1_3d$Y[,3],
       col=gmm$classification,
       type="s",radius=1)
