setwd("C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio")
getwd()
library("geomorph")

#select TPS file
TPS <- readland.tps("C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/Strictdata_Sept27.tps", specID ="ID")
TPS #contains all data with no removed landmarks

any(is.na(TPS))

#get dimensions
d<-dim(TPS)
p=d[1] #p number of landmarks
k=d[2] #k number of dimensions
n=d[3] #n number of specimens/individuals

#Perform Procrustes transformation
fishgpa <- gpagen (TPS, PrinAxes = TRUE) #GPA alignment
plot(TPS) #original coords
plot (fishgpa) #new coords with Procrustes

gdf_c<- fishgpa$coords
dim(gdf_c) #check this matches the original
dim(gdf_c)-d ##needs to be 0 0 0

#landmark linking
meanshape<- mshape(gdf_c)

landmarklinks<- define.links(meanshape)

#plots landmark coordinates for a set of specimens
plotAllSpecimens(gdf_c, links=landmarklinks)

#perform PCA
dev.off()
PCA_geomorph <- gm.prcomp(fishgpa$coords)
plot(PCA_geomorph , main = "PCA")
PCA_geomorph 

#Extract PCA values
shapecoords<-PCA_geomorph $shapes #min/max PCAs
PCAscores <-PCA$x #Component scores for each specimen

plot(PCA_geomorph $shapecoords$shapes.comp1, PCA_geomorph $shapecoords$shapes.comp2)
PCA_geomorph $rotation
PCA_geomorph $rotation[,1:2]
plot(PCA_geomorph $rotation[,1:2])

#PCA
plot(gm.prcomp(fishgpa$coords))

# Create data frame with trait start and end points
lmks <- data.frame(CP = c(6,8), BDP = c(5,9), BDA = c(3,10), PPFL = c(7,10), HD = c(11,21), HL=c(1,20),SL=c(1,16),EW=c(16,18),ML=c(14,15),
                   row.names = c("start", "end"))
#calculate linear distances
lineardists <- interlmkdist(TPS, lmks) #switched from gdf_c to TPS

#export linear distances to .csv
write.csv(lineardists, file = "C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/lineardistances_strictdata.csv") #should not have 2018

library(ggplot2)

#Import linear distances
dat<- read.csv("C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/lineardistances_strictdata_FL.csv")
dat
head(dat)

#remove rows with NA (FL was not measured due to bending in fish)
dat2 <- na.omit(dat) 
head(dat2)

#checking normality of raw data (size-corrected linear measurements) 
par(mfrow=c(3,3))
hist(dat2$CP)
hist(dat2$BDP)
hist(dat2$BDA)
hist(dat2$PPFL)
hist(dat2$HD)
hist(dat2$HL)
hist(dat2$SL)
hist(dat2$EW)
hist(dat2$ML)
dev.off() #turn off par function
#================================================
#ALLOMETRY GROWTH EQUATION
#================================================

#b for caudal peduncle (CP)
scatter.smooth(dat2$CP, dat2$ForkLength)
dat2$log.CP = log10(dat2$CP)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$log.CP)

#Fit the model summary (Mass.vs.CP)
Mass.vs.CP <- lm(log.FL~log.CP, data=dat2)
summary(Mass.vs.CP)

dat2$predicted.CP <- predict(Mass.vs.CP) # Save the predicted values
dat2$residuals.CP <- residuals(Mass.vs.CP)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for CP
beta_CP <-Mass.vs.CP[["coefficients"]][2]
dat2$beta_CP=beta_CP

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.CP = dat2$delta.length*dat2$beta_CP #b*(log10 Lm - log10 Li)

dat2$trans.CP = dat2$log.CP + dat2$delta.length.times.coefficient.CP  #log10 Mi (measured trait) + previous answer

dat2$trans.CP.FINAL=exp(dat2$trans.CP) #transform back from log to get final linear distance

#=================================
#b for body depth posterior (BDP)
scatter.smooth(dat2$BDP, dat2$ForkLength)
dat2$log.BDP = log10(dat2$BDP)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$log.BDP)

#Fit the model summary (Mass.vs.BDP)
Mass.vs.BDP <- lm(log.FL~log.BDP, data=dat2)
summary(Mass.vs.BDP)

dat2$predicted.BDP <- predict(Mass.vs.BDP) # Save the predicted values
dat2$residuals.BDP <- residuals(Mass.vs.BDP)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for BDP
beta_BDP <-Mass.vs.BDP[["coefficients"]][2]
dat2$beta_BDP=beta_BDP

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.BDP = dat2$delta.length*dat2$beta_BDP #b*(log10 Lm - log10 Li)

dat2$trans.BDP = dat2$log.BDP + dat2$delta.length.times.coefficient.BDP  #log10 Mi (measured trait) + previous answer

dat2$trans.BDP.FINAL=exp(dat2$trans.BDP) #transform back from log to get final linear distance
#=================================
#b for postpelvic fin length (BDA)
scatter.smooth(dat2$BDA, dat2$ForkLength)
dat2$log.BDA = log10(dat2$BDA)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$log.BDA)

#Fit the model summary (Mass.vs.BDA)
Mass.vs.BDA <- lm(log.FL~log.BDA, data=dat2)
summary(Mass.vs.BDA)

dat2$predicted.BDA <- predict(Mass.vs.BDA) # Save the predicted values
dat2$residuals.BDA <- residuals(Mass.vs.BDA)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for BDA
beta_BDA <-Mass.vs.BDA[["coefficients"]][2]
dat2$beta_BDA=beta_BDA

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.BDA = dat2$delta.length*dat2$beta_BDA #b*(log10 Lm - log10 Li)

dat2$trans.BDA = dat2$log.BDA + dat2$delta.length.times.coefficient.BDA #log10 Mi (measured trait) + previous answer

dat2$trans.BDA.FINAL=exp(dat2$trans.BDA) #transform back from log to get final linear distance

#=================================
#b for postpelvic fin length (PPFL)
scatter.smooth(dat2$PPFL, dat2$ForkLength)
dat2$log.PPFL = log10(dat2$PPFL)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$log.PPFL)

#Fit the model summary (Mass.vs.PPFL)
Mass.vs.PPFL <- lm(log.FL~log.PPFL, data=dat2)
summary(Mass.vs.PPFL)

dat2$predicted.PPFL <- predict(Mass.vs.PPFL) # Save the predicted values
dat2$residuals.PPFL <- residuals(Mass.vs.PPFL)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for PPFL
beta_PPFL <-Mass.vs.PPFL[["coefficients"]][2]
dat2$beta_PPFL=beta_PPFL

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.PPFL = dat2$delta.length*dat2$beta_PPFL #b*(log10 Lm - log10 Li)

dat2$trans.PPFL = dat2$log.PPFL + dat2$delta.length.times.coefficient.PPFL  #log10 Mi (measured trait) + previous answer

dat2$trans.PPFL.FINAL=exp(dat2$trans.PPFL) #transform back from log to get final linear distance


#=================================
#b for head depth (HD)
scatter.smooth(dat2$HD, dat2$ForkLength)
dat2$log.HD = log10(dat2$HD)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$log.HD)

#Fit the model summary (Mass.vs.HD)
Mass.vs.HD <- lm(log.FL~log.HD, data=dat2)
summary(Mass.vs.HD)

dat2$predicted.HD <- predict(Mass.vs.HD) # Save the predicted values
dat2$residuals.HD <- residuals(Mass.vs.HD)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for HD
beta_HD <-Mass.vs.HD[["coefficients"]][2]
dat2$beta_HD=beta_HD

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.HD = dat2$delta.length*dat2$beta_HD #b*(log10 Lm - log10 Li)

dat2$trans.HD = dat2$log.HD + dat2$delta.length.times.coefficient.HD  #log10 Mi (measured trait) + previous answer

dat2$trans.HD.FINAL=exp(dat2$trans.HD) #transform back from log to get final linear distance

#=================================
#b for head length (HL)
scatter.smooth(dat2$HD, dat2$ForkLength)
dat2$log.HL = log10(dat2$HL)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$log.HL)

#Fit the model summary (Mass.vs.HL)
Mass.vs.HL <- lm(log.FL~log.HL, data=dat2)
summary(Mass.vs.HL)

dat2$predicted.HL <- predict(Mass.vs.HL) # Save the predicted values
dat2$residuals.HL <- residuals(Mass.vs.HL)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for HL
beta_HL <-Mass.vs.HL[["coefficients"]][2]
dat2$beta_HL=beta_HL

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.HL = dat2$delta.length*dat2$beta_HL #b*(log10 Lm - log10 Li)

dat2$trans.HL = dat2$log.HL + dat2$delta.length.times.coefficient.HL  #log10 Mi (measured trait) + previous answer

dat2$trans.HL.FINAL=exp(dat2$trans.HL) #transform back from log to get final linear distance

#=================================
#b for snout length (SL)
scatter.smooth(dat2$SL, dat2$ForkLength)
dat2$log.SL = log10(dat2$SL)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$log.SL)

#Fit the model summary (Mass.vs.SL)
Mass.vs.SL <- lm(log.FL~log.SL, data=dat2)
summary(Mass.vs.SL)

dat2$predicted.SL <- predict(Mass.vs.SL) # Save the predicted values
dat2$residuals.SL <- residuals(Mass.vs.SL)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for SL
beta_SL <-Mass.vs.SL[["coefficients"]][2]
dat2$beta_SL=beta_SL

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.SL = dat2$delta.length*dat2$beta_SL #b*(log10 Lm - log10 Li)

dat2$trans.SL = dat2$log.SL + dat2$delta.length.times.coefficient.SL  #log10 Mi (measured trait) + previous answer

dat2$trans.SL.FINAL=exp(dat2$trans.SL) #transform back from log to get final linear distance

#=================================
#b for eye width (EW)
scatter.smooth(dat2$EW, dat2$ForkLength)
dat2$log.EW = log10(dat2$EW)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$log.EW)

#Fit the model summary (Mass.vs.EW)
Mass.vs.EW <- lm(log.FL~log.EW, data=dat2)
summary(Mass.vs.EW)

dat2$predicted.EW <- predict(Mass.vs.EW) # Save the predicted values
dat2$residuals.EW <- residuals(Mass.vs.EW)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for EW
beta_EW <-Mass.vs.EW[["coefficients"]][2]
dat2$beta_EW=beta_EW

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.EW = dat2$delta.length*dat2$beta_EW #b*(log10 Lm - log10 Li)

dat2$trans.EW = dat2$log.EW + dat2$delta.length.times.coefficient.EW  #log10 Mi (measured trait) + previous answer

dat2$trans.EW.FINAL=exp(dat2$trans.EW) #transform back from log to get final linear distance

#=================================
#b for maxilla length (ML)
scatter.smooth(dat2$ML, dat2$ForkLength)
dat2$log.ML = log10(dat2$ML)
dat2$log.FL = log10(dat2$ForkLength)
scatter.smooth(dat2$log.FL, dat2$ML)

#Fit the model summary (Mass.vs.ML)
Mass.vs.ML <- lm(log.FL~log.ML, data=dat2)
summary(Mass.vs.ML)

dat2$predicted.ML <- predict(Mass.vs.ML) # Save the predicted values
dat2$residuals.ML <- residuals(Mass.vs.ML)

#Mean ForkLength
dat2$mean.FL = mean(dat2$ForkLength)
dat2$log.mean.FL = log10(dat2$mean.FL)

#Calculate b for ML
beta_ML <-Mass.vs.ML[["coefficients"]][2]
dat2$beta_ML=beta_ML

#Create new column & put equation together
dat2$delta.length=dat2$log.mean.FL-dat2$log.FL  #log10 Lm (log Average FL) - log10 Li (Log FL)

dat2$delta.length.times.coefficient.ML = dat2$delta.length*dat2$beta_ML #b*(log10 Lm - log10 Li)

dat2$trans.ML = dat2$log.ML + dat2$delta.length.times.coefficient.ML  #log10 Mi (measured trait) + previous answer

dat2$trans.ML.FINAL=exp(dat2$trans.ML) #transform back from log to get final linear distance


#create dataframe with adjusted linear measurements
dat3 <-cbind.data.frame(dat2$X,dat2$trans.CP.FINAL, dat2$trans.BDP.FINAL, dat2$trans.BDA.FINAL, dat2$trans.PPFL.FINAL, dat2$trans.HD.FINAL, dat2$trans.HL.FINAL, dat2$trans.SL.FINAL, dat2$trans.EW.FINAL, dat2$trans.ML.FINAL) 
dat3

dat4 <- dat3[,-c(1,11)] #df with only the linear measurements and no ID (used for PCA)
nrow(dat4)
#================================================
#PCA on Linear Distances
#================================================
#PCA reference: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

library(FactoMineR)
library(factoextra)
library(cluster)

#Incorporate classifiers
#set factors
dat3$Year <- factor (dat2$Year, levels=c("2018", "2019"))
is.factor(dat3$Year)
head(dat3)

#PCA
# Interpret data https://www.youtube.com/watch?v=pks8m2ka7Pk&ab_channel=Fran%C3%A7oisHusson

#PCA 
dat4
res.pca <- PCA(dat4, graph = FALSE)
res.pca
summary(res.pca)

res.pca$ind$coord #individual coordinates

# Extract eigenvalues/variances
PCA_eigen_strictdata<- get_eig(res.pca)
#export cluster data to .csv
write.csv(PCA_eigen_strictdata, file = "C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/PCA_eigenvalues_strictdata.csv")

# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

# Extract the results for variables
var <- get_pca_var(res.pca)
var

# Coordinates of variables
var$coord

# Contribution of variables
var$contrib

# Graph of variables: default plot
fviz_pca_var(res.pca, col.var = "black")

                                #Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
#Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

#Contributions of variables to PC1 nice plot
Dim.1.plot <- data.frame(
  name=factor(c("Head Length", "Head Depth", "Snout Length", "Eye Width", "Maxilla Length", "Body Depth Anterior", "Body Depth Posterior", "Postpelvic Fin Length", "Caudal Peduncle Depth")),
  value=c(24.09, 19.01, 13.90, 13.54, 13.06, 8.48, 4.83, 2.35, 0.74)
)

ggplot(Dim.1.plot, aes(x=name,y=value))+
  theme_minimal()+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 25, hjust=1, size=8), axis.title.x=element_text(size=9, face="bold"), axis.title.y = element_text(size=9, face="bold"))+
  labs(x="Linear Measurement", y="Percent Contribution (%)")+
  scale_x_discrete(limits=c("Head Length", "Head Depth", "Snout Length", "Eye Width", "Maxilla Length", "Body Depth Anterior", "Body Depth Posterior", "Postpelvic Fin Length", "Caudal Peduncle Depth"))+
  geom_text(aes(label = value), vjust = -0.4, size=3.2)
  

#Extract the results for individuals
ind <- get_pca_ind(res.pca)
ind
ind$coord

#Create dataframe with fish ID and individual coords from PCA
indcoord <- cbind.data.frame(dat2$X, ind$coord)
head(indcoord)


min(indcoord$Dim.1, "dat2$X")
max(indcoord$Dim.1)
min(indcoord$Dim.2, "dat2$X")
max(indcoord$Dim.2)

#use silhouette method to determine no. of clusters
set.seed(10)
nrow(res.pca$ind$coord)
fviz_nbclust(
  res.pca$ind$coord,
  FUNcluster = kmeans,
  method = "silhouette",
  k.max = 8,
  verbose = interactive(),
  barfill = "steelblue",
  barcolor = "steelblue",
  linecolor = "steelblue",
  print.summary = TRUE,
)

#form clusters in PCA
set.seed(20)
k <- kmeans(res.pca$ind$coord, centers = 2, nstart = 25)
class(k)
k$cluster

#display cluster silhouette coefficients
#>0 means the observation is well grouped. The closer the coefficient is to 1, the better the observation is grouped.
#<0 means  the observation has been placed in the wrong cluster.
#=0 means the observation is between two clusters.
sil <- silhouette(k$cluster, dist(res.pca$ind$coord))
fviz_silhouette(sil)
dat3$Year


#Extract individuals to corresponding clusters
k$cluster<-as.factor(k$cluster)
is.factor(k$cluster)
ind_clusters <- cbind.data.frame(dat2$X, k$cluster, indcoord[,2:6])
ind_clusters

colnames(ind_clusters)[1] <- "ID"
colnames(ind_clusters)[2] <- "cluster"
as.factor(ind_clusters$cluster)
ind_clusters

#cluster plot
base <- ggplot (ind_clusters, aes(x=Dim.1,y=Dim.2, color=cluster)) +
geom_point(aes(x=Dim.1, y=Dim.2), size=2, pch=19)+
 labs(x="PC1 (38.0%)",y= "PC2 (26.8%)",col="")+  
  theme(panel.background = element_rect(fill = "transparent"))+
  theme(plot.background = element_rect(fill = "transparent"))+
  theme_classic()+
  scale_x_continuous(waiver())+
scale_y_continuous(waiver())+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=14, vjust=-2),
        axis.title.y = element_text(size=14, vjust=2),
        legend.text=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12))+
  stat_ellipse(aes(color=cluster), type="norm", alpha=1, lwd = 1.1)+
  scale_color_manual(labels=c("Cluster 1", "Cluster 2"), values = c("#276695", "#B47B08"))+ 
theme(plot.margin = margin(1,1,1,1,"cm"))

  

  
  
#export cluster data to .csv
write.csv(ind_clusters, file = "C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/clusters_strictdata.csv")


#=====================================
#Visualize PCA data with TPS grids
#resource: https://www.r-bloggers.com/2015/03/tips-tricks-7-plotting-pca-with-tps-grids/
library(ggrepel)

#procrustes coordinates only
#Extract PCA values
ind_clusters

#get ID names for graph label
ID_cleaned<- read.csv("C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/IDnames_cleaned.csv")
head(ID_cleaned)

ind_clusters <- cbind(ind_clusters, ID_cleaned[3])

#plot PCA data in ggplot
#par(fig=c(0.2, 0.7, 0.3, 0.7), new=TRUE)
dev.off()
plot3 <- ggplot(ind_clusters, aes(x=Dim.1, y=Dim.2, color=cluster))+
  geom_point(aes(x=Dim.1, y=Dim.2), size=2, pch=19) + 
  scale_color_manual(labels=c("Cluster 1", "Cluster 2"), values = c("#4682B4", "#B47846"))+
     labs(x="",y= "",col="Cluster")+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  theme_minimal()+
  scale_x_continuous(waiver())+
  scale_y_continuous(waiver())+
  theme(legend.position = c(1.1,0.75),
        legend.direction = "vertical",
        legend.title = element_blank())+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)+
  theme(axis.ticks = element_line(size = 1) , 
    axis.ticks.length = unit(.1, "cm"))+
  theme(plot.margin = margin(0, 4.5,0, 0,"cm"))+
  stat_ellipse(aes(group = cluster))

shift_axis_y <- function(plot3, y=0){
  g <- ggplotGrob(plot3)
  dummy <- data.frame(y=y)
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  plot3 + annotation_custom(grid::grobTree(ax, vp = grid::viewport(y=1, height=sum(ax$height))), 
                        ymax=y, ymin=y) +
    geom_hline(aes(yintercept=y), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank())
}

shift_axis_x <- function(plot3, x=0){
  g <- ggplotGrob(plot3)
  dummy <- data.frame(x=x)
  ax <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  plot3 + annotation_custom(grid::grobTree(ax, vp = grid::viewport(x=1, width = sum(ax$height))), 
                        xmax=x, xmin=x) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())
}

plot3<-shift_axis_y(plot3, y=0)
plot3<-shift_axis_x(plot3, x=0)
plot3

ggsave("cluster.png", plot3, height = 5, width = 7, dpi = 600)

ref <- mshape(fishgpa$coords)

GP1 <- gridPar(tar.pt.bg = "black", tar.link.lwd=1.5, grid.col = "grey")

#make TPS grids
plotRefToTarget(shapecoords$shapes.comp1$min, ref, gridPars=GP1, links=landmarklinks, mag=2) 
plotRefToTarget(shapecoords$shapes.comp1$max, ref, gridPars=GP1, links=landmarklinks, mag=2)
plotRefToTarget(shapecoords$shapes.comp2$min, ref, gridPars=GP1, links=landmarklinks, mag=2)
plotRefToTarget(shapecoords$shapes.comp2$max, ref, gridPars=GP1, links=landmarklinks, mag=2)

#label
#  geom_label_repel(aes(label = ID_cleaned),
box.padding   = 0.35, 
point.padding = 0.5,
segment.color = 'grey50')+

#==================
#MANOVA on clusters
#==================
library(mvnormtest)

#REMOVEUNTIL LINE 596
cluster_indpca2 <- cbind(cluster_indpca, dat3$Year)
cluster_indpca2

#check classes for each variable
lapply(cluster_indpca2, class)
cluster_indpca2$Year<-as.factor(cluster_indpca2[,10])
cluster_indpca2

cluster_indpca3 <- cluster_indpca2[ -c(7, 8, 10)] #remove duplicated columns
ls(cluster_indpca3)
cluster_indpca3

binddata <-cbind(cluster_indpca,cluster_indpca3[,2:6])

#check assumptions - normality, equal variances, independence, random sampling
#normality
shapiro.test(ind_clusters$Dim.1)
shapiro.test(ind_clusters$Dim.2)
shapiro.test(ind_clusters$Dim.3)
shapiro.test(ind_clusters$Dim.4)
shapiro.test(ind_clusters$Dim.5)
shapiro.test(ind_clusters$Dim.6)
shapiro.test(ind_clusters$Dim.7)
shapiro.test(ind_clusters$Dim.8)
shapiro.test(ind_clusters$Dim.9)

#equal variances
#Ho: variance-covariance matrix is same across groups
#Ha: variance-covariance matrix is not same across groups
boxM(ind_clusters[,3:11],ind_clusters$cluster)
boxM

manova <- manova(cbind(Dim.1, Dim.2, Dim.3, Dim.4, Dim.4, Dim.5) ~ k$cluster, data=bind_data)
manova
summary(manova)
summary.aov(manova)


#get dataframe with linear distances and cluster groups
ind_clusters.lineardist <- cbind.data.frame(dat3, k$cluster)
head(ind_clusters.lineardist)

#rename columns
colnames(ind_clusters.lineardist)[1] <- "ID"
colnames(ind_clusters.lineardist)[2] <- "trans.CP.FINAL"
colnames(ind_clusters.lineardist)[3] <- "trans.BDP.FINAL"
colnames(ind_clusters.lineardist)[4] <- "trans.BDA.FINAL"
colnames(ind_clusters.lineardist)[5] <- "trans.PPFL.FINAL"
colnames(ind_clusters.lineardist)[6] <- "trans.HD.FINAL"
colnames(ind_clusters.lineardist)[7] <- "trans.HL.FINAL"
colnames(ind_clusters.lineardist)[8] <- "trans.SL.FINAL"
colnames(ind_clusters.lineardist)[9] <- "trans.EW.FINAL"
colnames(ind_clusters.lineardist)[10] <- "trans.ML.FINAL"
colnames(ind_clusters.lineardist)[11] <- "cluster"

ind_clusters.lineardist

#ASSUMPTIONS for MANOVA
#Normality test
shapiro.test(ind_clusters.lineardist$trans.CP.FINAL) #p >0.05 to assume normality
shapiro.test(ind_clusters.lineardist$trans.BDP.FINAL)
shapiro.test(ind_clusters.lineardist$trans.BDA.FINAL)
shapiro.test(ind_clusters.lineardist$trans.PPFL.FINAL)
shapiro.test(ind_clusters.lineardist$trans.HD.FINAL)
shapiro.test(ind_clusters.lineardist$trans.HL.FINAL)
shapiro.test(ind_clusters.lineardist$trans.SL.FINAL)
shapiro.test(ind_clusters.lineardist$trans.EW.FINAL)
shapiro.test(ind_clusters.lineardist$trans.ML.FINAL) #not all are normal, must log data

hist(ind_clusters.lineardist$trans.CP.FINAL) #histograms
hist(ind_clusters.lineardist$trans.BDP.FINAL)
hist(ind_clusters.lineardist$trans.BDA.FINAL)
hist(ind_clusters.lineardist$trans.PPFL.FINAL)
hist(ind_clusters.lineardist$trans.HD.FINAL)
hist(ind_clusters.lineardist$trans.HL.FINAL)
hist(ind_clusters.lineardist$trans.SL.FINAL)
hist(ind_clusters.lineardist$trans.EW.FINAL)
hist(ind_clusters.lineardist$trans.ML.FINAL)

#equal variance test
leveneTest(trans.CP.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist) #we want p>0.05 so data is not significant (homogenous)
leveneTest(trans.BDP.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist)
leveneTest(trans.BDA.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist)
leveneTest(trans.PPFL.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist)
leveneTest(trans.HD.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist)
leveneTest(trans.HL.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist)
leveneTest(trans.SL.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist)
leveneTest(trans.EW.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist)
leveneTest(trans.ML.FINAL ~ as.factor(cluster), data = ind_clusters.lineardist)

manova2 <- manova(cbind(dat2$trans.CP.FINAL, dat2$trans.BDP.FINAL, dat2$trans.BDA.FINAL, dat2$trans.PPFL.FINAL, dat2$trans.HD.FINAL, dat2$trans.HL.FINAL, dat2$trans.SL.FINAL, dat2$trans.EW.FINAL, dat2$trans.ML.FINAL) ~ k$cluster, data=ind_clusters.lineardist)
manova2
summary(manova2) #statistically significant

#summary statistics for linear measurements
mean(ind_clusters.lineardist$trans.CP.FINAL)
mean(ind_clusters.lineardist$trans.BDP.FINAL)
mean(ind_clusters.lineardist$trans.BDA.FINAL)
mean(ind_clusters.lineardist$trans.PPFL.FINAL)
mean(ind_clusters.lineardist$trans.HD.FINAL)
mean(ind_clusters.lineardist$trans.HL.FINAL)
mean(ind_clusters.lineardist$trans.SL.FINAL)
mean(ind_clusters.lineardist$trans.EW.FINAL)
mean(ind_clusters.lineardist$trans.ML.FINAL)

std.error(ind_clusters.lineardist$trans.CP.FINAL)
std.error(ind_clusters.lineardist$trans.BDP.FINAL)
std.error(ind_clusters.lineardist$trans.BDA.FINAL)
std.error(ind_clusters.lineardist$trans.PPFL.FINAL)
std.error(ind_clusters.lineardist$trans.HD.FINAL)
std.error(ind_clusters.lineardist$trans.HL.FINAL)
std.error(ind_clusters.lineardist$trans.SL.FINAL)
std.error(ind_clusters.lineardist$trans.EW.FINAL)
std.error(ind_clusters.lineardist$trans.ML.FINAL)

min(ind_clusters.lineardist$trans.CP.FINAL)
min(ind_clusters.lineardist$trans.BDP.FINAL)
min(ind_clusters.lineardist$trans.BDA.FINAL)
min(ind_clusters.lineardist$trans.PPFL.FINAL)
min(ind_clusters.lineardist$trans.HD.FINAL)
min(ind_clusters.lineardist$trans.HL.FINAL)
min(ind_clusters.lineardist$trans.SL.FINAL)
min(ind_clusters.lineardist$trans.EW.FINAL)
min(ind_clusters.lineardist$trans.ML.FINAL)

max(ind_clusters.lineardist$trans.CP.FINAL)
max(ind_clusters.lineardist$trans.BDP.FINAL)
max(ind_clusters.lineardist$trans.BDA.FINAL)
max(ind_clusters.lineardist$trans.PPFL.FINAL)
max(ind_clusters.lineardist$trans.HD.FINAL)
max(ind_clusters.lineardist$trans.HL.FINAL)
max(ind_clusters.lineardist$trans.SL.FINAL)
max(ind_clusters.lineardist$trans.EW.FINAL)
max(ind_clusters.lineardist$trans.ML.FINAL)


#==================================================
#T-test on linear distances among different clusters
#==================================================
library(car)
library(dplyr)
library(ggpubr)

#Welch Two Sample t-test for each group for each linear measurement and cluster
#CP
t.test(ind_clusters.lineardist$trans.CP.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")
#BDP
t.test(ind_clusters.lineardist$trans.BDP.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")
#BDA
t.test(ind_clusters.lineardist$trans.BDA.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")
#PPFL
t.test(ind_clusters.lineardist$trans.PPFL.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")
#HD
t.test(ind_clusters.lineardist$trans.HD.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")
#HL
t.test(ind_clusters.lineardist$trans.HL.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")
#SL
t.test(ind_clusters.lineardist$trans.SL.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")
#EW
t.test(ind_clusters.lineardist$trans.EW.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")
#ML
t.test(ind_clusters.lineardist$trans.ML.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")

#boxplot for linear distances
CP <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.CP.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Caudal peduncle depth", x="", y="Value") +theme_bw()+
  theme(plot.title = element_text(size=10, face = "bold"), 
  axis.title.y = element_text(size=9, margin = margin(t = 0, r = 5, b = 0, l = 0)))+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

BDP <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.BDP.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Body depth posterior", x="", y="") +theme_bw()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

BDA <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.BDA.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Body depth anterior", x="", y="") +theme_bw()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

PPFL <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.PPFL.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Postpelvic fin length", x="", y="Value") +theme_bw()+
  theme(plot.title = element_text(size=10, face = "bold"),
  axis.title.y = element_text(size=9, margin = margin(t = 0, r = 5, b = 0, l = 0)))+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

HD <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.HD.FINAL,fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Head depth", x="", y="") +theme_bw()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

HL <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.HL.FINAL,fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Head length", x="", y="") +theme_bw()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

SL <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.SL.FINAL,fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Snout length", x="", y="Value") +theme_bw()+
  theme(plot.title = element_text(size=10, face = "bold"),
        axis.title.y = element_text(size=9, margin = margin(t = 0, r = 5, b = 0, l = 0)))+
theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

EW <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.EW.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Eye width", x="", y="") +theme_bw()+
  theme(plot.title = element_text(size=10, face = "bold"))+ 
        theme(plot.title = element_text(size=10, face="bold"))+
          scale_fill_manual(values = c("#0072B2", "#E69F00"))+
          theme(legend.position="none")+
          scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ML <- ggplot(ind_clusters.lineardist, aes(x=as.factor(cluster), y=trans.ML.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=1, width=0.5) + 
  labs(title="Maxilla length", x="", y="") +theme_bw()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
LM_boxplots <- ggarrange(CP, BDP, BDA, PPFL, HD, HL, SL, EW, ML, 
ncol=3, nrow=3)

ggsave("9LM_boxplots_highres.png",LM_boxplots, height = 12, width = 12, dpi = 600)


#basic statistics for linear measurements
lineardists <-cbind.data.frame(ind_clusters.lineardist, forklength)
colnames(lineardists)[13] <- "FL"
colnames(lineardists)[12] <- "cluster"
colnames(lineardists)[11] <- "year"
lineardists

mean(lineardists[lineardists$cluster == 1, "FL"])
mean(lineardists[lineardists$cluster == 2, "FL"])
mean(lineardists[lineardists$cluster == 1, "trans.CP.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.CP.FINAL"])
mean(lineardists[lineardists$cluster == 1, "trans.BDP.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.BDP.FINAL"])
mean(lineardists[lineardists$cluster == 1, "trans.BDA.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.BDA.FINAL"])
mean(lineardists[lineardists$cluster == 1, "trans.PPFL.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.PPFL.FINAL"])
mean(lineardists[lineardists$cluster == 1, "trans.HD.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.HD.FINAL"])
mean(lineardists[lineardists$cluster == 1, "trans.HL.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.HL.FINAL"])
mean(lineardists[lineardists$cluster == 1, "trans.SL.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.SL.FINAL"])
mean(lineardists[lineardists$cluster == 1, "trans.EW.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.EW.FINAL"])
mean(lineardists[lineardists$cluster == 1, "trans.ML.FINAL"])
mean(lineardists[lineardists$cluster == 2, "trans.ML.FINAL"])

t.test(lineardists[lineardists$cluster == 1, "FL"],lineardists[lineardists$cluster== 2, "FL"])

#REPEAT 
#boxplot for linear distances
CP <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.CP.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="CPD: Caudal peduncle depth", x="", y="Length (cm)") +theme_bw()+
  theme(plot.title = element_text(size=10, face = "bold"), 
        axis.title.y = element_text(size=9, margin = margin(t = 0, r = 5, b = 0, l = 0)))+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("mediumorchid4", "mediumorchid4"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

BDP <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.BDP.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="BDP: Body depth posterior", x="", y="") +theme_bw()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("blue3", "blue3"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

BDA <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.BDA.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="BDA: Body depth anterior", x="", y="") +theme_bw()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("cyan1", "cyan1"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

PPFL <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.PPFL.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="PPFL: Postpelvic fin length", x="", y="Length (cm)") +theme_bw()+
  theme(plot.title = element_text(size=10, face = "bold"),
        axis.title.y = element_text(size=9, margin = margin(t = 0, r = 5, b = 0, l = 0)))+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("dark green", "dark green"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

HD <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.HD.FINAL,fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="", x="", y="") +theme_classic()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("chartreuse2", "chartreuse2"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

HL <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.HL.FINAL,fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="", x="", y="") +theme_classic()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("yellow", "yellow"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

SL <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.SL.FINAL,fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="", x="", y="") +theme_classic()+
  theme(plot.title = element_text(size=10, face = "bold"),
        axis.title.y = element_text(size=9, margin = margin(t = 0, r = 5, b = 0, l = 0)))+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("orange", "orange"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

EW <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.EW.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="", x="", y="") +theme_classic()+
  theme(plot.title = element_text(size=10, face = "bold"))+ 
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("red2", "red2"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

ML <- ggplot(lineardists, aes(x=as.factor(cluster), y=trans.ML.FINAL, fill=cluster)) + 
  geom_boxplot(alpha=0.6, width=0.5) + 
  labs(title="", x="", y="") +theme_classic()+
  theme(plot.title = element_text(size=10, face="bold"))+
  scale_fill_manual(values = c("firebrick4", "firebrick4"))+
  theme(legend.position="none")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))

ggarrange(HD, HL, SL, EW, ML, 
          ncol=3, nrow=2)





#make table for PCA and t-test data
#PCA, Explained variance (%)
PCA_df <- data.frame (PCAxis = c("PC1","PC2","PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9"),
                      PercentVariance = c(38.0, 26.8, 10.3, 8.6, 7.1, 3.9, 2.6, 1.5, 1.3))
PCA_df
write.csv(PCA_df, file = "C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/PCAtable_strictdata.csv")


#table 2 with T-test results
t.test.df <- data.frame("LM" = c("CP", "BDP", "BDA", "PPFL", "HD", "HL", "SL", "EW", "ML"),
                "LM2"= c("Caudal peduncle", "Body depth posterior","Body depth anterior","Postpelvic fin length","Head depth","Head length","Snout length","Eye width", "Maxilla length"),
                "t"= c(-1.3,5.2,16.5,25.9,13.4,16.8,3.4,-2.3,4.9),
                "df" = c(44.7,46.0,45.5,45.6,45.9,47.2,47.6,45.9,48.5),
                "p-value" = c("0.19","p < .001","p < .001","p < .001","p < .001","p < .001","0.0014","0.02","p < .001"))
write.csv(t.test.df, file = "C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/t.test.table_strictdata.csv")

#boxplot of LMs divided by cluster
install.packages("reshape")
library("reshape")

ind_clusters.lineardist_reshapedv1 <- ind_clusters.lineardist[-c(11,12)]
ind_clusters.lineardist_reshapedv2 <- melt(ind_clusters.lineardist_reshapedv1 , id.vars="ID",
                                          measure.cars=c("trans.CP.FINAL","trans.BDP.FINAL","trans.BDA.FINAL","trans.PPFL.FINAL","trans.HD.FINAL", "trans.HL.FINAL","trans.SL.FINAL","trans.EW.FINAL","trans.ML.FINAL"))

colnames(ind_clusters.lineardist)[11] <- "year"
colnames(ind_clusters.lineardist)[12] <- "cluster"

cluster_repeat <-rep(ind_clusters.lineardist$cluster, 9)  

#repeat rows and bind
ind_clusters.lineardist_repeat <- cbind.data.frame(ind_clusters.lineardist_reshapedv2, cluster_repeat)


ggplot(ind_clusters.lineardist_repeat, aes(x=as.factor(variable), y=value, fill=as.factor(cluster_repeat))) +
  geom_boxplot()+
  labs(x="", y="Size-adjusted linear distance (mm)", col="")+
  theme_bw()









#==============================
#Discriminant Function Analysis
#==============================
library(MASS)
library(klaR)

#check classes for each variable
lapply(cluster_indpca3, class)

#check assumptions for LDA
#normality, equal variances, multicollinearity, independence
#Ho: variance-covariance matrix is same across groups
#Ha: variance-covariance matrix is not same across groups
boxM(cluster_indpca3[,2:6],cluster_indpca3$cluster)


#perform LDA
pca.lda <- lda(cluster ~ Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5, data=cluster_indpca3)
pca.lda

pca.lda.p <- predict(pca.lda, cluster_indpca3)  #get coords for LDA plot
pca.lda.p

#histogram base on LD1
ldahist(data=pca.lda.p$x[,1], g=cluster_indpca3$cluster)

#partition plot - classification of each Dimension bw clusters
partimat(cluster~Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5, data = cluster_indpca3, method = "lda")


#LDA for Head depth
pca.lda.hd <- lda(cluster~trans.HD.FINAL, data=ind_clusters.lineardist_renamed)
pca.lda.hd.p <- predict(pca.lda.hd, ind_clusters.lineardist_renamed)  #get coords for LDA plot
ldahist(data=pca.lda.hd.p$x[,1], g=ind_clusters.lineardist_renamed$cluster)

#========================
#Stable Isotope Analysis 
#========================
install.packages("expression")
library("expression")

#2018&2019 data
# data: clusters_strictdata_isotopes.csv 

SI<- read.csv("C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/clusters_strictdata_isotopes.csv")
head(SI)

#Explore data and normality
par(mfrow=c(2,2))
hist(SI$d13C_PLAS)
hist(SI$dN15_PLAS)
hist(SI$d13C_RBC)
hist(SI$dN15_RBC)

mean(SI$dN15_RBC, na.rm=TRUE)


#apply diet discrimination factor
#all 13C and 15N data going forward has been corrected
SI$d13C_PLAS <- SI$d13C_PLAS  - 1.2 #C PLAS
SI$d13C_PLAS
min(SI$d13C_PLAS,na.rm = T)
max(SI$d13C_PLAS,na.rm = T)

SI$d13C_RBC <- SI$d13C_RBC - 0.1 #C RBC
SI$d13C_RBC
min(SI$d13C_RBC,na.rm = T)
max(SI$d13C_RBC,na.rm = T)

SI$dN15_PLAS <- SI$dN15_PLAS - 0.9 #N PLAS
SI$dN15_PLAS
min(SI$dN15_PLAS,na.rm = T)
max(SI$dN15_PLAS,na.rm = T)

SI$dN15_RBC <- SI$dN15_RBC - 1.1  #N RBC
SI$dN15_RBC 
min(SI$dN15_RBC,na.rm = T)
max(SI$dN15_RBC,na.rm = T)

SI$X13c.value1= SI$X13c.value1  - 1.2 #C PLAS
SI$X13c.value1
SI$X13c.value2= SI$X13c.value2 - 0.1 #C RBC
SI$X13c.value2
SI$X15n.value1 = SI$X15n.value1 - 0.9 #N PLAS
SI$X15n.value1
SI$X15n.value2 = SI$X15n.value2 - 1.1  #N RBC
SI$X15n.value2
SI

#=====================
# d13C vs d15N PLASMA
#=====================
#resource for interpreting regression output: https://towardsdatascience.com/understanding-linear-regression-output-in-r-7a9cbda948b3


#Linear Regression PLASMA, cluster 1
fitPLASMA1 <- function(fitPLASMA){ 
  m <- lm(dN15_PLAS ~ d13C_PLAS, data=subset(SI, cluster_2018.2019==1));
    eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#Linear Regression PLASMA, cluster 2
fitPLASMA2 <- function(fitPLASMA){ 
  m <- lm(dN15_PLAS ~ d13C_PLAS, data=subset(SI, cluster_2018.2019==2));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#plot of plasma data
ggplot(data= SI, aes (x=d13C_PLAS, y=dN15_PLAS)) +
       geom_point(aes(color=factor(cluster_2018.2019)))+
  labs(title="PLASMA",
      x=delta^13~"C%",
       y= delta^15~"N%", 
       col="Cluster")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE,aes(color = factor(cluster_2018.2019)))+
    theme_bw()+
  geom_text(x = -18.8, y = 15.5, label = fitPLASMA1(SI), parse = TRUE, color="#F8766D")+   #red
  geom_text(x = -18.8, y = 15.4, label = fitPLASMA2(SI), parse = TRUE, color="#00BFC4")    #blue


#plot to verify line & equation
dev.off() #turn off par function
plot(dN15_PLAS ~ d13C_PLAS, data=subset(SI, cluster_2018.2019==1))
lm(dN15_PLAS ~ d13C_PLAS, data=subset(SI, cluster_2018.2019==1))
abline(29.2497, 0.6422)

plot(dN15_PLAS ~ d13C_PLAS, data=subset(SI, cluster_2018.2019==2))
lm(dN15_PLAS ~ d13C_PLAS, data=subset(SI, cluster_2018.2019==2))
abline(25.0081 ,  0.4416 )

#=====================
# d13C vs d15N RBC
#=====================
#Linear Regression RBC, cluster 1
fitRBC1 <- function(fitRBC){ 
  m <- lm(dN15_RBC ~ d13C_RBC, data=subset(SI, cluster_2018.2019==1));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#Linear Regression RBC, cluster 2
fitRBC2 <- function(fitRBC){ 
  m <- lm(dN15_RBC ~ d13C_RBC, data=subset(SI, cluster_2018.2019==2));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#plot of rbc data
ggplot(data= SI, aes (x=d13C_RBC, y=dN15_RBC)) +
  geom_point(aes(color=factor(cluster_2018.2019)))+
  labs(title="RBC",
      x=delta^13~"C%",
       y= delta^15~"N%", 
       col="Cluster")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE,aes(color = factor(cluster_2018.2019)))+
  theme_bw()+
  geom_text(x = -19.5, y = 15.6, label = fitRBC1(SI), parse = TRUE, color="#F8766D")+   #red
  geom_text(x = -19.5, y = 15.5, label = fitRBC2(SI), parse = TRUE, color="#00BFC4")    #blue



#plot to verify line & equation
dev.off() #turn off par function
plot(dN15_RBC ~ d13C_RBC, data=subset(SI, cluster_2018.2019==1))
lm(dN15_RBC ~ d13C_RBC, data=subset(SI, cluster_2018.2019==1))
abline(20.8518,  0.2767 )

plot(dN15_RBC ~ d13C_RBC, data=subset(SI, cluster_2018.2019==2))
lm(dN15_RBC ~ d13C_RBC, data=subset(SI, cluster_2018.2019==2))
abline(21.6508, 0.3018)

SI_edited<- read.csv("C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/clusters_strictdata_isotopes_edited_discriminated.csv")
head(SI_edited)

SI_edited_C<- subset(SI_edited[1:90, 64:66]) #use for boxplot only!
na.omit(SI_edited_C) 
as.numeric(SI_edited_C$X13c_value_lediscrim)

SI_edited_N <- subset(SI_edited[1:90, 66:68]) #use for boxplot only!
na.omit(SI_edited_N) 
as.numeric(SI_edited_C$X15n_value_lediscrim)

#boxplot comparisons for 13C and 15N for each tissue
boxplot_13C <- ggplot(SI_edited_C, aes(x=as.factor(X13c_tissue_lediscrim), y=as.numeric(X13c_value_lediscrim), fill=as.factor(cluster))) +
  geom_boxplot()+
  labs(x="", y=delta^13~"C%", col="")+
  theme_bw()+
  theme(plot.title = element_text(size=9, face = "bold"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)+
  scale_fill_manual(labels=c("Cluster 1", "Cluster 2"), values = c("#4682B4", "#B47846"))+
  scale_x_discrete(labels= c("Plasma", "RBC"))+
  theme(legend.title = element_blank(),
        legend.direction = "horizontal")

boxplot_15N <- ggplot(SI_edited_N, aes(x=as.factor(X15n_tissue_lediscrim), y=as.numeric(X15n_value_lediscrim), fill=as.factor(cluster))) +
  geom_boxplot()+
  labs(x="", y=delta^15~"N%", col="")+
  theme_bw()+
  theme(plot.title = element_text(size=9, face = "bold"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)+
  scale_fill_manual(labels=c("Cluster 1", "Cluster 2"), values = c("#4682B4", "#B47846"))+
  scale_x_discrete(labels= c("RBC", "Plasma"))+
  theme(legend.title = element_blank(),
        legend.direction = "horizontal")

library(ggpubr)
ggarrange(
  boxplot_13C, boxplot_15N, labels = c("A", "B"),
  common.legend = TRUE, legend = "bottom"
)

#boxplots, REARRANGED
#boxplot comparisons for 13C and 15N for each tissue (assume all data is here, 2018&2019)
boxplot_13C2 <- ggplot(SI_edited_C, aes(x=as.factor(cluster), y=as.numeric(X13c_value), fill=as.factor(X13c_tissue))) +
  geom_boxplot()+
  labs(x="", y=delta^13~"C%", col="")+
  theme_bw()+
  theme(plot.title = element_text(size=9, face = "bold"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)+
  scale_fill_manual(labels=c("Plasma", "RBC"), values = c("#4682B4", "#B47846"),name="Tissue type")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(legend.title = element_text(),
        legend.direction = "horizontal")

boxplot_15n2 <- ggplot(SI_edited_N, aes(x=as.factor(cluster), y=as.numeric(X15n_value), fill=as.factor(X15n_tissue))) +
  geom_boxplot()+
  labs(x="", y=delta^15~"N%", col="")+
  theme_bw()+
  theme(plot.title = element_text(size=9, face = "bold"))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)+
  scale_fill_manual(labels=c("Plasma", "RBC"), values = c("#4682B4", "#B47846"), name="Tissue type  ")+
  scale_x_discrete(labels= c("Cluster 1", "Cluster 2"))+
  theme(legend.title = element_text(),
        legend.direction = "horizontal")
       

ggarrange(
  boxplot_13C2, boxplot_15n2, labels = c("", ""),
  common.legend = TRUE, legend = "bottom"
)

#SI difference data
SI_difference_raw <- read.csv("C:/Users/Stephanie/Desktop/UNBC/Cod Data/Morphometrics/RStudio/clusters_strictdata_isotopes_edited_discriminated - COPY.csv")
head(SI_difference_raw)


#get difference of RBC-plasma for C and N
C_difference <- as.data.frame.vector(SI_difference_raw$C_RBC - SI_difference_raw$C_plasma)
C_difference
N_difference <- as.data.frame.vector(SI_difference_raw$N_RBC - SI_difference_raw$N_plasma)
N_difference

SI_difference <- cbind(SI_difference_raw, C_difference, N_difference)
head(SI_difference)

colnames(SI_difference)[5] <- "cluster"
colnames(SI_difference)[65] <- "C_difference"
colnames(SI_difference)[66] <- "N_difference"
SI_difference


#summarize data, Mean and standard errors
install.packages("plotrix")
library("plotrix")

#fork length basic stats
min(SI_difference[SI_difference$cluster ==2, "FL_mm"], na.rm=TRUE)
max(SI_difference[SI_difference$cluster ==2, "FL_mm"], na.rm=TRUE)
length(SI_difference[SI_difference$cluster ==2, "FL_mm"])
mean(SI_difference[SI_difference$cluster ==2, "FL_mm"], na.rm=TRUE)
sd(SI_difference[SI_difference$cluster ==2, "FL_mm"], na.rm=TRUE)

#stable isotope basic stats
mean(SI_difference[SI_difference$cluster ==2, "C_RBC"], na.rm=TRUE)
mean(SI_difference[SI_difference$cluster ==2, "C_plasma"], na.rm=TRUE)
mean(SI_difference[SI_difference$cluster ==2, "N_RBC"], na.rm=TRUE)
mean(SI_difference[SI_difference$cluster ==2, "N_plasma"], na.rm=TRUE)

std.error(SI_difference[SI_difference$cluster ==2, "C_RBC"], na.rm=TRUE)
std.error(SI_difference[SI_difference$cluster ==2, "C_plasma"], na.rm=TRUE)
std.error(SI_difference[SI_difference$cluster ==2, "N_RBC"], na.rm=TRUE)
std.error(SI_difference[SI_difference$cluster ==2, "N_plasma"], na.rm=TRUE)

min(SI_difference[SI_difference$cluster ==2, "C_RBC"], na.rm=TRUE)
min(SI_difference[SI_difference$cluster ==2, "C_plasma"], na.rm=TRUE)
min(SI_difference[SI_difference$cluster ==2, "N_RBC"], na.rm=TRUE)
min(SI_difference[SI_difference$cluster ==2, "N_plasma"], na.rm=TRUE)

max(SI_difference[SI_difference$cluster ==2, "C_RBC"], na.rm=TRUE)
max(SI_difference[SI_difference$cluster ==2, "C_plasma"], na.rm=TRUE)
max(SI_difference[SI_difference$cluster ==2, "N_RBC"], na.rm=TRUE)
max(SI_difference[SI_difference$cluster ==2, "N_plasma"], na.rm=TRUE)

sd(SI_difference[SI_difference$cluster ==2, "C_RBC"], na.rm=TRUE)
sd(SI_difference[SI_difference$cluster ==2, "C_plasma"], na.rm=TRUE)
sd(SI_difference[SI_difference$cluster ==2, "N_RBC"], na.rm=TRUE)
sd(SI_difference[SI_difference$cluster ==2, "N_plasma"], na.rm=TRUE)

length(SI_difference[SI_difference$cluster ==2, "C_RBC"])
length(SI_difference[SI_difference$cluster ==2, "C_plasma"])
length(SI_difference[SI_difference$cluster ==2, "N_RBC"])
length(SI_difference[SI_difference$cluster ==2, "N_plasma"])

length(SI_difference[SI_difference$cluster ==2, "C_difference"])
length(SI_difference[SI_difference$cluster ==2, "N_difference"])



# >0.05, distribution of data not significantly different from normal. = assume normality
shapiro.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 1, "C_RBC"])
shapiro.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 2, "C_RBC"])
shapiro.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 1, "C_plasma"])
shapiro.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 2, "C_plasma"])
shapiro.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 1, "N_RBC"])
shapiro.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 2, "N_RBC"])
shapiro.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 1, "N_plasma"])
shapiro.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 2, "N_plasma"])

t.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 1, "C_RBC"],SI_difference[SI_difference$C_RBC & SI_difference$cluster == 2, "C_RBC"])
t.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 1, "C_plasma"],SI_difference[SI_difference$C_RBC & SI_difference$cluster == 2, "C_plasma"])
t.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 1, "N_RBC"],SI_difference[SI_difference$C_RBC & SI_difference$cluster == 2, "N_RBC"])
t.test(SI_difference[SI_difference$C_RBC & SI_difference$cluster == 1, "N_plasma"],SI_difference[SI_difference$C_RBC & SI_difference$cluster == 2, "N_plasma"])

SI_difference_subset <- data.frame(SI_difference$ID,SI_difference$cluster,SI_difference$C_difference, SI_difference$N_difference)
SI_difference_subset

#equal variance test
SI_C_PLAS <- SI_edited_C[c(1:15, 33:49) , 2:3] #PLAS, C
SI_C_RBC <- SI_edited_C[c(16:32, 50:68), 2:3] #RBC, C
SI_N_PLAS <- SI_edited_N[c(1:15, 33:49), -2] #PLAS, N
SI_N_RBC <- SI_edited_N[c(16:32, 50:68), -2] #RBC, N

leveneTest(X13c_value_lediscrim ~ as.factor(cluster), data=SI_C_PLAS) #C
leveneTest(X13c_value_lediscrim ~ as.factor(cluster), data=SI_C_RBC)
leveneTest(X15n_value_lediscrim ~ as.factor(cluster), data=SI_N_PLAS) #N
leveneTest(X15n_value_lediscrim ~ as.factor(cluster), data=SI_N_RBC)    #all values >0.05, assumption of equal variances met
      

#t-test on plasma(C1) vs plasma(C2) and RBC(C1) vs RBC(C2)
t.test(SI_edited_C[SI_edited_C$X13c_tissue_lediscrim=="d13C_PLAS" & SI_edited_C$cluster == 1, "X13c_value_lediscrim"],SI_edited_C[SI_edited_C$X13c_tissue_lediscrim=="d13C_PLAS" & SI_edited_C$cluster == 2, "X13c_value_lediscrim"])
t.test(SI_edited_C[SI_edited_C$X13c_tissue_lediscrim=="d13C_RBC" & SI_edited_C$cluster == 1, "X13c_value_lediscrim"],SI_edited_C[SI_edited_C$X13c_tissue_lediscrim=="d13C_RBC" & SI_edited_C$cluster == 2, "X13c_value_lediscrim"])

t.test(SI_edited_N[SI_edited_N$X15n_tissue_lediscrim=="d15N_PLAS" & SI_edited_N$cluster == 1, "X15n_value_lediscrim"],SI_edited_N[SI_edited_N$X15n_tissue_lediscrim=="d15N_PLAS" & SI_edited_N$cluster == 2, "X15n_value_lediscrim"])
t.test(SI_edited_N[SI_edited_N$X15n_tissue_lediscrim=="d15N_RBC" & SI_edited_N$cluster == 1, "X15n_value_lediscrim"],SI_edited_N[SI_edited_N$X15n_tissue_lediscrim=="d15N_RBC" & SI_edited_N$cluster == 2, "X15n_value_lediscrim"])

#t-test on RBC vs PLASMA for each cluster
t.test(SI_edited_C[SI_edited_C$X13c_tissue_lediscrim=="d13C_PLAS" & SI_edited_C$cluster == 1, "X13c_value_lediscrim"],SI_edited_C[SI_edited_C$X13c_tissue_lediscrim=="d13C_RBC" & SI_edited_C$cluster == 1, "X13c_value_lediscrim"])
t.test(SI_edited_C[SI_edited_C$X13c_tissue_lediscrim=="d13C_PLAS" & SI_edited_C$cluster == 2, "X13c_value_lediscrim"],SI_edited_C[SI_edited_C$X13c_tissue_lediscrim=="d13C_RBC" & SI_edited_C$cluster == 2, "X13c_value_lediscrim"])

t.test(SI_edited_N[SI_edited_N$X15n_tissue_lediscrim=="d15N_PLAS" & SI_edited_N$cluster == 1, "X15n_value_lediscrim"],SI_edited_N[SI_edited_N$X15n_tissue_lediscrim=="d15N_RBC" & SI_edited_N$cluster == 1, "X15n_value_lediscrim"])
t.test(SI_edited_N[SI_edited_N$X15n_tissue_lediscrim=="d15N_PLAS" & SI_edited_N$cluster == 2, "X15n_value_lediscrim"],SI_edited_N[SI_edited_N$X15n_tissue_lediscrim=="d15N_RBC" & SI_edited_N$cluster == 2, "X15n_value_lediscrim"])



#multi tissue
#get fork length data
forklength <-as.data.frame(dat[,11])
forklength <-na.omit(forklength)
head(forklength)



#C
#Linear Regression multi-tissue carbon, cluster 1
fitCarbon_1 <- function(fitCarbon){ 
  m <- lm(C_difference ~ FL, data=subset(SI_difference, cluster==1));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#Linear Regression multi-tissue carbon, cluster 2
fitCarbon_2 <- function(fitCarbon){ 
  m <- lm(C_difference ~ FL, data=subset(SI_difference, cluster==2));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                            list(a = format(unname(coef(m)[1]), digits = 3),
                                 b = format(unname(coef(m)[2]), digits = 3),
                                 r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#plot to verify line & equation
dev.off() #turn off par function
plot(C_difference ~ FL, data=subset(SI_difference, cluster==1))
model_Carbon_1 <- lm(C_difference ~ FL, data=subset(SI_difference, cluster==1))
str(summary(model_Carbon_1))  #get R^2 and P-values for carbon_c1 linear model
pVal_model_C_1 <- anova(model_Carbon_1)$'Pr(>F)'[1]

#plot to verify line & equation
dev.off() #turn off par function
plot(C_difference ~ FL, data=subset(SI_difference, cluster==2))
model_Carbon_2 <- lm(C_difference ~ FL, data=subset(SI_difference, cluster==2))
str(summary(model_Carbon_2))  #get R^2 and P-values for carbon_c2 linear model
pVal_model_C_2 <- anova(model_Carbon_2)$'Pr(>F)'[1]

#plot of multi tissue carbon vs FL data, linear regression
C_difference_plot <- ggplot(data= SI_difference, aes (x=FL, y=C_difference)) +
  geom_point(aes(color=factor(cluster)))+
  scale_color_manual(values = c("1" = "#4682B4", "2" = "#B47846"))+
  labs(x="Fork length (mm)",
       y= "RBC - Plasma (???)", 
       col="Cluster")+
  theme_bw()+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE, aes(color = factor(cluster)))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  geom_text(x = 30, y = -2.7, label = "", parse = TRUE, color="black")+   #cluster 1
  geom_text(x = 30, y = -2.9, label = "", parse = TRUE, color="#999999")    #cluster 2


#N
#Linear Regression multi-tissue nitrogen, cluster 1
fitNitrogen_1 <- function(fitNitrogen){ 
  m <- lm(N_difference ~ FL, data=subset(SI_difference, cluster==1));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}


#Linear Regression multi-tissue nitrogen, cluster 2
fitNitrogen_2 <- function(fitNitrogen){ 
  m <- lm(N_difference ~ FL, data=subset(SI_difference, cluster==2));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#Residual plots for multi-tissues
#Carbon
SI_difference_na <- na.omit(SI_difference)
nrow(SI_difference_na)
C_difference_lm <- lm(C_difference~FL, data=SI_difference_na)

#t.test for difference bw groups - multi tissue
t.test(SI_difference_na[SI_difference_na$cluster == 1, "C_difference"],SI_difference_na[SI_difference_na$cluster == 2, "C_difference"])
t.test(SI_difference_na[SI_difference_na$cluster == 1, "N_difference"],SI_difference_na[SI_difference_na$cluster == 2, "N_difference"])




C_residuals <- resid(C_difference_lm)
length(C_residuals) #31
plot(C_residuals)
SI_carbon_dataframe = cbind.data.frame(SI_difference_na, C_residuals)
plot(C_residuals ~ FL, data=SI_carbon_dataframe)
C_FL <- SI_carbon_dataframe$FL
lm(C_residuals ~ C_FL)
abline(-4.005e-16,1.217e-17)  #for cluster 1 & 2 combined. subset to get individuals

#C_cluster 1
SI_carbon_dataframe_c1<- subset(SI_carbon_dataframe, cluster ==1)
SI_carbon_c1_plot <- ggplot(SI_carbon_dataframe_c1, aes(x=FL, y=C_residuals)) +  #nicer plot
  geom_point(color="#4682B4")+
  geom_smooth(method=lm, se=FALSE, color ="#4682B4")+
  theme_bw()+
  labs(x="Fork length (mm)",
       y= "Residuals of diet switch (RBC - Plasma (???))")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))

#C_cluster 2
SI_carbon_dataframe_c2<- subset(SI_carbon_dataframe, cluster ==2)
SI_carbon_c2_plot <- ggplot(SI_carbon_dataframe_c2, aes(x=FL, y=C_residuals)) +  #nicer plot
  geom_point(color="#B47846")+
  geom_smooth(method=lm, se=FALSE, color = "#B47846")+
  theme_bw()+
  labs(x="Fork length (mm)",
       y= "")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))

ggarrange(
  SI_carbon_c1_plot, SI_carbon_c2_plot, labels = c("", "B"),
  common.legend = TRUE, legend = "right"
)



#Nitrogen
SI_difference_na <- na.omit(SI_difference)
nrow(SI_difference_na)
N_difference_lm <- lm(N_difference~FL, data=SI_difference_na)

N_residuals <- resid(N_difference_lm)
length(N_residuals) #31
plot(N_residuals)
SI_nitrogen_dataframe = cbind.data.frame(SI_difference_na, N_residuals)
plot(N_residuals ~ FL, data=SI_nitrogen_dataframe)
N_FL <- SI_nitrogen_dataframe$FL
lm(N_residuals ~ N_FL)
abline(-6.582e-17 ,1.751e-18)  #for cluster 1 & 2 combined. subset to get individuals

#N_cluster 1
SI_nitrogen_dataframe_c1<- subset(SI_nitrogen_dataframe, cluster ==1)
SI_nitrogen_c1_plot <- ggplot(SI_nitrogen_dataframe_c1, aes(x=FL, y=N_residuals)) +  #nicer plot
  geom_point(color="#4682B4")+
  geom_smooth(method=lm, se=FALSE, color ="#4682B4")+
  theme_bw()+
  labs(x="Fork length (mm)",
       y= "Residuals of diet switch (RBC - Plasma (???))")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))

#N_cluster 2
SI_nitrogen_dataframe_c2<- subset(SI_nitrogen_dataframe, cluster ==2)
SI_nitrogen_c2_plot <- ggplot(SI_nitrogen_dataframe_c2, aes(x=FL, y=N_residuals)) +  #nicer plot
  geom_point(color="#B47846")+
  geom_smooth(method=lm, se=FALSE, color = "#B47846")+
  theme_bw()+
  labs(x="Fork length (mm)",
       y= "")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))

ggarrange(
  SI_nitrogen_c1_plot, SI_nitrogen_c2_plot, labels = c("A", "B"),
  common.legend = TRUE, legend = "right"
)


#plot to verify line & equation
dev.off() #turn off par function
plot(N_difference ~ FL, data=subset(SI_difference, cluster==1))
model_Nitrogen_1 <- lm(N_difference ~ FL, data=subset(SI_difference, cluster==1))
str(summary(model_Nitrogen_1))  #get R^2 and P-values for carbon_c1 linear model
pVal_model_N_1 <- anova(model_Nitrogen_1)$'Pr(>F)'[1]

#plot to verify line & equation
dev.off() #turn off par function
plot(N_difference ~ FL, data=subset(SI_difference, cluster==2))
model_Nitrogen_2 <- lm(N_difference ~ FL, data=subset(SI_difference, cluster==2))
str(summary(model_Nitrogen_2))  #get R^2 and P-values for carbon_c2 linear model
pVal_model_N_2 <- anova(model_Nitrogen_2)$'Pr(>F)'[1]

#plot of multi-tissue nitrogen data, linear regression
N_difference_plot <- ggplot(data= SI_difference, aes (x=FL, y=N_difference)) +
  geom_point(aes(color=factor(cluster)))+
  scale_color_manual(values = c("1" = "#4682B4", "2" = "#B47846"))+
  labs(x="Fork length (mm)",
       y= "", 
       col="Cluster")+
  theme(plot.title = element_text(size=12, face="bold")) +
  geom_smooth(method = "lm", se = FALSE, aes(color = factor(cluster)))+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  geom_text(x = 30, y = -2.7, label = "", parse = TRUE, color="black")+   #cluster 1
  geom_text(x = 30, y = -2.9, label = "", parse = TRUE, color="#999999")    #cluster 2
  
ggarrange(
  C_difference_plot, N_difference_plot, labels = c("A", "B"),
  common.legend = TRUE, legend = "right"
)


#C
#Linear Regression 13C , cluster 1
fitCarbon_1 <- function(fitCarbon){ 
  m <- lm(C_difference ~ FL, data=subset(SI_difference, cluster==1));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#Linear Regression multi-tissue carbon, cluster 1
fitCarbon_2 <- function(fitCarbon){ 
  m <- lm(C_difference ~ FL, data=subset(SI_difference, cluster==2));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#plot to verify line & equation
dev.off() #turn off par function
plot(C_difference ~ FL, data=subset(SI_difference, cluster==1))
model_Carbon_1 <- lm(C_difference ~ FL, data=subset(SI_difference, cluster==1))
str(summary(model_Carbon_1))  #get R^2 and P-values for carbon_c1 linear model
pVal_model_C_1 <- anova(model_Carbon_1)$'Pr(>F)'[1]

#plot to verify line & equation
dev.off() #turn off par function
plot(C_difference ~ FL, data=subset(SI_difference, cluster==2))
model_Carbon_2 <- lm(C_difference ~ FL, data=subset(SI_difference, cluster==2))
str(summary(model_Carbon_2))  #get R^2 and P-values for carbon_c2 linear model
pVal_model_C_2 <- anova(model_Carbon_2)$'Pr(>F)'[1]

#plot of multi tissue carbon vs FL data, linear regression
C_difference_plot <- ggplot(data= SI_difference, aes (x=FL, y=C_difference)) +
  geom_point(aes(color=factor(cluster)))+
  scale_color_manual(values = c("1" = "#4682B4", "2" = "#B47846"))+
  labs(x="Fork length (mm)",
       y= "RBC - Plasma (???)", 
       col="Cluster")+
  theme_bw()+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE, aes(color = factor(cluster)))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  geom_text(x = 30, y = -2.7, label = "", parse = TRUE, color="black")+   #cluster 1
  geom_text(x = 30, y = -2.9, label = "", parse = TRUE, color="#999999")    #cluster 2



#Diet switch (RBC-plasma) for C and N 
boxplot_13c2 <- ggplot(SI_difference, aes(x=as.factor(cluster), y=C_difference))+
    geom_boxplot(fill=c("#276695", "#B47B08"))+
    labs(x="", y="Habitat switch (RBC-plasma)", col="")+
    theme_bw()+
        theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),)+
    scale_x_discrete(labels=c("Cluster 1", "Cluster 2")) +
    theme(legend.title = element_blank(),
          legend.direction = "horizontal")+
  theme(plot.margin = margin(.3,.3,0,.3,"cm"))
    
boxplot_15n2 <- ggplot(SI_difference, aes(x=as.factor(cluster), y=N_difference)) +
  geom_boxplot(fill=c("#276695", "#B47B08"))+
  labs(x="", y="Trophic switch (RBC-plasma)", col="")+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)+
  scale_x_discrete(labels=c("Cluster 1", "Cluster 2")) +
  theme(legend.title = element_blank(),
        legend.direction = "horizontal")+
  theme(plot.margin = margin(.3,.3,0,.3,"cm"))

ggarrange(boxplot_13c2, boxplot_15n2, common.legend = TRUE, legend = "bottom")

#t.test of differences bw groups
t.test(SI_difference[SI_difference$cluster == 1, "C_difference"],SI_difference[SI_difference$cluster == 2, "C_difference"])
t.test(SI_difference[SI_difference$cluster == 1, "N_difference"],SI_difference[SI_difference$cluster == 2, "N_difference"])
length(SI_difference[SI_difference$cluster == 2, "C_difference"])

#mean for multi tissue
SI_difference_na <- na.omit(SI_difference)
C_mtissue_na <- data.frame(SI_difference_na$cluster, SI_difference_na$C_difference)
mean(C_mtissue_na [C_mtissue_na$SI_difference_na.cluster == 1, "SI_difference_na.C_difference"])
mean(C_mtissue_na [C_mtissue_na$SI_difference_na.cluster == 2, "SI_difference_na.C_difference"])

N_mtissue_na <- data.frame(SI_difference_na$cluster, SI_difference_na$N_difference)
mean(N_mtissue_na [N_mtissue_na$SI_difference_na.cluster == 1, "SI_difference_na.N_difference"])
mean(N_mtissue_na [N_mtissue_na$SI_difference_na.cluster == 2, "SI_difference_na.N_difference"])

       
#anova - variation bw groups
SI_difference_year <- cbind.data.frame(SI_difference, year)

anova_C_difference <- aov(C_difference ~ cluster + year, data = SI_difference) # does cluster or yr impact diet switch values for 13C?
anova_N_difference <- aov(N_difference ~ cluster + year, data = SI_difference) # does cluster or yr impact diet switch values for 15N?
summary(anova_C_difference)
summary(anova_N_difference)

#remove outliers
#no outliers in 13c

#anova for N_difference no outlier
anova_N_difference_nooutlier <- aov(N_difference ~ cluster, data = N_difference_nooutlier)
summary(anova_N_difference_nooutlier)
#remove outliers from N and run anova
colnames(N_difference)[1] <- "SI_difference"
max(N_difference$SI_difference, na.rm=TRUE) #outliers
min(N_difference$SI_difference, na.rm=TRUE)

N_difference
N_difference_nooutlier <- SI_difference[-c(12,38),]
N_difference_nooutlier

#plot
boxplot_15n_nooutlier<- ggplot(N_difference_nooutlier, aes(x=as.factor(cluster), y=N_difference)) +
  geom_boxplot(fill=c("#276695", "#B47B08"))+
  labs(x="", y="", col="")+
  theme_bw()+
  theme(plot.title = element_text(size=9, face = "bold"),
        axis.text.x = element_text(size=14, face="bold", vjust=-2),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14, vjust=4))+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)+
  scale_x_discrete(labels=c("Cluster 1", "Cluster 2")) +
  theme(legend.title = element_blank(),
        legend.direction = "horizontal")+
  theme(plot.margin = margin(.3,.3,0,.3,"cm"))

#anova for N_difference no outlier
anova_N_difference_nooutlier <- aov(N_difference ~ cluster, data = N_difference_nooutlier)
summary(anova_N_difference_nooutlier)




#boxplot
library(ggplot2)
x <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
y<- c( 2, 2, 3, 3, 3, 3, 4, 5, 5, 5, 5, 6, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5, 5, 6)
example <- cbind.data.frame(x, y)
examplebox <- ggplot (example, aes(x=as.factor(x), y=y, fill=x))+
  geom_boxplot(alpha=0.3, width=0.3, fill=c("grey"))+
  theme_classic()+
  labs(x="", y="Score", col="")+
  theme(legend.position="none")

#===========================  
#GLM to calculate IS
#===========================
library(lme4)
SI #test for normality above w shapiro-wilk
is.factor(SI$X) #DELETE START
is.numeric(SI$X)

#SI_difference <- subset(SI_difference[(1:59),])
SI_difference
year <-rep(c("2018", "2019"),times=c(19,26))
SI_difference <- cbind.data.frame(SI_difference, year)


SI_difference= na.omit(SI_difference)
nrow(SI_difference) #45

SI_difference
#LM
qqnorm(SI_difference$C_plasma, pch = 1, frame = FALSE)
qqline(SI_difference$C_plasma, col = "steelblue", lwd = 2)

SI
model_C<- lm (formula = d13C_PLAS ~ d13C_RBC, data = SI)
anova(model_C)
plot(model_C)
summary(model_C)

#BIC (mean sum of sq)
BIC_C = (2*11.176)/(2*(31-1))

#WIC mean sum of squares of the error
#sum of sq of error = sum of sq of residual error/deg free
WIC_C = 43.914/(31*(2-1))
BIC_C /WIC_C 

#WIC for each individual (Var)
SI$d13C_PLAS <- na.omit(SI$d13C_PLAS)
mean(SI$d13C_PLAS, na.rm = TRUE)
prob_C_PLAS <- dnorm(SI$d13C_PLAS, mean=-22.105, sd=1.333117) #PLASMA ind variance
E_C_PLAS = weighted.mean(x, p)
E_C_PLAS = weighted.mean(SI$d13C_PLAS, prob_C_PLAS, na.rm=TRUE)
E_C_PLAS = abs(E_C_PLAS)
ind_var_C_PLAS <- E*((x-mean(x))^2)
ind_var_C_PLAS <- E_C_PLAS*((SI$d13C_PLAS - mean(SI$d13C_PLAS, na.rm = TRUE))^2)

mean(SI$d13C_RBC, na.rm=TRUE)
sd(SI$d13C_RBC, na.rm=TRUE)
prob_C_RBC <- dnorm(SI$d13C_RBC, mean=-20.90639, sd=0.6883984) #RBC ind variance
E_C_RBC = weighted.mean(x, p)
E_C_RBC = weighted.mean(SI$d13C_RBC, prob_C_RBC, na.rm = TRUE)
E_C_RBC = abs(E_C_RBC)
ind_var_C_RBC <- E*((x-mean(x))^2)
ind_var_C_RBC <- E_C_RBC*((SI$d13C_RBC - mean(SI$d13C_RBC, na.rm=TRUE))^2)


plot(SI$d13C_PLAS,SI$d13C_RBC)
#add variaces from dependent variables
#Var(X+Y)=Var(X)+Var(Y)+2Cov(X,Y)
d13C_PLAS_na<-na.omit(SI$d13C_PLAS)
d13C_RBC_na<-na.omit(SI$d13C_RBC)
ind_var_C_RBC_na<-na.omit(ind_var_C_RBC)
ind_var_C_PLAS_na<-na.omit(ind_var_C_PLAS)

ind_var_C <- ind_var_C_PLAS + ind_var_C_RBC + 2*(cov(SI$d13C_PLAS, SI$d13C_RBC, use = "pairwise.complete.obs",method = "pearson"))

#IS for each individual
BIC/WIC(ind)
IS_C_ind <-BIC_C/ind_var_C
hist(IS_C_ind)  #check normality (response var in linear regression)
shapiro.test(IS_C_ind) # p<.05 - not normal

logtrans_indIS_C<- log10(IS_C_ind)
hist(logtrans_indIS_C)
shapiro.test(logtrans_indIS_C)

#multiple linear regression - multi tissue
model_C_tissue <- lm (formula = ind_var_C ~ cluster_2018.2019 + Year,data = SI)
anova(model_C_mtissue)
summary(model_C_mtissue)
plot(model_C)


#IS from package
library("RInSp")
SI_C <- subset(SI_difference[,c(2,5,56,62)])
SI_C$C_plasma = abs(SI_C$C_plasma)
SI_C$C_RBC = abs(SI_C$C_RBC)
SI_C

SI_C_C1 <- subset(SI_C, cluster==1)
SI_C_C2 <- subset(SI_C, cluster==2)

RIS_C <- import.RInSp(SI_C, col.header=TRUE, row.names = 1, info.cols= 1:2,
                    subset.column = 0, subset.rows = NA, data.type= "double")


RIS_C_values <- WTcMC(RIS_C, replicates = 999, weight = "equal", print.ris=TRUE)
sumMC.RInSp(RIS_C_values)

#C1
RIS_C_C1 <- import.RInSp(SI_C_C1, col.header=TRUE, row.names = 1, info.cols= 1:2,
                      subset.column = 0, subset.rows = NA, data.type= "double")

RIS_C_C1_values <- WTcMC(RIS_C_C1, replicates = 999, weight = "equal", print.ris=TRUE)
sumMC.RInSp(RIS_C_C1_values)

#C2
RIS_C_C2 <- import.RInSp(SI_C_C2, col.header=TRUE, row.names = 1, info.cols= 1:2,
                         subset.column = 0, subset.rows = NA, data.type= "double")

RIS_C_C2_values <- WTcMC(RIS_C_C2, replicates = 999, weight = "equal", print.ris=TRUE)
sumMC.RInSp(RIS_C_C2_values)

#==================================
#IS measurements Grouped by cluster
#==================================
SI
model_C_pop <- lm (formula = d13C_PLAS ~ d13C_RBC, data = SI)
anova(model_C_pop)
summary(model_C_pop)
str((model_C_pop))

#cluster 1
SI_C1 <- subset(SI, cluster_2018.2019==1)
model_C_C1 <- lm (formula = d13C_PLAS ~ d13C_RBC, data = SI_C1)
anova(model_C_C1)
plot(model_C_C1)
summary(model_C_C1)
str(model_C_C1)

#cluster 2
SI_C2 <- subset(SI, cluster_2018.2019==2)
model_C_C2 <- lm (formula = d13C_PLAS ~ d13C_RBC, data = SI_C2)
anova(model_C_C2)
plot(model_C_C2)
summary(model_C_C2)
str(model_C_C2)

#graphs for linear regression of rbc vs plasma for each cluster
#Linear Regression  cluster 1
fitcarbon_c1 <- function(fitcluster1_c){ 
  m <- lm(d13C_PLAS ~ d13C_RBC, data=SI_C1);
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#plot of carbon - cluster 1 data
linearregplot_carbon_c1 <- ggplot(data= SI_C1, aes (x=d13C_RBC, y=d13C_PLAS)) +
  geom_point()+
  labs(title="",
       x="RBC"~delta^13~"C(???)",
       y="Plasma"~delta^13~"C(???)")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  theme(axis.title.x = element_text(margin = unit(c(1.5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 1.5, 0, 0), "mm")))
   #scale_x_continuous(limits=c(-24, -20))

#plot of carbon - cluster 2 data
linearregplot_carbon_c2 <- ggplot(data= SI_C2, aes (x=d13C_RBC, y=d13C_PLAS)) +
  geom_point()+
  labs(title="",
       x="RBC"~delta^13~"C(???)",
       y="")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  theme(axis.title.x = element_text(margin = unit(c(1.5, 0, 0, 0), "mm")))

#plot of carbon - full population data
linearregplot_carbon_pop <- ggplot(data= SI, aes (x=d13C_RBC, y=d13C_PLAS)) +
  geom_point()+
  labs(title="",
       x="RBC"~delta^13~"C(???)",
       y="")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  theme(axis.title.x = element_text(margin = unit(c(1.5, 0, 0, 0), "mm")))


linearregression_carbon_ggarrange_disrim_LE <- ggarrange(linearregplot_carbon_c1, linearregplot_carbon_c2, linearregplot_carbon_pop,
          nrow = 1, ncol = 3)

ggsave("linearregression_carbon_ggarrange_disrim_LE", linearregression_carbon_ggarrange_disrim, height = 5, width = 7, dpi = 600)

#BIC (mean sum of sq)
BIC_C_C1 = (2*1.008)/(2*(19-1)) #C1
BIC_C_C2 = (2*9.9183)/(2*(12-1)) #C2

#WIC mean sum of squares of the error
#sum of sq of error = sum of sq of residual error/deg free
WIC_C_C1 = 32.787/(19*(2-1)) #C1
WIC_C_C2 = 5.2942/(12*(2-1)) #C2

#F ratio (BIC:WIC) , model output
IS_pop_C_C1 = BIC_C_C1/WIC_C_C1 #C1
IS_pop_C_C1
IS_pop_C_C2 = BIC_C_C2/WIC_C_C2 #C2
IS_pop_C_C2

#Get WIC for each individual
#var(x)
mean(SI_C1$d13C_RBC, na.rm=TRUE) #c1
sd(SI_C1$d13C_RBC, na.rm=TRUE)
prob_C_c1_RBC <- dnorm(SI_C1$d13C_RBC, mean=-21.05381, sd=0.6096103) #RBC ind variance
E_C_c1_RBC = weighted.mean(x, p)
E_C_c1_RBC = weighted.mean(SI_C1$d13C_RBC, prob_C_c1_RBC)
E_C_c1_RBC = abs(E_C_c1_RBC)
ind_var_C_c1_RBC <- E*((x-mean(x))^2)
ind_var_C_c1_RBC <- E_C_c1_RBC*((SI_C1$d13C_RBC - mean(SI_C1$d13C_RBC))^2)

#IS for each individual
BIC/WIC(ind)
IS_C_c1_ind <-BIC_C_C1/ind_var_C_c1_RBC #c1, RBC

mean(SI_C2$d13C_RBC, na.rm=TRUE) #c2
sd(SI_C2$d13C_RBC, na.rm=TRUE)
prob_C_c2_RBC <- dnorm(SI_C2$d13C_RBC, mean=-20.7, sd=0.7586831) #RBC ind variance
E_C_c2_RBC = weighted.mean(x, p)
E_C_c2_RBC = weighted.mean(SI_C2$d13C_RBC, prob_C_c2_RBC)
E_C_c2_RBC = abs(E_C_c2_RBC)
ind_var_C_c2_RBC <- E*((x-mean(x))^2)
ind_var_C_c2_RBC <- E_C_c2_RBC*((SI_C2$d13C_RBC - mean(SI_C2$d13C_RBC))^2)

#IS for each individual
BIC/WIC(ind)
IS_C_c2_ind <-BIC_C_C2/ind_var_C_c2_RBC #c1, RBC


#=================================
#Compare IS values for both methods
#=================================
#make into dataframe
Index_name <- rep(c("BIC", "WIC", "IS"),6)
length(Index_name)
Index_group <- c(rep("Population",3),rep("Cluster 1",3), rep("Cluster 2",3),rep("Population",3),rep("Cluster 1",3), rep("Cluster 2",3))
Index_value <- c(0.3725333, 1.416581, 0.2629807,0.056,1.725632,0.03245189,0.9016636, 0.4421053,2.04373919,0.7525363, 0.7173667, 1.049025972, 0.613145, 0.9762869, 0.628037721, 0.7328559, 0.3074097, 2.383971293)
Index_method <- c(rep("GLM",9),rep("RInSp",9))
IS <- cbind.data.frame(Index_name, Index_group, Index_value, Index_method)
IS$Index_name <- factor(IS$Index_name ,                                   
                  levels = c("BIC", "WIC", "IS"))

Index_plot_C <- ggplot(IS, aes(x =Index_name, y = Index_value, fill = Index_group))+
  geom_bar(position='dodge', stat='identity')+
  facet_grid(. ~ Index_method)+
  xlab("") +
  ylab("Index value") +
  scale_fill_manual(values=c("#0072B2", "#E69F00", "#F0E442"))+
  theme_grey()+
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(face="bold"),
        strip.text.x = element_text(size=10, face= "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))

#IS_graph - poster
Element <- c("Carbon", "Carbon", "Carbon", "Nitrogen", "Nitrogen", "Nitrogen")
Group <- c("Cluster 1", "Cluster 2", "Population", "Cluster 1", "Cluster 2", "Population")
IS_value <- c(0.0324518785001669,2.04373918958401,0.262980584943607,0.15344,0.2952996,0.1716046)
IS_graph <- cbind.data.frame(Element, Group, IS_value)

IS_plot <- ggplot (IS_graph, aes(x=Element, y=IS_value, fill = Group))+
  geom_bar(position='dodge', stat='identity', width=.7)+
  ylab("Index") +
  scale_fill_manual(values=c("#276695", "#B47B08", "#FED57A"))+
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(size=14, face= "bold", vjust = -2),
          axis.title.x=element_blank(),
        axis.title.y = element_text(size=14, vjust = 4),
        axis.text.y = element_text(size=12),
        legend.text=element_text(size=11))+
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))+
  theme(plot.margin = margin(1,1,1,1, "cm"))
  



#individual IS scores
IS_ind_C_df <- cbind(Index$name, IS_C_ind, as.factor(SI_difference$cluster))

ind_cplot <- ggplot(IS_ind_C_df, aes(x =IS_C_ind, fill = as.factor(SI_difference$cluster)))+
  geom_histogram(position="dodge", bins = 5)+
  labs(title="Carbon",x="Specialization index value", y = "Count")+
    scale_fill_manual(labels = c("Cluster 1", "Cluster 2"), values=c("#0072B2", "#E69F00"))+
  labs(fill="")+
  scale_x_continuous(breaks = seq(0.00126062,0.24381498, by = .04), labels= c("0-0.039","0.04-0.079","0.08-0.119","0.12-0.159","0.16-0.199","0.2-0.239","0.24-0.279"))+
  theme_classic()+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)))+
  theme(axis.title.x = element_text(vjust = -0.75))

#Nitrogen
#GLMM
SI_difference
qqnorm(SI_difference$N_plasma, pch = 1, frame = FALSE)
qqline(SI_difference$N_plasma, col = "steelblue", lwd = 2)


model_N<- lm (formula = dN15_PLAS ~ dN15_RBC, data = SI)
summary(model_N)
anova(model_N)
plot(model_N)

#BIC (mean sum of sq)
BIC_N = (2*15.013)/(2*(31-1))

#WIC mean sum of squares of the error
#sum of sq of error = sum of sq of residual error/deg free
WIC_N = 14.423/(31*(2-1))

#IS
BIC_N/WIC_N 

#WIC foe each individual (Var)
mean(SI$dN15_PLAS, na.rm=TRUE)
sd(SI$dN15_PLAS, na.rm=TRUE)
prob_N_PLAS <- dnorm(SI$dN15_PLAS, mean=15.09094, sd=0.9952836) #PLASMA ind variance
E_N_PLAS = weighted.mean(x, p)
E_N_PLAS = weighted.mean(SI$dN15_PLAS, prob_N_PLAS, na.rm=TRUE)
E_N_PLAS = abs(E_C_PLAS)
ind_var_N_PLAS <- E*((x-mean(x))^2)
ind_var_N_PLAS <- E_N_PLAS*((SI$dN15_PLAS - mean(SI$dN15_PLAS, na.rm = TRUE))^2)

mean(SI$dN15_RBC, na.rm=TRUE)
sd(SI$dN15_RBC, na.rm=TRUE)
prob_N_RBC <- dnorm(SI$dN15_RBC, mean=15.18306, sd=0.7923954) #PLASMA ind variance
E_N_RBC = weighted.mean(x, p)
E_N_RBC = weighted.mean(SI$dN15_RBC, prob_N_RBC, na.rm=TRUE)
E_N_RBC = abs(E_C_RBC)
ind_var_N_RBC <- E*((x-mean(x))^2)
ind_var_N_RBC <- E_N_RBC*((SI$dN15_RBC - mean(SI$dN15_RBC, na.rm = TRUE))^2)


#add variaces from dependent variables
#Var(X+Y)=Var(X)+Var(Y)+2Cov(X,Y)
d15N_PLAS_na<-na.omit(SI$d15N_PLAS)
d15N_RBC_na<-na.omit(SI$d15N_RBC)
ind_var_N_RBC_na<-na.omit(ind_var_N_RBC)
ind_var_N_PLAS_na<-na.omit(ind_var_N_PLAS)

ind_var_N <- ind_var_N_PLAS + ind_var_N_RBC + 2*(cov(SI$dN15_PLAS, SI$dN15_RBC, use = "pairwise.complete.obs",method = "pearson"))
                                          
#IS for each individual
BIC/WIC(ind)
IS_N_ind <-BIC_N/ind_var_N
hist(IS_N_ind)  #check normality (response var in linear regression)
shapiro.test(IS_N_ind) # p<.05 - not normal

logtrans_indIS_N <- log10(IS_N_ind)
hist(logtrans_indIS_N)
shapiro.test(logtrans_indIS_N)

#IS from package
library("RInSp")
SI_N <- subset(SI_difference[,c(2,5,57,63)])
SI_N$N_plasma = abs(SI_N$N_plasma)
SI_N$N_RBC = abs(SI_N$N_RBC)
SI_N

SI_N_C1 <- subset(SI_N, cluster==1)
SI_N_C2 <- subset(SI_N, cluster==2)

RIS_N <- import.RInSp(SI_N, col.header=TRUE, row.names = 1, info.cols= 1:2,
                      subset.column = 0, subset.rows = NA, data.type= "double")

RIS_N_values <- WTcMC(RIS_N, replicates = 999, weight = "equal", print.ris=TRUE)
sumMC.RInSp(RIS_N_values)

#C1
RIS_N_C1 <- import.RInSp(SI_N_C1, col.header=TRUE, row.names = 1, info.cols= 1:2,
                         subset.column = 0, subset.rows = NA, data.type= "double")

RIS_N_C1_values <- WTcMC(RIS_N_C1, replicates = 999, weight = "equal", print.ris=TRUE)
sumMC.RInSp(RIS_N_C1_values)

#C2
RIS_N_C2 <- import.RInSp(SI_N_C2, col.header=TRUE, row.names = 1, info.cols= 1:2,
                         subset.column = 0, subset.rows = NA, data.type= "double")

RIS_N_C2_values <- WTcMC(RIS_N_C2, replicates = 999, weight = "equal", print.ris=TRUE)
sumMC.RInSp(RIS_N_C2_values)

#==================================
#IS measurements Grouped by cluster
#==================================
SI
model_N_pop <- lm (formula = dN15_PLAS ~ dN15_RBC , data = SI)
anova(model_N_pop)
summary(model_N_pop)
str(model_N_pop)

#cluster 1
SI_C1 <- subset(SI, cluster_2018.2019==1)
model_N_C1 <- lm (formula = dN15_PLAS ~ dN15_RBC, data = SI_C1)
anova(model_N_C1)
plot(model_N_C1)
summary(model_N_C1)
str(model_N_C1)

#cluster 2
model_N_C2 <- lm (formula = dN15_PLAS ~ dN15_RBC, data = SI_C2)
anova(model_N_C2)
plot(model_N_C2)
summary(model_N_C2)
str(model_N_C2)

#plot of nitrogen - cluster 1 data
linearregplot_nitrogen_c1 <- ggplot(data= SI_C1, aes (x=dN15_RBC, y= dN15_PLAS)) +
  geom_point()+
  labs(title="",
       x="RBC"~delta^15~"N(???)",
       y="Plasma"~delta^15~"N(???)")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  theme(axis.title.x = element_text(margin = unit(c(1.5, 0, 0, 0), "mm")),
    axis.title.y = element_text(margin = unit(c(0, 1.5, 0, 0), "mm")))

#plot of nitrogen - cluster 2 data
linearregplot_nitrogen_c2 <- ggplot(data= SI_C2, aes (x=dN15_RBC, y= dN15_PLAS)) +
  geom_point()+
  labs(title="",
       x="RBC"~delta^15~"N(???)",
       y="")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  theme(axis.title.x = element_text(margin = unit(c(1.5, 0, 0, 0), "mm")))

#plot of nitrogen - full population data
linearregplot_nitrogen_pop <- ggplot(data= SI, aes (x=dN15_RBC, y= dN15_PLAS)) +
  geom_point()+
  labs(title="",
       x="RBC"~delta^15~"N(???)",
       y="")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line=element_line(colour = "black"))+
  theme(axis.title.x = element_text(margin = unit(c(1.5, 0, 0, 0), "mm")))

ggarrange(linearregplot_nitrogen_c1, linearregplot_nitrogen_c2, linearregplot_nitrogen_pop, 
          ncol=3, nrow=1)

ggarrange(linearregplot_carbon_c1, linearregplot_carbon_c2, linearregplot_carbon_pop, 
          linearregplot_nitrogen_c1, linearregplot_nitrogen_c2, linearregplot_nitrogen_pop,
          ncol=3, nrow=2, labels = c("a", "b", "c", "d", "e", "f"))

#BIC (mean sum of sq)
BIC_N_C1 = (2*12.92)/(2*(19-1)) #C1
BIC_N_C2 = (2*1.0231 )/(2*(12-1)) #C2

#WIC mean sum of squares of the error
#sum of sq of error = sum of sq of residual error/deg free
WIC_N_C1 = 8.40/(19*(2-1)) #C1
WIC_N_C2 = 4.5953/(12*(2-1)) #C2

#F ratio (BIC:WIC) , model output
IS_pop_N_C1 = BIC_N_C1/WIC_N_C1 #C1
IS_pop_N_C1
IS_pop_N_C2 = BIC_N_C2/WIC_N_C2 #C2
IS_pop_N_C2

#Get WIC for each individual (see above from ALL clusters)
#add variances


#=================================
#Compare IS values for both methods
#=================================
#make into dataframe
Index_name_N <- rep(c("BIC", "WIC", "IS"),6)
length(Index_name_N)
Index_group_N<- c(rep("Population",3),rep("Cluster 1",3), rep("Cluster 2",3),rep("Population",3),rep("Cluster 1",3), rep("Cluster 2",3))
Index_value_N <- c(0.5004333,0.4652581,1.075603627,0.7177778,0.4421053,1.623544889,0.09288182,0.3826833,0.242711976,0.6768753,0.1211499,5.587089218,0.8684226,0.1250868,6.942559886,0.2695778,0.1149163,2.345862162)
Index_method_N <- c(rep("GLM",9),rep("RInSp",9))
IS_N <- cbind.data.frame(Index_name_N, Index_group_N, Index_value_N, Index_method_N)
IS_N$Index_name_N <- factor(IS_N$Index_name_N ,                                    # Change ordering manually
                        levels = c("BIC", "WIC","IS"))

Index_plot_N <- ggplot(IS_N, aes(x =Index_name_N, y = Index_value_N, fill = Index_group_N))+
  geom_bar(position='dodge', stat='identity')+
    xlab("") +
  ylab("Index value") +
  scale_fill_manual(values=c("#0072B2", "#E69F00", "#F0E442"))+
    facet_grid(. ~ Index_method_N)+
theme_grey()+
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(face="bold"),
        strip.text.x = element_text(size=10, face= "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))

ggarrange(Index_plot_C, Index_plot_N, common.legend=TRUE, legend = "right",
          nrow=2, ncol=1,labels = c("a", "b"))

#binding C and N into one dataframe for graphing GLMS only (no RInSP), Jacks suggestion
Index_name_C_N <- rep(c("BIC", "WIC", "IS"),6)
Index_group_C_N<- c(rep("Population",3),rep("Cluster 1",3), rep("Cluster 2",3),rep("Population",3),rep("Cluster 1",3), rep("Cluster 2",3))
Index_value_C_N <- c(0.3725333, 1.416581, 0.2629807,0.056,1.725632,0.03245189,0.9016636, 0.4421053,2.04373919, 0.5004333,0.4652581,1.075603627,0.7177778,0.4421053,1.623544889,0.09288182,0.3826833,0.242711976)
Index_method_C_N<- c(rep("Carbon", times = 9),
                       rep("Nitrogen", times = 9))

IS_C_N <- cbind.data.frame(Index_name_C_N, Index_group_C_N, Index_value_C_N, Index_method_C_N)
IS_C_N$Index_name_C_N <- factor(IS_C_N$Index_name_C_N ,                                    # Change ordering manually
                            levels = c("BIC", "WIC","IS"))


Index_plot_C_N <- ggplot(IS_C_N, aes(x =Index_name_C_N, y = Index_value_C_N, fill = Index_group_C_N))+
  geom_bar(position='dodge', stat='identity')+
  xlab("") +
  ylab("Index value") +
  scale_fill_manual(values=c("#0072B2", "#E69F00", "#F0E442"))+
  facet_grid(. ~ Index_method_C_N)+
  theme_grey()+
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(face="bold"),
        strip.text.x = element_text(size=10, face= "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))




#individual IS scores
IS_ind_N_df <- cbind.data.frame(SI_cleaned$name, IS_N_ind, as.factor(SI_cleaned$cluster))

range(IS_N_ind)

ind_nplot <- ggplot(IS_ind_N_df, aes(x =IS_N_ind, fill = as.factor(SI_cleaned$cluster)))+
  geom_histogram(position="dodge",bins = 10)+
  labs(title="Nitrogen",x="Specialization index value", y = "Count")+
  scale_fill_manual(labels = c("Cluster 1", "Cluster 2"), values=c("#0072B2", "#E69F00"))+
  labs(fill="")+
  scale_x_continuous(breaks = seq(0.003862924,0.312114710, by = .04))+
  theme_classic()+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)))+
  theme(axis.title.x = element_text(vjust = -0.75))


#plot RBC VS PLASMA
#Linear Regression CARBON , cluster 1
fit_C_PLASMA1 <- function(fitPLASMA){ 
  m <- lm(N_plasma ~ d13C_PLAS, data=subset(SI_difference, cluster==1));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#Linear Regression PLASMA, cluster 2
fitPLASMA2 <- function(fitPLASMA){ 
  m <- lm(N_plasma ~ C_plasma, data=subset(SI_difference, cluster==2));
  eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}

#plot of plasma data
ggplot(data= SI_difference, aes (x=C_plasma, y=N_plasma)) +
  geom_point(aes(color=factor(cluster)))+
  labs(title="PLASMA",
       x=delta^13~"C%",
       y= delta^15~"N%", 
       col="Cluster")+
  theme(plot.title = element_text(face="bold")) +
  geom_smooth(method = "lm", se = FALSE,aes(color = factor(cluster)))+
  theme_bw()+
  geom_text(x = -18.8, y = 15.5, label = fitPLASMA1(SI_difference), parse = TRUE, color="#F8766D")+   #red
  geom_text(x = -18.8, y = 15.4, label = fitPLASMA2(SI_difference), parse = TRUE, color="#00BFC4")    #blue


#plot to verify line & equation
dev.off() #turn off par function
plot(N_plasma ~ C_plasma, data=subset(SI_difference, cluster==1))
lm(N_plasma ~ C_plasma, data=subset(SI_difference, cluster==1))
abline(28.8197, 0.6244)

plot(N_plasma ~ C_plasma, data=subset(SI_difference, cluster==2))
lm(N_plasma ~ C_plasma, data=subset(SI_difference, cluster==2))
abline(24.1446, 0.3992)

#which factors influence individual IS values
#assumptions for anova:
#normal response var, equal variance
logtrans_indIS_C
logtrans_indIS_N
SI_difference_IS <-  cbind.data.frame(SI, forklength, ind_var_C_RBC, ind_var_C_PLAS, ind_var_N_RBC, ind_var_N_RBC, IS_C_ind, logtrans_indIS_C, IS_N_ind, logtrans_indIS_N)

shapiro.test(SI_difference_IS$FL)
shapiro.test(SI_difference_IS$IS_C_ind)
shapiro.test(SI_difference_IS$IS_N_ind)
shapiro.test(SI_difference_IS[SI_difference_IS$cluster == 1, "IS_C_ind"])
shapiro.test(SI_difference_IS[SI_difference_IS$cluster == 2, "IS_C_ind"])


#check for multicollinearity
library("corrplot")
multi.coll.test <- subset(SI_difference_IS[,c(3,5,73)]) 
corrplot(cor(multi.coll.test), method = "number") 

colnames(SI_difference_IS)[5] <- "cluster"
colnames(SI_difference_IS)[73] <- "FL"

#multiple linear regression
model_C_IS <- lm(logtrans_indIS_C ~ cluster + FL + Year, data = SI_difference_IS)
summary(model_C_IS)
AIC(model_C_IS)
10^0.34497  #backtransform estimates #cluster
10^-0.03404 #FL
10^-0.36154  #year 
10^0.20687  #backtransform std error
10^0.01671
10^0.23304

m1 <- glm(logtrans_indIS_C ~ cluster + FL + Year, data = SI_difference_IS, family = gaussian)
summary(m1)

m2 <- glm(logtrans_indIS_C ~ cluster +  Year, data = SI_difference_IS, family = gaussian)
summary(m2)

m3 <- glm(logtrans_indIS_C ~ cluster + FL, data = SI_difference_IS, family = gaussian)
summary(m3)

m4 <- glm(logtrans_indIS_C ~ FL + Year, data = SI_difference_IS, family = gaussian)
summary(m4)

AIC(m1)
AIC(m2)
AIC(m3)
AIC(m4)



# Calculate the AIC score for the initial model
initial_aic_C <- AIC(glm_model_C)

# Model optimization using AIC
optimized_model_C <- glm_model_C
optimized_aic_C <- initial_aic_C

# Iterate through different models by adding or removing predictors
predictors_C <- c("cluster", "fork_length", "year")

for (i in 1:length(predictors_C)) {
  candidate_model_C <- update(optimized_model_C, paste(". ~ . -", predictors_C[i]))
  candidate_aic_C <- AIC(candidate_model_C)
  
  if (candidate_aic_C < optimized_aic_C) {
    optimized_model_C <- candidate_model_C
    optimized_aic_C <- candidate_aic_C
  }
}

# Print the optimized model and AIC score
print(optimized_model_C)
print(optimized_aic_C)
 
res_C <- resid(glm_model_C)   #visually inspect residuals for normality
plot(fitted(glm_model_C), res_C)
abline(0,0)
qqnorm(res_C)
qqline(res_C) 
plot(density(res_C))

model_N_IS <- lm(logtrans_indIS_N ~ cluster + FL + Year, data = SI_difference_IS)
summary(model_N_IS)
10^0.40969   #backtransform estimates #cluster
10^-0.01936  #FL
10^-0.13152   #year 
10^0.18949  #backtransform std error
10^0.01531
10^0.21347

m1_nitrogen <- glm(logtrans_indIS_N ~ cluster + FL + Year, data = SI_difference_IS, family = gaussian)
summary(m1_nitrogen)
m2_nitrogen <- glm(logtrans_indIS_N ~ cluster + Year, data = SI_difference_IS, family = gaussian)
summary(m2_nitrogen)
m3_nitrogen <- glm(logtrans_indIS_N ~ cluster + FL, data = SI_difference_IS, family = gaussian)
summary(m3_nitrogen)
m4_nitrogen <- glm(logtrans_indIS_N ~ FL + Year, data = SI_difference_IS, family = gaussian)
summary(m4_nitrogen)
m5_nitrogen <- glm(logtrans_indIS_N ~ cluster, data = SI_difference_IS, family = gaussian)
summary(m5_nitrogen)

m5_nitrogen_lm <- lm(logtrans_indIS_N ~ cluster, data= SI_difference_IS)
summary(m5_nitrogen_lm)


10^0.3785    #backtransform estimates #cluster
10^0.1643  #cluster SE

AIC(m1_nitrogen)
AIC(m2_nitrogen)
AIC(m3_nitrogen)
AIC(m4_nitrogen)
AIC(m5_nitrogen)

glm_model_N <- glm(logtrans_indIS_N ~ cluster + FL + Year, data = SI_difference_IS, family = gaussian)
summary(glm_model_N)
# Calculate the AIC score for the initial model
initial_aic_N <- AIC(glm_model_N)

# Model optimization using AIC
optimized_model_N<- glm_model_N
optimized_aic_N <- initial_aic_N

# Iterate through different models by adding or removing predictors
predictors_N <- c("cluster", "FL", "Year")

for (i in 1:length(predictors_N)) {
  candidate_model_N <- update(optimized_model_N, paste(". ~ . -", predictors_N[i]))
  candidate_aic_N <- AIC(candidate_model_N)
  
  if (candidate_aic_N < optimized_aic_N) {
    optimized_model_N <- candidate_model_N
    optimized_aic_N <- candidate_aic_N
  }
}

# Print the optimized model and AIC score
print(optimized_model_N)
print(optimized_aic_N)

res_N <- resid(glm_model_N)   #visually inspect residuals for normality
plot(fitted(glm_model_N), res_N)
abline(0,0)
qqnorm(res_N)
qqline(res_N) 
plot(density(res_N))


SI_difference_IS
min(SI_difference_IS$IS_C_ind, na.rm = TRUE)
sd(SI_difference_IS$IS_C_ind, na.rm = TRUE)
min(SI_difference_IS$IS_N_ind, na.rm = TRUE)
sd(SI_difference_IS$IS_N_ind, na.rm = TRUE)

hist(SI_difference_IS$IS_C_ind)
hist(SI_difference_IS$IS_N_ind)
length(which(SI_difference_IS$IS_C_ind < 0.05))  #out of 45
length(which(SI_difference_IS$IS_N_ind < 0.05))

