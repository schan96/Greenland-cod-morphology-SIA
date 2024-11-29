#R version 4.1.1

#packages required:
#geomorph
#FactoMineR
#factoextra
#cluster
#mvnormtest
#lme4

**#Calculate 9 linear distances from landmarks (n=21)**

#load packages

library("geomorph")

#select TPS file

TPS <- readland.tps("Landmarks.tps", specID="ID")

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

landmarklinks<- define.links(meanshape) #create landmark links before proceeding

#plots landmark coordinates for a set of specimens

plotAllSpecimens(gdf_c, links=landmarklinks)

#perform PCA

dev.off()

PCA_geomorph <- gm.prcomp(fishgpa$coords)

plot(PCA_geomorph , main = "PCA")

PCA_geomorph

#Create data frame with trait start and end points

lmks <- data.frame(CP = c(6,8), BDP = c(5,9), BDA = c(3,10), PPFL = c(7,10), HD = c(11,21), HL=c(1,20),SL=c(1,16),EW=c(16,18),ML=c(14,15),
                   row.names = c("start", "end"))

#calculate linear distances

lineardists <- interlmkdist(TPS, lmks) #switched from gdf_c to TPS

#Import linear distances

dat <- read.csv("lineardistances_FL.csv")

#remove rows with NA (FL was not measured based on selection criteria)

dat2 <- na.omit(dat)

**#Allometry growth equation**

#CP

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

#BDP

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

#BDA

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

#PPFL

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

#HD

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

#HL

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

#SL

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

#EW

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

#ML

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

**#PCA on Linear Distances**

#loadpackages

library(FactoMineR)

library(factoextra)

library(cluster)

#Incorporate classifiers

#set factors

dat3$Year <- factor (dat2$Year, levels=c("2018", "2019"))

is.factor(dat3$Year)

head(dat3)

#PCA

dat4

res.pca <- PCA(dat4, graph = FALSE)

res.pca$ind$coord #individual coordinates

#visualize variances

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

#Extract the results for individuals

ind <- get_pca_ind(res.pca)

ind$coord

#Create dataframe with fish ID and individual coords from PCA

indcoord <- cbind.data.frame(dat2$X, ind$coord)

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

#Extract individuals to corresponding clusters

k$cluster<-as.factor(k$cluster)

ind_clusters <- cbind.data.frame(dat2$X, k$cluster, indcoord[,2:6])

colnames(ind_clusters)[1] <- "ID"

colnames(ind_clusters)[2] <- "cluster"

as.factor(ind_clusters$cluster)

ind_clusters

**#MANOVA on clusters**

library(mvnormtest)

manova <- manova(cbind(dat2$trans.CP.FINAL, dat2$trans.BDP.FINAL, dat2$trans.BDA.FINAL, dat2$trans.PPFL.FINAL, dat2$trans.HD.FINAL, dat2$trans.HL.FINAL, dat2$trans.SL.FINAL, dat2$trans.EW.FINAL, dat2$trans.ML.FINAL) ~ k$cluster, data=ind_clusters.lineardist)

**#T-test on linear distances among different clusters**

#example of t.test on each of the 9 linear distances (caudal peduncle) for cluster 1 vs 2

t.test(ind_clusters.lineardist$trans.CP.FINAL, as.numeric(ind_clusters.lineardist$cluster), p.adjust.method="bonferroni")

**#Stable isotope analysis**

SI <- read.csv("clusters_isotopes.csv")

#apply discrimination factors

SI$X13c.value1= SI$X13c.value1  - 1.2 #C PLAS

SI$X13c.value2= SI$X13c.value2 - 0.1 #C RBC

SI$X15n.value1 = SI$X15n.value1 - 0.9 #N PLAS

SI$X15n.value2 = SI$X15n.value2 - 1.1  #N RBC

#import formatted SI data

SI_difference_raw <- read.csv("clusters_isotopes_discriminated.csv")

#calculate habitat and trophic switch values (plasma-RBC)

C_difference <- as.data.frame.vector(SI_difference_raw$C_RBC - SI_difference_raw$C_plasma)

N_difference <- as.data.frame.vector(SI_difference_raw$N_RBC - SI_difference_raw$N_plasma)

SI_difference <- cbind(SI_difference_raw, C_difference, N_difference)

colnames(SI_difference)[5] <- "cluster"

colnames(SI_difference)[65] <- "C_difference"

colnames(SI_difference)[66] <- "N_difference"

#test for difference between clusters for C and N (plasma vs rbc)

t.test(SI_difference[SI_difference$cluster == 1, "C_plasma"],SI_difference[SI_difference$cluster == 2, "C_plasma"])

t.test(SI_difference[SI_difference$cluster == 1, "N_plasma"],SI_difference[SI_difference$cluster == 2, "N_rbc"])

#test for difference between clusters for C_difference and N_difference (plasma vs rbc)

t.test(SI_difference[SI_difference$cluster == 1, "C_difference"],SI_difference[SI_difference$cluster == 2, "C_difference"])

t.test(SI_difference[SI_difference$cluster == 1, "N_difference"],SI_difference[SI_difference$cluster == 2, "N_difference"])

**#calculation of BIC, WIC  for clusters**

library(lme4)

SIA_C_clus<-read.csv("SIA_C_clus.csv")

SIA_N_clus<-read.csv("SIA_N_clus.csv")

#these models include fixed and random factors to try and remove the variance explained by 'length'

SIA_N_clus_1=SIA_N_clus[SIA_N_clus$cluster_2018.2019==1,]

SIA_N_clus_2=SIA_N_clus[!SIA_N_clus$cluster_2018.2019==1,]

WIC_N1<-lmer(X15n_value_lediscrim ~ FL_mm + X15n_tissue_lediscrim + (1|ID), data = SIA_N_clus_1)

WIC_N2<-lmer(X15n_value_lediscrim ~ FL_mm + X15n_tissue_lediscrim + (1|ID), data = SIA_N_clus_2)

SIA_C_clus_1=SIA_C_clus[SIA_C_clus$cluster_2018.2019==1,]

SIA_C_clus_2=SIA_C_clus[!SIA_C_clus$cluster_2018.2019==1,]

WIC_C1<-lmer(X13c_value_lediscrim ~ FL_mm + X13c_tissue_lediscrim + (1|ID), data = SIA_C_clus_1)

WIC_C2<-lmer(X13c_value_lediscrim ~ FL_mm + X13c_tissue_lediscrim + (1|ID), data = SIA_C_clus_2)

summary(WIC_N1)

#performing a summary() on the model fits provides you with the relevant information under Random effects (intercept and residual)

**#calculation of individual WIC**

SIA_C_pair<-read.csv("SIA_C_pair.csv")

SIA_N_pair<-read.csv("SIA_N_pair.csv")

#Nitrogen

mod_n<-(lm(dN_p ~ dN_rbc, data = SIA_N_pair))  # Fit the model

SIA_N_pair$predpred <- predict(mod_n)   # Save the predicted values

SIA_N_pair$WIC <- abs(residuals(mod_n)) #difference bw predicted and the actual data, in absolute value (degree to which they differ)

substrRight <- function(SIA_N_pair, n) {
  substr(SIA_N_pair, nchar(SIA_N_pair) - n + 1, nchar(SIA_N_pair))
}

SIA_N_pair$Yr<-substrRight(SIA_N_pair$ID,4)

hist(SIA_N_pair$WIC)

log10_WIC_N<- log10(SIA_N_pair$WIC)

hist(log10_WIC_N)

shapiro.test(log10_WIC_N)

group <- as.factor(SIA_N_pair$cluster_2018.2019)  #equal variances

group2 <- as.factor(SIA_N_pair$Yr)

leveneTest(FL_mm ~ group, data = SIA_N_pair)

leveneTest(FL_mm ~ group2, data = SIA_N_pair)

SIA_N_pair_log <- cbind.data.frame(SIA_N_pair, log10_WIC_N)

summary(glm(log10_WIC_N~cluster_2018.2019+Yr+FL_mm, data=SIA_N_pair_log)) #based on AIC

10^0.144744    #backtransform estimates #cluster

#Carbon

mod_c<-(lm(dC_p ~ dC_rbc, data = SIA_C_pair))  # Fit the model

SIA_C_pair$predpred <- predict(mod_c)   # Save the predicted values

SIA_C_pair$WIC <- abs(residuals(mod_c)) #difference bw predicted and the actual data, in absolute value (degree to which they differ)

substrRight <- function(SIA_C_pair, n) {
  substr(SIA_C_pair, nchar(SIA_C_pair) - n + 1, nchar(SIA_C_pair))
}

SIA_C_pair$Yr<-substrRight(SIA_C_pair$ID,4)

hist(SIA_C_pair$WIC) #no need to transform. looks normal

shapiro.test(log10_WIC_C)

hist(SIA_C_pair$WIC)

log10_WIC_C<- log10(SIA_C_pair$WIC)

hist(log10_WIC_C)

shapiro.test(log10_WIC_C)

group <- as.factor(SIA_C_pair$cluster_2018.2019)  #equal variances

group2 <- as.factor(SIA_C_pair$Yr)

leveneTest(FL_mm ~ group, data = SIA_C_pair)

leveneTest(FL_mm ~ group2, data = SIA_C_pair)

SIA_C_pair_log <- cbind.data.frame(SIA_C_pair, log10_WIC_C)

summary(glm(log10_WIC_C~cluster_2018.2019+Yr+FL_mm, data=SIA_C_pair_log)) #based on AIC

10^-0.3092933   #backtransform estimates #cluster

summary(glm(WIC~cluster_2018.2019+Yr+FL_mm, data=SIA_C_pair))
