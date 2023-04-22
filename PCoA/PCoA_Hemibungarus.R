#Code for generating Principal Coordinate Analysis for Hemibungarus

#Set working directory
getwd()
setwd("C:/Users/Justin/Documents/Publications/Hemibungarus_NewRecord_2019/Quant_Stats/Hemibungarus_PCoA/")

#Load packages
library(vegan)
library(labdsv)
library(ggplot2)
library(dplyr)


Hemibungarus_Adults <- read.csv("Hemibungarus_Numericc.csv", header = TRUE, na.strings = c("","NA"))

testpcoa <- Hemibungarus_Adults[c(6:8,11:16,18:31)]

#calculate distance matrix with Gower transform in Vegan

Gowerdist <- vegdist(testpcoa,method="gower")

#This spits out an error about missing values; this can be fixed with:
Gowerdist <- vegdist(testpcoa,method="gower", na.rm= TRUE)

#error, must be numeric
sapply(Hemibungarus_Adults, class)  #lets you know which columns are factors (i.e., a no-no as factors are not allowed)

Gowerdist

# Documentation is a little sparse, From doc re: na.rm "logical. Should missing values (including NaN) be omitted from the calculations?"
# Need to check and see if this is solved by ignoring columns. Can do by making dataset omitting columns with missing data to see if identical results.
#Note: na.rm does not remove columns, only cells!

#Run PCoA in labdsv saving first four dimensions
pcotest <- pco(Gowerdist,k=4)

#Save output as .csv

write.csv(pcotest$points,'pcotestPOINTS.csv')

#Visualize plot

#Don't forget to append a species column to the outputted .csv (also, can turn the species into numbers to make a gradient color scheme)

pcoaTESTplot <- read.csv("pcotestPOINTS.csv")       
pcoaTESTplot.4sp <- read.csv("pcotestPOINTS_4sp.csv") #3 species; can also subsitute with the file 'pcotestPOINTS_4sp.csv'
ggplot(pcoaTESTplot.4sp, aes(x=V1, y=V2, color=Species)) + geom_point()

ggplot(pcoaTESTplot.4sp, aes(x=V1, y=V2, color=Species)) + geom_point() + geom_text(aes(label=Species),hjust=0, vjust=0) +
  scale_color_manual(values=c("#7ECEFA", "#93E07F", "#B971F0", "#E0A618")) 


#The next chunk of code is for  experimenting with drawing convex hulls
#calculating a convex polygon hull for each data group
hull_fr<-pcoaTESTplot.4sp%>%
  group_by(Species)%>%
  slice(chull(V1,V2))
#now create base scatterplot
fr.pco.plot<-ggplot(pcoaTESTplot.4sp,aes(x=V1,y=V2,colour=Species))+
  geom_point()  +
  scale_color_manual(values=c("#7ECEFA", "#93E07F", "#B971F0", "#E0A618"))
#put in a fill group on the plot and overlay the hulls
fr.pco.plot+aes(fill=factor(Species))+geom_polygon(data=hull_fr,alpha=0.3)
#This block works for convex hulls
fr.pco.plot+aes(fill=factor(Species))+geom_polygon(data=hull_fr,alpha=0.3) +
  scale_fill_manual(values=c("#7ECEFA", "#93E07F", "#B971F0", "#E0A618"))


# Generating a 3D, rotating video of your PCA (note: must have more than 1 point and multiple groups)
library(scatterplot3d)
library(rgl)
library(ggplot2)
library(ggfortify)
#library(ggbiplot)
library(readr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(heplots)
library(candisc)
library(rpart)
library(rpart.plot)
library(corrplot)
library(magrittr)
library(pca3d)
library(scatterplot3d)
library(plot3D)
library(cluster)

# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp1_3<-c("#E69F00", "#56B4E9", "#009E73")
cbp1_4<-c("#E69F00", "#56B4E9", "#009E73", "#999999")
cbp1_5<-c("#E69F00", "#56B4E9", "#009E73", "#999999", "#CC79A7")
cbp1_4b<-c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp1_6<-c("#E69F00", "#56B4E9", "#009E73", "#999999", "#D55E00", "#CC79A7")

scatter3d(V1~V2+V3|Species,data=pcoaTESTplot.3sp,
          surface=TRUE, sphere.size=0.8,ellipsoid=TRUE, level=0.5, 
          surface.col=c("#E69F00", "#56B4E9", "#009E73", "#999999", "#D55E00", "#CC79A7"))

#No ellipsoids
scatter3d(x = pcoaTESTplot.3sp$V1, y = pcoaTESTplot$V2, z = pcoaTESTplot$V3,
          point.col = "blue", groups = as.factor(pcoaTESTplot.3sp$Species), grid = FALSE, surface = FALSE)


#With ellipsoids : this isn't working, but may be due to how many points per grouping we have
scatter3d(x = pcoaTESTplot.3sp$V1, y = pcoaTESTplot$V2, z = pcoaTESTplot$V3,
          point.col = "blue", groups = as.factpr(pcoaTESTplot.3sp$Species), ellipsoid = TRUE, grid = TRUE, surface = TRUE)


#No ellipsoids 3 species version
scatter3d(x = pcoaTESTplot.3sp$V1, y = pcoaTESTplot.3sp$V2, z = pcoaTESTplot.3sp$V3,
          point.col = "blue", groups = as.factor(pcoaTESTplot.3sp$Species), grid = FALSE, surface = FALSE)

pcoaTESTplot.3sp <- read.csv("pcotestPOINTS_3sp.csv")
#With ellipsoids 3 species version
scatter3d(x = pcoaTESTplot.3sp$V1, y = pcoaTESTplot.3sp$V2, z = pcoaTESTplot.3sp$V3,
          point.col = "blue", groups = as.factor(pcoaTESTplot.3sp$Species), ellipsoid = TRUE, grid = TRUE, surface = FALSE)

#With ellipsoids 4 species version: This version works because before cf. geminaulis and cf. calligaster had too few points to plot as an ellipse

pcoaTESTplot.4sp <- read.csv("pcotestPOINTS_4sp.csv")
scatter3d(x = pcoaTESTplot.4sp$V1, y = pcoaTESTplot.4sp$V2, z = pcoaTESTplot.4sp$V3,
          point.col = "blue", groups = as.factor(pcoaTESTplot.4sp$Species), ellipsoid = TRUE, grid = TRUE, surface = FALSE)

#Rotation Gif
library(tidyverse)
library(dplyr)
library(rgl)
library(ggpubr)
library(heplots)
library(candisc)
library(rpart)
library(rpart.plot)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(magrittr)
library(pca3d)
library(scatterplot3d)
library(plot3D)
library(cluster)
library(sf)
library(sp)
library(GISTools)

scatter3d(V2~V1+V3|as.factor(Species),data=pcoaTESTplot.4sp, fogtype ="none", fov=20,  box=TRUE, xlab="V1",ylab="V2",zlab="V3",
          bg.col= "black",
          axis.col=c("white", "white", "white"), axis.scales=FALSE, axis.ticks=TRUE, text.col="pink",
          surface=FALSE, fill=TRUE, grid=TRUE, sphere.size=0.8, ellipsoid=TRUE, level=0.5,
          surface.col=c("#7ECEFA", "#93E07F", "#B971F0", "#E0A618"), font="10")


grid3d(c("x", "y", "z"), at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)

#You must leave the temporary window from line 165 (scatter3d) open and the images WILL be the size of the window; maximize the window for best results and THEN run the export code below
movie3d(spin3d(axis = c(0, 1, 0)), duration = 12,
        dir = "C:/Users/Justin/Documents/Publications/Hemibungarus_NewRecord_2019/Quant_Stats/Hemibungarus_PCoA/", convert=TRUE)


##duration = 12 gave me a full rotation for my data


