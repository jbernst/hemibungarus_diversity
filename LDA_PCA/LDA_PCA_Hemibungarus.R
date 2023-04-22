### --- Script to analyze phenotypic data for Hemibungarus species --- ###
### --- Author of original code: J.A. Soto-Centeno (2019) --- ###
### --- Modified by Justin Bernstein: 24 April 2020 --- ###

## Main script sections:
## 1. Packages and dependencies
## 2. Preliminary data management & imputation (mice)
## 3. Data wrangling to create multiple partitions for LDA & PCA

setwd("C:/Users/Justin/Documents/Coding_Practice_R/LDA/HemiLDA_CF/")

# 1. load packages
# data imputation - to correct data matrices if you have missing values: statistical amputation
library(mice)
library(ade4)
# data wrangling and plotting
library(ggbiplot)
library(tidyverse)
library(gridExtra)
# machine learning based linear discriminant analyses
library(MASS)
library(car)
library(caret)

# 1. Import the data
# read the Hemibungarus data
data  <- read.csv(file = "Hemibungarus_AN_CF.csv", header = TRUE) 

# examine the dimensions of the imported data.frame
dim(data)

# examine the header of the imported data.frame
head(data)

# summarize the imported data.frame
summary(data)
View(data)

# 2. Data imputation using MICE: works with continuous data really well
# NOTE: to reduce bias in the resulting imputation, select a priori partitions of the dataset
# (e.g. by country) to ensure consistent measurements & a low missing data ratio. Data imputation
# is unreliable for partitions with >60% missing data (see Penone et al. 2014 - Met. Ecol. Evol.).
# MICE will automatically discard (i.e. not fill in) those partitions with a large proportion of 
# missing data.

# Perform multiple data imputation for each partition in the data.frame. The range of continuous 
# data colums and rows must be specified for each partition ******************NOT NEEDED IF NO MISSING DATA
part1 <- mice(data[1:9,7:24])   #These are partitioned by geographically proximate locality (for phenotypic similarity) 
part2 <- mice(data[10:25,7:24]) # ... etc. for each partition with missing cells

# obtain the imputed data matrix for each partition
complete1 <- complete(part1) #The complete function will give the complete matrix (the previous is computing the imputation)
complete2 <- complete(part2) #... etc. for each partition

### Used custom dplyr scripts to compile all of these matrices into one *.csv data.frame ###

# 3. Morphological analyses (LDA & PCA)

# use dplyr::filter to make different data partitions (comment and uncomment as needed)
# A. male vs female partition to test for sexual dimorphism within the group
#mf <- data %>% filter(Sex == "F" | Sex == "M") #LDA works with 3 groups or more, so this may get fuzzy

#data wrangle to create groups based on species, or even geography (e.g., island vs. mainland), and other sets
# B. morphological comparison of the current species of Hemibungarus: mcclungi vs gemianulis vs calligaster vs cf. mcclungi
csv <- data %>% filter(species == "calligaster" | species == "mcclungi" | species == "gemianulis" | species =="cf_mcclungi")        

## 3B. full comparison of recognized forms (4 groups here)
# log transform the morphological variables to normalize them
csv_ln <- log(csv[, c(6:25)])  #log transform (only a certain number of variables): this is all rows, and columns 6-25; we are normalizing here 

# run PCA, examine PC contributions & rotation
csv_pca <- prcomp(csv_ln, center = T, scale. = T)
summary(csv_pca)
csv_pca$rotation

#create scree plot of pca proportion of variance
ggscreeplot(csv_pca)

# create pca biplot
# set the theme of the plot
theme_set(theme_bw())
# set raw plot
csv_fig <- ggbiplot(csv_pca, obs.scale = 1, var.scale = 1, groups = csv$species, ellipse = T,  #obs.scale and var.scale will help plot scale to the axes so it looks nicer. not scaling the data. it is grouped by the original dataframe's groupings (Ssp here); ellipse = T gives data ellipses. Circle = F cuz the circles don't really mean anything. var.axes is the arrows (it does PC1 by default but can specify more); size = size of dota; alphe = transparency of dots (good for ones that are super close together)
                    circle = F, var.axes = F, size = 2.5, alpha = 0.75)          #ggbiplot PCA plot by ggplot; axes here will determine axis labels in line 129+
csv_fig
    ## change var.axes = T to see the character vectors
    ## ellipses represent 68% Gausian data ellipses ~ standard deviation & size shows the variance
# set final plot
csv_fig + scale_color_manual(values = c("#7ECEFA", "#E0A618", "#B971F0", "#93E07F"),
                             labels = c("calligaster", "mcclungi", "gemianulis","cf_mcclungi")) +
  labs(x = "PC1 (30.5%)", y = "PC2 (16.7%)") +  #these are the axes (outside the color code above)   #remember to change the axes names for different analyses! 
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),  #sizes of axes 
        legend.title = element_blank(), legend.text = element_text(face = "italic", size = 15))  #make fonts italics

# OPTIONAL: make a marginal density plot
# source: http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
# 3xtract csv_pca$x <- as.data.frame(csv_pca)

# run LDA
# create a data.frame with island names and the log transformed characters 
# extract species variable from the calligaster, mcclungi, cf mcclungi, gemianulis data.frame  
csv_names <- csv[4]                                                           

# combine csv_names with csv_ln
csv_nam_ln <- cbind(csv_names, csv_ln)

### write to new d.f  - this prevents errors when doing the confusion matrix 
write.csv(csv_nam_ln, file = "csv_nam_ln.csv", row.names = F)

# reopen the file
csv_nam_ln <- read.csv(file = "csv_nam_ln.csv", header = T)  #The LDA on line 156 said variable 18 was constant; so I removed it in the csv_name_ln.csv. If the line 144 code is rerun, it will add Variable 18 back in

# make scatterplots of measurements
# LDA assumes variables are normally distributed. scatterplots help you see that
scatterplotMatrix(csv_nam_ln[2:9]) # do it in smaller groups so it is readable 
scatterplotMatrix(csv_nam_ln[10:15])                                          
scatterplotMatrix(csv_nam_ln[16:21])                                          

#Some of the variables are not normally distributed; these should be removed later
#Note: This was done in the publication, but this current script still includes these variables.

# run linear discriminant analysis
# ~ . in the argument means that we'll use all other variables than Ssp or Country as covariates
csv.lda <- lda(species ~ ., data = csv_nam_ln) # ignore the error that "group new is empty" "Ssp~. says take Ssp from the current dataframe; fix error by above code

csv_nam_ln2 <- csv_nam_ln[, -c(18)] #This will get rid of column 18

csv.lda #this is what you are going to pull your values out of
## display the results of lda. note the proportion of trace
## the values presented by the lda function is the % separation achieved by each discriminant function

# confusion matrix of LDA machine learning classification
# assess the prediciton accuracy of the LDA from the model to the actual data [using package caret]
csv.lda.predict <- train(species ~ ., method = "lda", data = csv_nam_ln)        

# calculate a confusion matrix to see the accuracy of classification. Conf. matrix tells you have well your classification model works. It has your prediction and the references. Not bad but 1 calligaster got classified as a cf_mcclungi
confusionMatrix(as.factor(csv_nam_ln$species), predict(csv.lda.predict, csv_nam_ln))       

# to interpret Kappa, see here: https://stats.stackexchange.com/questions/82162/cohens-kappa-in-plain-english


# create stacked histograms of the LDA values for the groups
csv.lda.values <- predict(csv.lda)
ldahist(csv.lda.values$x[ ,1], g = csv_nam_ln$species, #x is going to be the individual LDA values we are looking for
        type = "histogram", col = "light gray") # [ ,1] indicates LD1, change number for the LD you want to plot; can change type to density for density plots
## note: if "margins are too large", expand the plot viewer & do this:
par(mar = c(4,3,1,1)) # tweak each side of the plot to make it viewable

# create a scatterplot of LD to visualize discrimination
# convert LD data into a data.frame for cinereus/semotus/villosissimus
csv_LD <- data.frame(type = csv_nam_ln[,1], lda = csv.lda.values$x)
### step-by-step ggplot for LDA with all features
# set the plot theme
theme_set(theme_bw())
# basic plot for LDA
ggplot(csv_LD, aes(lda.LD1, lda.LD2, colour = type)) +
  geom_point()
# add rug
ggplot(csv_LD, aes(lda.LD1, lda.LD2, colour = type)) + 
  geom_point() + 
  geom_rug()
# change point size
ggplot(csv_LD, aes(lda.LD1, lda.LD2, colour = type)) + 
  geom_point(size = 2.5, alpha = 0.65) + 
  geom_rug()
# add 68% ellipse to points. this ellipse are proportional to standard dev & size indicates magnitude of variance
csvFig1 <- ggplot(csv_LD, aes(lda.LD1, lda.LD2, colour = type)) + 
  geom_point(size = 2.5, alpha = 0.65) + stat_ellipse(level = 0.68) +
  geom_rug()
# manually change the colors
csvFig1 + scale_color_manual(values = c("#7ECEFA", "#E0A618", "#B971F0", "#93E07F"))   # color order is calligaster gemianulis mcclungi 


# manually change the legend labels, axis labels, then set their sizes
csvFig1a <- csvFig1 + scale_color_manual(values = c("#7ECEFA", "#E0A618", "#B971F0", "#93E07F"),                
                            labels = c("H. calligaster", "H. mcclungi", "H. gemianulis", "H. cf. mcclungi")) +  
  labs(x = "LD1", y = "LD2") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.position = "none")
        legend.title = element_blank()

# OPTIONAL: make a marginal density plot
# source: http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
csvFig1density <- ggplot(csv_LD, aes(lda.LD1, fill = type)) + 
  geom_density(alpha = 0.55) + 
  scale_fill_manual(values = c("#7ECEFA", "#E0A618", "#B971F0", "#93E07F"),
                     labels = c("H. calligaster", "H. mcclungi", "H. gemianulis", "H. cf. mcclungi")) + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 18))
    ## removed the legend and x axis title to make cleaner when plotting against the scatterplot

# arrange both LD scatterplot & density plot of LD1 using gridExtra::grid.arrange
grid.arrange(csvFig1a, csvFig1density,
             ncol = 1, nrow = 2, heights = c(4, 1.4))

##################Chagne the Density plot (Will cause overwriting of the above Figures)###################
view(csv.lda.values$x)

csv_LD2 <- data.frame(type = csv_nam_ln[,1], lda = csv.lda.values$x[,2])
View(csv_LD2)

csv_LD <- data.frame(type = csv_nam_ln[,1], lda = csv.lda.values$x)
### step-by-step ggplot for LDA with all features
# set the plot theme
theme_set(theme_bw())
# basic plot for LDA
ggplot(csv_LD, aes(lda.LD2, lda.LD1, colour = type)) +
  geom_point()
# add rug
ggplot(csv_LD, aes(lda.LD2, lda.LD1, colour = type)) + 
  geom_point() + 
  geom_rug()
# change point size
ggplot(csv_LD, aes(lda.LD2, lda.LD1, colour = type)) + 
  geom_point(size = 2.5, alpha = 0.65) + 
  geom_rug()
# add 68% ellipse to points. this ellipse are proportional to standard dev & size indicates magnitude of variance
csvFig1 <- ggplot(csv_LD, aes(lda.LD2, lda.LD1, colour = type)) + 
  geom_point(size = 2.5, alpha = 0.65) + stat_ellipse(level = 0.68) +
  geom_rug()

# manually change the colors
csvFig1 + scale_color_manual(values = c("#7ECEFA", "#E0A618", "#B971F0", "#93E07F"))   #color order is calligaster gemianulis mcclungi


csvFig1a <- csvFig1 + scale_color_manual(values = c("#7ECEFA", "#E0A618", "#B971F0", "#93E07F"),                
                                         labels = c("H. calligaster", "H. mcclungi", "H. gemianulis", "H. cf. mcclungi")) +          
  labs(x = "LD2", y = "LD1") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.position = "none")
legend.title = element_blank()

csvFig2density <- ggplot(csv_LD2, aes(lda, fill = type)) + 
  geom_density(alpha = 0.55) + 
  scale_fill_manual(values = c("#7ECEFA", "#E0A618", "#B971F0", "#93E07F"),
                    labels = c("H. calligaster", "H. mcclungi", "H. gemianulis", "H. cf. mcclungi")) + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 18))


grid.arrange(csvFig1a, csvFig2density,
             ncol = 1, nrow = 2, heights = c(4, 1.4))


