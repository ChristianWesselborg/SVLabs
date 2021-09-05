#Author Christian Wesselborg
#used for extracting which clusters of a GROMACS simululation should be used for comparison
#only clusters with >10% of the total models should be used

#clear memeory
rm(list = ls())

#load data
list <- read.csv("c:/users/wessec/documents/research/rpi gb/pdbfiles/dimersim/topbot/clust-size.xvg", sep = "", skip = 18, header = F)

colnames(list) <- c("cluster", "size")

#calculate

cutoff <- sum(list$size) #* 0.1

plot(list$size ~ list$cluster)
