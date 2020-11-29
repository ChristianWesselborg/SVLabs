#install.packages("ggplot2")
#install.packages("stringr")
#install.packages("gtools")
library(ggplot2)
library(stringr)
library(gtools)

home <- "g:/research/belfortlocalfiles/data/zdock_docking/"
masterlist <- dir(home)
sublist <- grep("output", masterlist)
sublist <- masterlist[sublist]

reslist <- NULL
liglist <- NULL
for (i in sublist){
  temp <- read.csv(paste(home, sublist[3], sep = ""), skip = 1, header = F, sep = " ", stringsAsFactors = F)
  
  for (j in 1:20){
    resextract <- strsplit(temp[j, 5], split = ",")
    reslist <- c(reslist, resextract[[1]])
    ligextract <- strsplit(temp[j, 9], split = ",")
    liglist <- c(liglist, ligextract[[1]])
  }
}

test <- c("B23", "D17", "C14")
sort <- str_sort(test, numeric = TRUE)


hold <- str_sort(liglist, numeric = T)
hold <- liglist[order(nchar(liglist), liglist)]
hold <- mixedsort(liglist)

#Wild Type 1
frameWT <- data.frame(liglist, NA)
colnames(frameWT) <- c("liglist", "color")
frameWT[frameWT$liglist == "ALA2U", 2] <- "red"

ggplot() +
  geom_bar(aes(x = frameWT$liglist, fill = frameWT$color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  ggtitle("WT")

#A2T
frameT <- data.frame(liglist, NA)
colnames(frameT) <- c("liglist", "color")
frameT[frameT$liglist == "VAL2U", 2] <- "red"

ggplot() +
  geom_bar(aes(x = frameT$liglist, fill = frameT$color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  ggtitle("A2T")

#A2V 2
frameV <- data.frame(liglist, NA)
colnames(frameV) <- c("liglist", "color")
frameV[frameV$liglist == "VAL2U", 2] <- "red"

ggplot() +
  geom_bar(aes(x = frameV$liglist, fill = frameV$color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  ggtitle("A2V")

##### 07-28-2020 #####
#Organize zdock runs for 1-200 models
library(plyr)
library(ggplot2)

rm(list = ls())

##dim_lateral
#A2T
rm(list = ls())

#path <- "g:/research/belfortlocalfiles/zdock output/transfer/dim_lateralA2T/results/20200726_153713_output.txt"
path <- "c:/users/wessec/documents/research/rpi gb/zdock output/transfer/dim_lateralA2T/results/20200726_153713_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_lateralA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_lateralreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAT) <- c("4pe5_reslist", "receptor_freq")
dim_lateralligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAT <- rbind.fill(dim_lateralreceptorAT, pereslist)
dim_lateralligandAT <- rbind.fill(dim_lateralligandAT, naoreslist)

dim_lateralligandAT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAT[dim_lateralreceptorAT$`4pe5_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAT[dim_lateralreceptorAT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAT$`2nao_reslist` <- factor(dim_lateralligandAT$`2nao_reslist`, levels = dim_lateralligandAT$`2nao_reslist`)

dim_lateralligandAT$rel_freq <- dim_lateralligandAT$ligand_freq / max(dim_lateralligandAT$ligand_freq)
dim_lateralreceptorAT$rel_freq <- dim_lateralreceptorAT$receptor_freq / max(dim_lateralreceptorAT$receptor_freq)

#plot A2T ligand
dim_lateralligandAT$color <- dim_lateralligandAT$`2nao_reslist` == "THR2U" | dim_lateralligandAT$`2nao_reslist` == "THR2X"

ggplot(data = dim_lateralligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Lateral Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralreceptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Lateral Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Lateral Dimer Ligand Interaction")

ggplot(data = dim_lateralreceptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Lateral Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_lateralligandAT, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralligandAT.csv")

write.csv(dim_lateralreceptorAT, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralreceptorAT.csv")

#A2V
rm(list = ls())

path <- "g:/research/belfortlocalfiles/zdock output/transfer/dim_lateralA2V/results/20200726_155055_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_lateralA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_lateralreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAV) <- c("4pe5_reslist", "receptor_freq")
dim_lateralligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAV) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAV <- rbind.fill(dim_lateralreceptorAV, pereslist)
dim_lateralligandAV <- rbind.fill(dim_lateralligandAV, naoreslist)

dim_lateralligandAV$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAV$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAV[dim_lateralreceptorAV$`4pe5_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAV[dim_lateralreceptorAV$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAV$`2nao_reslist` <- factor(dim_lateralligandAV$`2nao_reslist`, levels = dim_lateralligandAV$`2nao_reslist`)

dim_lateralligandAV$rel_freq <- dim_lateralligandAV$ligand_freq / max(dim_lateralligandAV$ligand_freq)
dim_lateralreceptorAV$rel_freq <- dim_lateralreceptorAV$receptor_freq / max(dim_lateralreceptorAV$receptor_freq)

#plot A2V ligand
dim_lateralligandAV$color <- dim_lateralligandAV$`2nao_reslist` == "VAL2U" | dim_lateralligandAV$`2nao_reslist` == "VAL2X"

ggplot(data = dim_lateralligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Lateral Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralreceptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Lateral Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Lateral Dimer Ligand Interaction")

ggplot(data = dim_lateralreceptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Lateral Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_lateralligandAV, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralligandAV.csv")

write.csv(dim_lateralreceptorAV, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralreceptorAV.csv")

#WT
rm(list = ls())

path <- "g:/research/belfortlocalfiles/zdock output/transfer/dim_lateralWT/results/20200726_160055_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_lateralWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_lateralreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorWT) <- c("4pe5_reslist", "receptor_freq")
dim_lateralligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandWT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorWT <- rbind.fill(dim_lateralreceptorWT, pereslist)
dim_lateralligandWT <- rbind.fill(dim_lateralligandWT, naoreslist)

dim_lateralligandWT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorWT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorWT[dim_lateralreceptorWT$`4pe5_reslist` == i, ]$receptor_freq <- dim_lateralreceptorWT[dim_lateralreceptorWT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandWT$`2nao_reslist` <- factor(dim_lateralligandWT$`2nao_reslist`, levels = dim_lateralligandWT$`2nao_reslist`)

dim_lateralligandWT$rel_freq <- dim_lateralligandWT$ligand_freq / max(dim_lateralligandWT$ligand_freq)
dim_lateralreceptorWT$rel_freq <- dim_lateralreceptorWT$receptor_freq / max(dim_lateralreceptorWT$receptor_freq)

#plot WT ligand
dim_lateralligandWT$color <- dim_lateralligandWT$`2nao_reslist` == "ALA2U" | dim_lateralligandWT$`2nao_reslist` == "ALA2X"

ggplot(data = dim_lateralligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Lateral Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralreceptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Lateral Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Lateral Dimer Ligand Interaction")

ggplot(data = dim_lateralreceptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Lateral Dimer Receptor Interaction")
#save file for b-factor coloring

write.csv(dim_lateralligandWT, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralligandWT.csv")

write.csv(dim_lateralreceptorWT, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralreceptorWT.csv")

##dim_stacked
#A2T
rm(list = ls())

path <- "g:/research/belfortlocalfiles/zdock output/transfer/dim_stackedA2T/results/20200726_161034_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_stackedA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_stackedreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorAT) <- c("4pe5_reslist", "receptor_freq")
dim_stackedligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedligandAT) <- c("2nao_reslist", "ligand_freq")

dim_stackedreceptorAT <- rbind.fill(dim_stackedreceptorAT, pereslist)
dim_stackedligandAT <- rbind.fill(dim_stackedligandAT, naoreslist)

dim_stackedligandAT$ligand_freq <- 0
for (i in llist){
  dim_stackedligandAT[dim_stackedligandAT$`2nao_reslist` == i, ]$ligand_freq <- dim_stackedligandAT[dim_stackedligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_stackedreceptorAT$receptor_freq <- 0
for (i in rlist){
  dim_stackedreceptorAT[dim_stackedreceptorAT$`4pe5_reslist` == i, ]$receptor_freq <- dim_stackedreceptorAT[dim_stackedreceptorAT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandAT$`2nao_reslist` <- factor(dim_stackedligandAT$`2nao_reslist`, levels = dim_stackedligandAT$`2nao_reslist`)

dim_stackedligandAT$rel_freq <- dim_stackedligandAT$ligand_freq / max(dim_stackedligandAT$ligand_freq)
dim_stackedreceptorAT$rel_freq <- dim_stackedreceptorAT$receptor_freq / max(dim_stackedreceptorAT$receptor_freq)

#plot A2T ligand
dim_stackedligandAT$color <- dim_stackedligandAT$`2nao_reslist` == "THR2U" | dim_stackedligandAT$`2nao_reslist` == "THR2V"

ggplot(data = dim_stackedligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Stacked Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedreceptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Stacked Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Stacked Dimer Ligand Interaction")

ggplot(data = dim_stackedreceptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Stacked Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_stackedligandAT, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedligandAT.csv")

write.csv(dim_stackedreceptorAT, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedreceptorAT.csv")

#A2V
rm(list = ls())

path <- "g:/research/belfortlocalfiles/zdock output/transfer/dim_stackedA2V/results/20200726_162053_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_stackedA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_stackedreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorAV) <- c("4pe5_reslist", "receptor_freq")
dim_stackedligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedligandAV) <- c("2nao_reslist", "ligand_freq")

dim_stackedreceptorAV <- rbind.fill(dim_stackedreceptorAV, pereslist)
dim_stackedligandAV <- rbind.fill(dim_stackedligandAV, naoreslist)

dim_stackedligandAV$ligand_freq <- 0
for (i in llist){
  dim_stackedligandAV[dim_stackedligandAV$`2nao_reslist` == i, ]$ligand_freq <- dim_stackedligandAV[dim_stackedligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_stackedreceptorAV$receptor_freq <- 0
for (i in rlist){
  dim_stackedreceptorAV[dim_stackedreceptorAV$`4pe5_reslist` == i, ]$receptor_freq <- dim_stackedreceptorAV[dim_stackedreceptorAV$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandAV$`2nao_reslist` <- factor(dim_stackedligandAV$`2nao_reslist`, levels = dim_stackedligandAV$`2nao_reslist`)

dim_stackedligandAV$rel_freq <- dim_stackedligandAV$ligand_freq / max(dim_stackedligandAV$ligand_freq)
dim_stackedreceptorAV$rel_freq <- dim_stackedreceptorAV$receptor_freq / max(dim_stackedreceptorAV$receptor_freq)

#plot A2V ligand
dim_stackedligandAV$color <- dim_stackedligandAV$`2nao_reslist` == "VAL2U" | dim_stackedligandAV$`2nao_reslist` == "VAL2V"

ggplot(data = dim_stackedligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Stacked Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedreceptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Stacked Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Stacked Dimer Ligand Interaction")

ggplot(data = dim_stackedreceptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Stacked Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_stackedligandAV, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedligandAV.csv")

write.csv(dim_stackedreceptorAV, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedreceptorAV.csv")

#WT
rm(list = ls())

path <- "g:/research/belfortlocalfiles/zdock output/transfer/dim_stackedWT/results/20200726_163036_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_stackedWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_stackedreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorWT) <- c("4pe5_reslist", "receptor_freq")
dim_stackedligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedligandWT) <- c("2nao_reslist", "ligand_freq")

dim_stackedreceptorWT <- rbind.fill(dim_stackedreceptorWT, pereslist)
dim_stackedligandWT <- rbind.fill(dim_stackedligandWT, naoreslist)

dim_stackedligandWT$ligand_freq <- 0
for (i in llist){
  dim_stackedligandWT[dim_stackedligandWT$`2nao_reslist` == i, ]$ligand_freq <- dim_stackedligandWT[dim_stackedligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_stackedreceptorWT$receptor_freq <- 0
for (i in rlist){
  dim_stackedreceptorWT[dim_stackedreceptorWT$`4pe5_reslist` == i, ]$receptor_freq <- dim_stackedreceptorWT[dim_stackedreceptorWT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandWT$`2nao_reslist` <- factor(dim_stackedligandWT$`2nao_reslist`, levels = dim_stackedligandWT$`2nao_reslist`)

dim_stackedligandWT$rel_freq <- dim_stackedligandWT$ligand_freq / max(dim_stackedligandWT$ligand_freq)
dim_stackedreceptorWT$rel_freq <- dim_stackedreceptorWT$receptor_freq / max(dim_stackedreceptorWT$receptor_freq)

#plot WT ligand
dim_stackedligandWT$color <- dim_stackedligandWT$`2nao_reslist` == "ALA2U" | dim_stackedligandWT$`2nao_reslist` == "ALA2V"

ggplot(data = dim_stackedligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Stacked Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedreceptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Stacked Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Stacked Dimer Ligand Interaction")

ggplot(data = dim_stackedreceptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Stacked Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_stackedligandWT, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedligandWT.csv")

write.csv(dim_stackedreceptorWT, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedreceptorWT.csv")

##mono
#A2T
rm(list = ls())

path <- "g:/research/belfortlocalfiles/zdock output/transfer/mono_A2T/results/20200727_150906_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/mono_A2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

mono_receptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorAT) <- c("4pe5_reslist", "receptor_freq")
mono_ligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_ligandAT) <- c("2nao_reslist", "ligand_freq")

mono_receptorAT <- rbind.fill(mono_receptorAT, pereslist)
mono_ligandAT <- rbind.fill(mono_ligandAT, naoreslist)

mono_ligandAT$ligand_freq <- 0
for (i in llist){
  mono_ligandAT[mono_ligandAT$`2nao_reslist` == i, ]$ligand_freq <- mono_ligandAT[mono_ligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

mono_receptorAT$receptor_freq <- 0
for (i in rlist){
  mono_receptorAT[mono_receptorAT$`4pe5_reslist` == i, ]$receptor_freq <- mono_receptorAT[mono_receptorAT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

mono_ligandAT$`2nao_reslist` <- factor(mono_ligandAT$`2nao_reslist`, levels = mono_ligandAT$`2nao_reslist`)

mono_ligandAT$rel_freq <- mono_ligandAT$ligand_freq / max(mono_ligandAT$ligand_freq)
mono_receptorAT$rel_freq <- mono_receptorAT$receptor_freq / max(mono_receptorAT$receptor_freq)

#plot A2T ligand
mono_ligandAT$color <- mono_ligandAT$`2nao_reslist` == "THR2U"

ggplot(data = mono_ligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Monomer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_receptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Monomer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_ligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Monomer Ligand Interaction")

ggplot(data = mono_receptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Monomer Receptor Interaction")

#save file for b-factor coloring

write.csv(mono_ligandAT, "g:/research/belfortlocalfiles/data/bfactor/mono_ligandAT.csv")

write.csv(mono_receptorAT, "g:/research/belfortlocalfiles/data/bfactor/mono_receptorAT.csv")

#A2V
rm(list = ls())

path <- "g:/research/belfortlocalfiles/zdock output/transfer/mono_A2V/results/20200726_164140_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/mono_A2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

mono_receptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorAV) <- c("4pe5_reslist", "receptor_freq")
mono_ligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_ligandAV) <- c("2nao_reslist", "ligand_freq")

mono_receptorAV <- rbind.fill(mono_receptorAV, pereslist)
mono_ligandAV <- rbind.fill(mono_ligandAV, naoreslist)

mono_ligandAV$ligand_freq <- 0
for (i in llist){
  mono_ligandAV[mono_ligandAV$`2nao_reslist` == i, ]$ligand_freq <- mono_ligandAV[mono_ligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

mono_receptorAV$receptor_freq <- 0
for (i in rlist){
  mono_receptorAV[mono_receptorAV$`4pe5_reslist` == i, ]$receptor_freq <- mono_receptorAV[mono_receptorAV$`4pe5_reslist` == i, ]$receptor_freq + 1
}

mono_ligandAV$`2nao_reslist` <- factor(mono_ligandAV$`2nao_reslist`, levels = mono_ligandAV$`2nao_reslist`)

mono_ligandAV$rel_freq <- mono_ligandAV$ligand_freq / max(mono_ligandAV$ligand_freq)
mono_receptorAV$rel_freq <- mono_receptorAV$receptor_freq / max(mono_receptorAV$receptor_freq)

#plot A2V ligand
mono_ligandAV$color <- mono_ligandAV$`2nao_reslist` == "VAL2U"

ggplot(data = mono_ligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Monomer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_receptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Monomer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_ligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Monomer Ligand Interaction")

ggplot(data = mono_receptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Monomer Receptor Interaction")

#save file for b-factor coloring

write.csv(mono_ligandAV, "g:/research/belfortlocalfiles/data/bfactor/mono_ligandAV.csv")

write.csv(mono_receptorAV, "g:/research/belfortlocalfiles/data/bfactor/mono_receptorAV.csv")

#WT
rm(list = ls())

path <- "g:/research/belfortlocalfiles/zdock output/transfer/mono_WT/results/20200726_165135_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/mono_WT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

mono_receptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorWT) <- c("4pe5_reslist", "receptor_freq")
mono_ligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_ligandWT) <- c("2nao_reslist", "ligand_freq")

mono_receptorWT <- rbind.fill(mono_receptorWT, pereslist)
mono_ligandWT <- rbind.fill(mono_ligandWT, naoreslist)

mono_ligandWT$ligand_freq <- 0
for (i in llist){
  mono_ligandWT[mono_ligandWT$`2nao_reslist` == i, ]$ligand_freq <- mono_ligandWT[mono_ligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

mono_receptorWT$receptor_freq <- 0
for (i in rlist){
  mono_receptorWT[mono_receptorWT$`4pe5_reslist` == i, ]$receptor_freq <- mono_receptorWT[mono_receptorWT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

mono_ligandWT$`2nao_reslist` <- factor(mono_ligandWT$`2nao_reslist`, levels = mono_ligandWT$`2nao_reslist`)

mono_ligandWT$rel_freq <- mono_ligandWT$ligand_freq / max(mono_ligandWT$ligand_freq)
mono_receptorWT$rel_freq <- mono_receptorWT$receptor_freq / max(mono_receptorWT$receptor_freq)

#plot WT ligand
mono_ligandWT$color <- mono_ligandWT$`2nao_reslist` == "ALA2U"

ggplot(data = mono_ligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Monomer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_receptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Monomer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_ligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Monomer Ligand Interaction")

ggplot(data = mono_receptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Monomer Receptor Interaction")

#save file for b-factor coloring

write.csv(mono_ligandWT, "g:/research/belfortlocalfiles/data/bfactor/mono_ligandWT.csv")

write.csv(mono_receptorWT, "g:/research/belfortlocalfiles/data/bfactor/mono_receptorWT.csv")

##### 10-20-2020 #####
#Organize zdock runs for 1-200 models
#rerun data because I didn't transfer it to my laptop :(
#blocked 4pe5
library(plyr)
library(ggplot2)

rm(list = ls())

##dim_lateral
#A2T
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/dim_latAT/20201020_163436_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_lateralreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAT) <- c("4pe5_reslist", "receptor_freq")
dim_lateralligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAT <- rbind.fill(dim_lateralreceptorAT, pereslist)
dim_lateralligandAT <- rbind.fill(dim_lateralligandAT, naoreslist)

dim_lateralligandAT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAT[dim_lateralreceptorAT$`4pe5_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAT[dim_lateralreceptorAT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAT$`2nao_reslist` <- factor(dim_lateralligandAT$`2nao_reslist`, levels = dim_lateralligandAT$`2nao_reslist`)

dim_lateralligandAT$rel_freq <- dim_lateralligandAT$ligand_freq / max(dim_lateralligandAT$ligand_freq)
dim_lateralreceptorAT$rel_freq <- dim_lateralreceptorAT$receptor_freq / max(dim_lateralreceptorAT$receptor_freq)

#plot A2T ligand
dim_lateralligandAT$color <- dim_lateralligandAT$`2nao_reslist` == "THR2U" | dim_lateralligandAT$`2nao_reslist` == "THR2X"

ggplot(data = dim_lateralligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Lateral Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralreceptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Lateral Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Lateral Dimer Ligand Interaction")

ggplot(data = dim_lateralreceptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Lateral Dimer Receptor Interaction")

#save file for b-factor coloring (don't need atm but here is path name)

write.csv(dim_lateralligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/rerun/dim_latligAT.csv")

write.csv(dim_lateralreceptorAT, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralreceptorAT.csv")

#A2V
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/dim_latAV/20201020_165406_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

#need to change below paths if intending to run
naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_lateralA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_lateralreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAV) <- c("4pe5_reslist", "receptor_freq")
dim_lateralligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAV) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAV <- rbind.fill(dim_lateralreceptorAV, pereslist)
dim_lateralligandAV <- rbind.fill(dim_lateralligandAV, naoreslist)

dim_lateralligandAV$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAV$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAV[dim_lateralreceptorAV$`4pe5_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAV[dim_lateralreceptorAV$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAV$`2nao_reslist` <- factor(dim_lateralligandAV$`2nao_reslist`, levels = dim_lateralligandAV$`2nao_reslist`)

dim_lateralligandAV$rel_freq <- dim_lateralligandAV$ligand_freq / max(dim_lateralligandAV$ligand_freq)
dim_lateralreceptorAV$rel_freq <- dim_lateralreceptorAV$receptor_freq / max(dim_lateralreceptorAV$receptor_freq)

#plot A2V ligand
dim_lateralligandAV$color <- dim_lateralligandAV$`2nao_reslist` == "VAL2U" | dim_lateralligandAV$`2nao_reslist` == "VAL2X"

ggplot(data = dim_lateralligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Lateral Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralreceptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Lateral Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Lateral Dimer Ligand Interaction")

ggplot(data = dim_lateralreceptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Lateral Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_lateralligandAV, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralligandAV.csv")

write.csv(dim_lateralreceptorAV, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralreceptorAV.csv")

#WT
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/dim_latWT/20201020_170409_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_lateralWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_lateralreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorWT) <- c("4pe5_reslist", "receptor_freq")
dim_lateralligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandWT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorWT <- rbind.fill(dim_lateralreceptorWT, pereslist)
dim_lateralligandWT <- rbind.fill(dim_lateralligandWT, naoreslist)

dim_lateralligandWT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorWT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorWT[dim_lateralreceptorWT$`4pe5_reslist` == i, ]$receptor_freq <- dim_lateralreceptorWT[dim_lateralreceptorWT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandWT$`2nao_reslist` <- factor(dim_lateralligandWT$`2nao_reslist`, levels = dim_lateralligandWT$`2nao_reslist`)

dim_lateralligandWT$rel_freq <- dim_lateralligandWT$ligand_freq / max(dim_lateralligandWT$ligand_freq)
dim_lateralreceptorWT$rel_freq <- dim_lateralreceptorWT$receptor_freq / max(dim_lateralreceptorWT$receptor_freq)

#plot WT ligand
dim_lateralligandWT$color <- dim_lateralligandWT$`2nao_reslist` == "ALA2U" | dim_lateralligandWT$`2nao_reslist` == "ALA2X"

ggplot(data = dim_lateralligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Lateral Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralreceptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Lateral Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_lateralligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Lateral Dimer Ligand Interaction")

ggplot(data = dim_lateralreceptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Lateral Dimer Receptor Interaction")
#save file for b-factor coloring

write.csv(dim_lateralligandWT, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralligandWT.csv")

write.csv(dim_lateralreceptorWT, "g:/research/belfortlocalfiles/data/bfactor/dim_lateralreceptorWT.csv")

##dim_stacked
#A2T
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/dim_stackAT/20201020_171341_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_stackedA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_stackedreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorAT) <- c("4pe5_reslist", "receptor_freq")
dim_stackedligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedligandAT) <- c("2nao_reslist", "ligand_freq")

dim_stackedreceptorAT <- rbind.fill(dim_stackedreceptorAT, pereslist)
dim_stackedligandAT <- rbind.fill(dim_stackedligandAT, naoreslist)

dim_stackedligandAT$ligand_freq <- 0
for (i in llist){
  dim_stackedligandAT[dim_stackedligandAT$`2nao_reslist` == i, ]$ligand_freq <- dim_stackedligandAT[dim_stackedligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_stackedreceptorAT$receptor_freq <- 0
for (i in rlist){
  dim_stackedreceptorAT[dim_stackedreceptorAT$`4pe5_reslist` == i, ]$receptor_freq <- dim_stackedreceptorAT[dim_stackedreceptorAT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandAT$`2nao_reslist` <- factor(dim_stackedligandAT$`2nao_reslist`, levels = dim_stackedligandAT$`2nao_reslist`)

dim_stackedligandAT$rel_freq <- dim_stackedligandAT$ligand_freq / max(dim_stackedligandAT$ligand_freq)
dim_stackedreceptorAT$rel_freq <- dim_stackedreceptorAT$receptor_freq / max(dim_stackedreceptorAT$receptor_freq)

#plot A2T ligand
dim_stackedligandAT$color <- dim_stackedligandAT$`2nao_reslist` == "THR2U" | dim_stackedligandAT$`2nao_reslist` == "THR2V"

ggplot(data = dim_stackedligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Stacked Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedreceptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Stacked Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Stacked Dimer Ligand Interaction")

ggplot(data = dim_stackedreceptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Stacked Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_stackedligandAT, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedligandAT.csv")

write.csv(dim_stackedreceptorAT, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedreceptorAT.csv")

#A2V
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/dim_stackAV/20201020_172324_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_stackedA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_stackedreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorAV) <- c("4pe5_reslist", "receptor_freq")
dim_stackedligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedligandAV) <- c("2nao_reslist", "ligand_freq")

dim_stackedreceptorAV <- rbind.fill(dim_stackedreceptorAV, pereslist)
dim_stackedligandAV <- rbind.fill(dim_stackedligandAV, naoreslist)

dim_stackedligandAV$ligand_freq <- 0
for (i in llist){
  dim_stackedligandAV[dim_stackedligandAV$`2nao_reslist` == i, ]$ligand_freq <- dim_stackedligandAV[dim_stackedligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_stackedreceptorAV$receptor_freq <- 0
for (i in rlist){
  dim_stackedreceptorAV[dim_stackedreceptorAV$`4pe5_reslist` == i, ]$receptor_freq <- dim_stackedreceptorAV[dim_stackedreceptorAV$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandAV$`2nao_reslist` <- factor(dim_stackedligandAV$`2nao_reslist`, levels = dim_stackedligandAV$`2nao_reslist`)

dim_stackedligandAV$rel_freq <- dim_stackedligandAV$ligand_freq / max(dim_stackedligandAV$ligand_freq)
dim_stackedreceptorAV$rel_freq <- dim_stackedreceptorAV$receptor_freq / max(dim_stackedreceptorAV$receptor_freq)

#plot A2V ligand
dim_stackedligandAV$color <- dim_stackedligandAV$`2nao_reslist` == "VAL2U" | dim_stackedligandAV$`2nao_reslist` == "VAL2V"

ggplot(data = dim_stackedligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Stacked Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedreceptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Stacked Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Stacked Dimer Ligand Interaction")

ggplot(data = dim_stackedreceptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Stacked Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_stackedligandAV, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedligandAV.csv")

write.csv(dim_stackedreceptorAV, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedreceptorAV.csv")

#WT
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/dim_stackWT/20201020_175704_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/dim_stackedWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

dim_stackedreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorWT) <- c("4pe5_reslist", "receptor_freq")
dim_stackedligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedligandWT) <- c("2nao_reslist", "ligand_freq")

dim_stackedreceptorWT <- rbind.fill(dim_stackedreceptorWT, pereslist)
dim_stackedligandWT <- rbind.fill(dim_stackedligandWT, naoreslist)

dim_stackedligandWT$ligand_freq <- 0
for (i in llist){
  dim_stackedligandWT[dim_stackedligandWT$`2nao_reslist` == i, ]$ligand_freq <- dim_stackedligandWT[dim_stackedligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_stackedreceptorWT$receptor_freq <- 0
for (i in rlist){
  dim_stackedreceptorWT[dim_stackedreceptorWT$`4pe5_reslist` == i, ]$receptor_freq <- dim_stackedreceptorWT[dim_stackedreceptorWT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandWT$`2nao_reslist` <- factor(dim_stackedligandWT$`2nao_reslist`, levels = dim_stackedligandWT$`2nao_reslist`)

dim_stackedligandWT$rel_freq <- dim_stackedligandWT$ligand_freq / max(dim_stackedligandWT$ligand_freq)
dim_stackedreceptorWT$rel_freq <- dim_stackedreceptorWT$receptor_freq / max(dim_stackedreceptorWT$receptor_freq)

#plot WT ligand
dim_stackedligandWT$color <- dim_stackedligandWT$`2nao_reslist` == "ALA2U" | dim_stackedligandWT$`2nao_reslist` == "ALA2V"

ggplot(data = dim_stackedligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Stacked Dimer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedreceptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Stacked Dimer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = dim_stackedligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Stacked Dimer Ligand Interaction")

ggplot(data = dim_stackedreceptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Stacked Dimer Receptor Interaction")

#save file for b-factor coloring

write.csv(dim_stackedligandWT, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedligandWT.csv")

write.csv(dim_stackedreceptorWT, "g:/research/belfortlocalfiles/data/bfactor/dim_stackedreceptorWT.csv")

##mono
#A2T
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/mono_AT/20201020_180648_output.txt"

temp <- read.csv(path, skip = 6, header = F, sep = "", stringsAsFactors = F) #odd fix skipping columns (all invalid anyway)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/mono_A2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

mono_receptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorAT) <- c("4pe5_reslist", "receptor_freq")
mono_ligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_ligandAT) <- c("2nao_reslist", "ligand_freq")

mono_receptorAT <- rbind.fill(mono_receptorAT, pereslist)
mono_ligandAT <- rbind.fill(mono_ligandAT, naoreslist)

mono_ligandAT$ligand_freq <- 0
for (i in llist){
  mono_ligandAT[mono_ligandAT$`2nao_reslist` == i, ]$ligand_freq <- mono_ligandAT[mono_ligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

mono_receptorAT$receptor_freq <- 0
for (i in rlist){
  mono_receptorAT[mono_receptorAT$`4pe5_reslist` == i, ]$receptor_freq <- mono_receptorAT[mono_receptorAT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

mono_ligandAT$`2nao_reslist` <- factor(mono_ligandAT$`2nao_reslist`, levels = mono_ligandAT$`2nao_reslist`)

mono_ligandAT$rel_freq <- mono_ligandAT$ligand_freq / max(mono_ligandAT$ligand_freq)
mono_receptorAT$rel_freq <- mono_receptorAT$receptor_freq / max(mono_receptorAT$receptor_freq)

#plot A2T ligand
mono_ligandAT$color <- mono_ligandAT$`2nao_reslist` == "THR2U"

ggplot(data = mono_ligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Monomer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_receptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Monomer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_ligandAT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Monomer Ligand Interaction")

ggplot(data = mono_receptorAT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2T Monomer Receptor Interaction")

#save file for b-factor coloring

write.csv(mono_ligandAT, "g:/research/belfortlocalfiles/data/bfactor/mono_ligandAT.csv")

write.csv(mono_receptorAT, "g:/research/belfortlocalfiles/data/bfactor/mono_receptorAT.csv")

#A2V
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/mono_AV/20201020_181729_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/mono_A2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

mono_receptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorAV) <- c("4pe5_reslist", "receptor_freq")
mono_ligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_ligandAV) <- c("2nao_reslist", "ligand_freq")

mono_receptorAV <- rbind.fill(mono_receptorAV, pereslist)
mono_ligandAV <- rbind.fill(mono_ligandAV, naoreslist)

mono_ligandAV$ligand_freq <- 0
for (i in llist){
  mono_ligandAV[mono_ligandAV$`2nao_reslist` == i, ]$ligand_freq <- mono_ligandAV[mono_ligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

mono_receptorAV$receptor_freq <- 0
for (i in rlist){
  mono_receptorAV[mono_receptorAV$`4pe5_reslist` == i, ]$receptor_freq <- mono_receptorAV[mono_receptorAV$`4pe5_reslist` == i, ]$receptor_freq + 1
}

mono_ligandAV$`2nao_reslist` <- factor(mono_ligandAV$`2nao_reslist`, levels = mono_ligandAV$`2nao_reslist`)

mono_ligandAV$rel_freq <- mono_ligandAV$ligand_freq / max(mono_ligandAV$ligand_freq)
mono_receptorAV$rel_freq <- mono_receptorAV$receptor_freq / max(mono_receptorAV$receptor_freq)

#plot A2V ligand
mono_ligandAV$color <- mono_ligandAV$`2nao_reslist` == "VAL2U"

ggplot(data = mono_ligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Monomer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_receptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Monomer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_ligandAV) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Monomer Ligand Interaction")

ggplot(data = mono_receptorAV) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "A2V Monomer Receptor Interaction")

#save file for b-factor coloring

write.csv(mono_ligandAV, "g:/research/belfortlocalfiles/data/bfactor/mono_ligandAV.csv")

write.csv(mono_receptorAV, "g:/research/belfortlocalfiles/data/bfactor/mono_receptorAV.csv")

#WT
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked/mono_WT/20201020_183103_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp)

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/mono_WT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("g:/research/belfortlocalfiles/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "4pe5_reslist"

mono_receptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorWT) <- c("4pe5_reslist", "receptor_freq")
mono_ligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_ligandWT) <- c("2nao_reslist", "ligand_freq")

mono_receptorWT <- rbind.fill(mono_receptorWT, pereslist)
mono_ligandWT <- rbind.fill(mono_ligandWT, naoreslist)

mono_ligandWT$ligand_freq <- 0
for (i in llist){
  mono_ligandWT[mono_ligandWT$`2nao_reslist` == i, ]$ligand_freq <- mono_ligandWT[mono_ligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

mono_receptorWT$receptor_freq <- 0
for (i in rlist){
  mono_receptorWT[mono_receptorWT$`4pe5_reslist` == i, ]$receptor_freq <- mono_receptorWT[mono_receptorWT$`4pe5_reslist` == i, ]$receptor_freq + 1
}

mono_ligandWT$`2nao_reslist` <- factor(mono_ligandWT$`2nao_reslist`, levels = mono_ligandWT$`2nao_reslist`)

mono_ligandWT$rel_freq <- mono_ligandWT$ligand_freq / max(mono_ligandWT$ligand_freq)
mono_receptorWT$rel_freq <- mono_receptorWT$receptor_freq / max(mono_receptorWT$receptor_freq)

#plot WT ligand
mono_ligandWT$color <- mono_ligandWT$`2nao_reslist` == "ALA2U"

ggplot(data = mono_ligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = rel_freq, fill = color), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Monomer Ligand Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_receptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = rel_freq), show.legend = F) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Monomer Receptor Interaction", x = "Residue ID", y = "Relative Frequency")

ggplot(data = mono_ligandWT) +
  geom_col(aes(x = `2nao_reslist`, y = ligand_freq, fill = color)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Monomer Ligand Interaction")

ggplot(data = mono_receptorWT) +
  geom_col(aes(x = `4pe5_reslist`, y = receptor_freq)) +
  theme(axis.text.x = element_text(angle = 270)) +
  labs(title = "WT Monomer Receptor Interaction")

#save file for b-factor coloring

write.csv(mono_ligandWT, "g:/research/belfortlocalfiles/data/bfactor/mono_ligandWT.csv")

write.csv(mono_receptorWT, "g:/research/belfortlocalfiles/data/bfactor/mono_receptorWT.csv")

##### 11-08-2020 #####
#Organize zdock runs for 1-200 models
#thing to note: it was easier for me to copy paste the code blocks and change minimal details for each section for that reason:
#data frames are all named dim_lateral* or dim_receptor*. does not change the final exported .csv files which are named individually
#blocked 6IRA zdock data
#blocked 6IRA
library(plyr)
library(ggplot2)

rm(list = ls())

##dim_lateral
#A2T
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/dim_lateralA2T/20201101_094231_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAT) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAT <- rbind.fill(dim_lateralreceptorAT, irareslist)
dim_lateralligandAT <- rbind.fill(dim_lateralligandAT, naoreslist)

dim_lateralligandAT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAT[dim_lateralreceptorAT$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAT[dim_lateralreceptorAT$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAT$`2nao_reslist` <- factor(dim_lateralligandAT$`2nao_reslist`, levels = dim_lateralligandAT$`2nao_reslist`)

dim_lateralligandAT$rel_freq <- dim_lateralligandAT$ligand_freq / max(dim_lateralligandAT$ligand_freq)
dim_lateralreceptorAT$rel_freq <- dim_lateralreceptorAT$receptor_freq / max(dim_lateralreceptorAT$receptor_freq)

#Add color column
dim_lateralligandAT$color <- dim_lateralligandAT$`2nao_reslist` == "THR2U" | dim_lateralligandAT$`2nao_reslist` == "THR2X"


#save file for b-factor coloring

write.csv(dim_lateralligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandAT.csv")

write.csv(dim_lateralreceptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralreceptorAT.csv")

#A2V
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/dim_lateralA2V/20201101_095526_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAV) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAV) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAV <- rbind.fill(dim_lateralreceptorAV, irareslist)
dim_lateralligandAV <- rbind.fill(dim_lateralligandAV, naoreslist)

dim_lateralligandAV$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAV$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAV[dim_lateralreceptorAV$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAV[dim_lateralreceptorAV$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAV$`2nao_reslist` <- factor(dim_lateralligandAV$`2nao_reslist`, levels = dim_lateralligandAV$`2nao_reslist`)

dim_lateralligandAV$rel_freq <- dim_lateralligandAV$ligand_freq / max(dim_lateralligandAV$ligand_freq)
dim_lateralreceptorAV$rel_freq <- dim_lateralreceptorAV$receptor_freq / max(dim_lateralreceptorAV$receptor_freq)

#Add color column
dim_lateralligandAV$color <- dim_lateralligandAV$`2nao_reslist` == "VAL2U" | dim_lateralligandAV$`2nao_reslist` == "VAL2X"


#save file for b-factor coloring

write.csv(dim_lateralligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandAV.csv")

write.csv(dim_lateralreceptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralreceptorAV.csv")

#WT
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/dim_lateralWT/20201101_100810_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorWT) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandWT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorWT <- rbind.fill(dim_lateralreceptorWT, irareslist)
dim_lateralligandWT <- rbind.fill(dim_lateralligandWT, naoreslist)

dim_lateralligandWT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorWT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorWT[dim_lateralreceptorWT$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorWT[dim_lateralreceptorWT$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandWT$`2nao_reslist` <- factor(dim_lateralligandWT$`2nao_reslist`, levels = dim_lateralligandWT$`2nao_reslist`)

dim_lateralligandWT$rel_freq <- dim_lateralligandWT$ligand_freq / max(dim_lateralligandWT$ligand_freq)
dim_lateralreceptorWT$rel_freq <- dim_lateralreceptorWT$receptor_freq / max(dim_lateralreceptorWT$receptor_freq)

#Add color column
dim_lateralligandWT$color <- dim_lateralligandWT$`2nao_reslist` == "ALA2U" | dim_lateralligandWT$`2nao_reslist` == "ALA2X"


#save file for b-factor coloring

write.csv(dim_lateralligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandWT.csv")

write.csv(dim_lateralreceptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralreceptorWT.csv")

##dim_stacked
#A2T
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/dim_stackedA2T/20201101_101909_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAT) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAT <- rbind.fill(dim_lateralreceptorAT, irareslist)
dim_lateralligandAT <- rbind.fill(dim_lateralligandAT, naoreslist)

dim_lateralligandAT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAT[dim_lateralreceptorAT$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAT[dim_lateralreceptorAT$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAT$`2nao_reslist` <- factor(dim_lateralligandAT$`2nao_reslist`, levels = dim_lateralligandAT$`2nao_reslist`)

dim_lateralligandAT$rel_freq <- dim_lateralligandAT$ligand_freq / max(dim_lateralligandAT$ligand_freq)
dim_lateralreceptorAT$rel_freq <- dim_lateralreceptorAT$receptor_freq / max(dim_lateralreceptorAT$receptor_freq)

#Add color column
dim_lateralligandAT$color <- dim_lateralligandAT$`2nao_reslist` == "THR2U" | dim_lateralligandAT$`2nao_reslist` == "THR2V"


#save file for b-factor coloring

write.csv(dim_lateralligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandAT.csv")

write.csv(dim_lateralreceptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedreceptorAT.csv")

#A2V
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/dim_stackedA2V/20201101_103048_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAV) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAV) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAV <- rbind.fill(dim_lateralreceptorAV, irareslist)
dim_lateralligandAV <- rbind.fill(dim_lateralligandAV, naoreslist)

dim_lateralligandAV$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAV$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAV[dim_lateralreceptorAV$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAV[dim_lateralreceptorAV$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAV$`2nao_reslist` <- factor(dim_lateralligandAV$`2nao_reslist`, levels = dim_lateralligandAV$`2nao_reslist`)

dim_lateralligandAV$rel_freq <- dim_lateralligandAV$ligand_freq / max(dim_lateralligandAV$ligand_freq)
dim_lateralreceptorAV$rel_freq <- dim_lateralreceptorAV$receptor_freq / max(dim_lateralreceptorAV$receptor_freq)

#Add color column
dim_lateralligandAV$color <- dim_lateralligandAV$`2nao_reslist` == "VAL2U" | dim_lateralligandAV$`2nao_reslist` == "VAL2V"


#save file for b-factor coloring

write.csv(dim_lateralligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandAV.csv")

write.csv(dim_lateralreceptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedreceptorAV.csv")

#WT
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/dim_stackedWT/20201102_110653_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorWT) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandWT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorWT <- rbind.fill(dim_lateralreceptorWT, irareslist)
dim_lateralligandWT <- rbind.fill(dim_lateralligandWT, naoreslist)

dim_lateralligandWT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorWT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorWT[dim_lateralreceptorWT$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorWT[dim_lateralreceptorWT$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandWT$`2nao_reslist` <- factor(dim_lateralligandWT$`2nao_reslist`, levels = dim_lateralligandWT$`2nao_reslist`)

dim_lateralligandWT$rel_freq <- dim_lateralligandWT$ligand_freq / max(dim_lateralligandWT$ligand_freq)
dim_lateralreceptorWT$rel_freq <- dim_lateralreceptorWT$receptor_freq / max(dim_lateralreceptorWT$receptor_freq)

#Add color column
dim_lateralligandWT$color <- dim_lateralligandWT$`2nao_reslist` == "ALA2U" | dim_lateralligandWT$`2nao_reslist` == "ALA2V"


#save file for b-factor coloring

write.csv(dim_lateralligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandWT.csv")

write.csv(dim_lateralreceptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedreceptorWT.csv")

##mono
#A2T
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/mono_A2T/20201102_111950_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_A2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAT) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAT <- rbind.fill(dim_lateralreceptorAT, irareslist)
dim_lateralligandAT <- rbind.fill(dim_lateralligandAT, naoreslist)

dim_lateralligandAT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAT[dim_lateralligandAT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAT[dim_lateralreceptorAT$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAT[dim_lateralreceptorAT$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAT$`2nao_reslist` <- factor(dim_lateralligandAT$`2nao_reslist`, levels = dim_lateralligandAT$`2nao_reslist`)

dim_lateralligandAT$rel_freq <- dim_lateralligandAT$ligand_freq / max(dim_lateralligandAT$ligand_freq)
dim_lateralreceptorAT$rel_freq <- dim_lateralreceptorAT$receptor_freq / max(dim_lateralreceptorAT$receptor_freq)

#Add color column
dim_lateralligandAT$color <- dim_lateralligandAT$`2nao_reslist` == "THR2U"


#save file for b-factor coloring

write.csv(dim_lateralligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandAT.csv")

write.csv(dim_lateralreceptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_receptorAT.csv")

#A2V
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/mono_A2V/20201102_113201_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_A2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAV) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAV) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorAV <- rbind.fill(dim_lateralreceptorAV, irareslist)
dim_lateralligandAV <- rbind.fill(dim_lateralligandAV, naoreslist)

dim_lateralligandAV$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAV$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAV[dim_lateralreceptorAV$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAV[dim_lateralreceptorAV$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAV$`2nao_reslist` <- factor(dim_lateralligandAV$`2nao_reslist`, levels = dim_lateralligandAV$`2nao_reslist`)

dim_lateralligandAV$rel_freq <- dim_lateralligandAV$ligand_freq / max(dim_lateralligandAV$ligand_freq)
dim_lateralreceptorAV$rel_freq <- dim_lateralreceptorAV$receptor_freq / max(dim_lateralreceptorAV$receptor_freq)

#Add color column
dim_lateralligandAV$color <- dim_lateralligandAV$`2nao_reslist` == "VAL2U"


#save file for b-factor coloring

write.csv(dim_lateralligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandAV.csv")

write.csv(dim_lateralreceptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_receptorAV.csv")

#WT
rm(list = ls())

path <- "c:/users/wessec/documents/research/rpi gb/zdockoutput/blocked_6ira/mono_WT/20201102_114845_output.txt"

temp <- read.csv(path, skip = 1, header = F, sep = "", stringsAsFactors = F)
temp[, 2] <- gsub("invalid", NA, temp[, 2])
temp <- na.omit(temp) #see how many invalids here

rlist <- NULL
llist <- NULL
for (j in 1:nrow(temp)){
  resextract <- strsplit(temp[j, 2], split = ",")
  rlist <- c(rlist, resextract[[1]])
  ligextract <- strsplit(temp[j, 3], split = ",")
  llist <- c(llist, ligextract[[1]])
}
rm(ligextract, resextract)

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_WT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
irareslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6ira_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(irareslist) <- "6ira_reslist"

dim_lateralreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorWT) <- c("6ira_reslist", "receptor_freq")
dim_lateralligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandWT) <- c("2nao_reslist", "ligand_freq")

dim_lateralreceptorWT <- rbind.fill(dim_lateralreceptorWT, irareslist)
dim_lateralligandWT <- rbind.fill(dim_lateralligandWT, naoreslist)

dim_lateralligandWT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorWT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorWT[dim_lateralreceptorWT$`6ira_reslist` == i, ]$receptor_freq <- dim_lateralreceptorWT[dim_lateralreceptorWT$`6ira_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandWT$`2nao_reslist` <- factor(dim_lateralligandWT$`2nao_reslist`, levels = dim_lateralligandWT$`2nao_reslist`)

dim_lateralligandWT$rel_freq <- dim_lateralligandWT$ligand_freq / max(dim_lateralligandWT$ligand_freq)
dim_lateralreceptorWT$rel_freq <- dim_lateralreceptorWT$receptor_freq / max(dim_lateralreceptorWT$receptor_freq)

#Add color column
dim_lateralligandWT$color <- dim_lateralligandWT$`2nao_reslist` == "ALA2U"


#save file for b-factor coloring

write.csv(dim_lateralligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandWT.csv")

write.csv(dim_lateralreceptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_receptorWT.csv")
