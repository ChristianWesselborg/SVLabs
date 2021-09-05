#Organize cluspro runs
library(plyr)
library(ggplot2)

rm(list = ls())

##### 4pe5 runs#####
##blocked
##dim_lateral
#A2T
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "dimlat_AT/20200906_163018_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

#save file for b-factor coloring

write.csv(dim_lateralligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_lateralligandAT.csv")

write.csv(dim_lateralreceptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_lateralreceptorAT.csv")

#A2V
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "dimlat_AV/20200906_163701_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
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

dim_lateralligandAV$color <- dim_lateralligandAV$`2nao_reslist` == "VAL2U" | dim_lateralligandAV$`2nao_reslist` == "VAL2X"

#save file for b-factor coloring

write.csv(dim_lateralligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_lateralligandAV.csv")

write.csv(dim_lateralreceptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_lateralreceptorAV.csv")

#WT
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "dimlat_WT/20200906_163926_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
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

dim_lateralligandWT$color <- dim_lateralligandWT$`2nao_reslist` == "ALA2U" | dim_lateralligandWT$`2nao_reslist` == "ALA2X"

#save file for b-factor coloring

write.csv(dim_lateralligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_lateralligandWT.csv")

write.csv(dim_lateralreceptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_lateralreceptorWT.csv")

##dim_stacked
#A2T
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "dimstack_AT/20200906_164128_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
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

dim_stackedligandAT$color <- dim_stackedligandAT$`2nao_reslist` == "THR2U" | dim_stackedligandAT$`2nao_reslist` == "THR2V"

#save file for b-factor coloring

write.csv(dim_stackedligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_stackedligandAT.csv")

write.csv(dim_stackedreceptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_stackedreceptorAT.csv")

#A2V
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "dimstack_AV/20200906_164634_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
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

dim_stackedligandAV$color <- dim_stackedligandAV$`2nao_reslist` == "VAL2U" | dim_stackedligandAV$`2nao_reslist` == "VAL2V"

#save file for b-factor coloring

write.csv(dim_stackedligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_stackedligandAV.csv")

write.csv(dim_stackedreceptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_stackedreceptorAV.csv")

#WT
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "dimstack_WT/20200906_164821_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
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

dim_stackedligandWT$color <- dim_stackedligandWT$`2nao_reslist` == "ALA2U" | dim_stackedligandWT$`2nao_reslist` == "ALA2V"

#save file for b-factor coloring

write.csv(dim_stackedligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_stackedligandWT.csv")

write.csv(dim_stackedreceptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_dim_stackedreceptorWT.csv")

##mono
#A2T
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "mono_AT/20200906_165007_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_A2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
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

mono_ligandAT$color <- mono_ligandAT$`2nao_reslist` == "THR2U"

#save file for b-factor coloring

write.csv(mono_ligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_mono_ligandAT.csv")

write.csv(mono_receptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_mono_receptorAT.csv")

#A2V
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "mono_AV/20200906_165201_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_A2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
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

mono_ligandAV$color <- mono_ligandAV$`2nao_reslist` == "VAL2U"

#save file for b-factor coloring

write.csv(mono_ligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_mono_ligandAV.csv")

write.csv(mono_receptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_mono_receptorAV.csv")

#WT
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/blocked/"

basename <- "mono_WT/20200906_165351_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_WT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/4pe5_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
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

mono_ligandWT$color <- mono_ligandWT$`2nao_reslist` == "ALA2U"

#save file for b-factor coloring

write.csv(mono_ligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_mono_ligandWT.csv")

write.csv(mono_receptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/redo/clus_mono_receptorWT.csv")

##### 6IRA runs #####
##dim_lateral
#A2T
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "dimlat_AT/20201004_212246_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

dim_lateralreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAT) <- c("6IRA_reslist", "receptor_freq")
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
  dim_lateralreceptorAT[dim_lateralreceptorAT$`6IRA_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAT[dim_lateralreceptorAT$`6IRA_reslist` == i, ]$receptor_freq + 1
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

write.csv(dim_lateralligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_lateralligandAT.csv")

write.csv(dim_lateralreceptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_lateralreceptorAT.csv")

#A2V
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "dimlat_AV/20201004_213007_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

dim_lateralreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorAV) <- c("6IRA_reslist", "receptor_freq")
dim_lateralligandAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandAV) <- c("6IRA_reslist", "ligand_freq")

dim_lateralreceptorAV <- rbind.fill(dim_lateralreceptorAV, pereslist)
dim_lateralligandAV <- rbind.fill(dim_lateralligandAV, naoreslist)

dim_lateralligandAV$ligand_freq <- 0
for (i in llist){
  dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandAV[dim_lateralligandAV$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorAV$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorAV[dim_lateralreceptorAV$`6IRA_reslist` == i, ]$receptor_freq <- dim_lateralreceptorAV[dim_lateralreceptorAV$`6IRA_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandAV$`2nao_reslist` <- factor(dim_lateralligandAV$`2nao_reslist`, levels = dim_lateralligandAV$`2nao_reslist`)

dim_lateralligandAV$rel_freq <- dim_lateralligandAV$ligand_freq / max(dim_lateralligandAV$ligand_freq)
dim_lateralreceptorAV$rel_freq <- dim_lateralreceptorAV$receptor_freq / max(dim_lateralreceptorAV$receptor_freq)

dim_lateralligandAV$color <- dim_lateralligandAV$`2nao_reslist` == "VAL2U" | dim_lateralligandAV$`2nao_reslist` == "VAL2X"

#save file for b-factor coloring

write.csv(dim_lateralligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_lateralligandAV.csv")

write.csv(dim_lateralreceptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_lateralreceptorAV.csv")

#WT
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "dimlat_WT/20201004_213131_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_lateralWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

dim_lateralreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralreceptorWT) <- c("6IRA_reslist", "receptor_freq")
dim_lateralligandWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_lateralligandWT) <- c("6IRA_reslist", "ligand_freq")

dim_lateralreceptorWT <- rbind.fill(dim_lateralreceptorWT, pereslist)
dim_lateralligandWT <- rbind.fill(dim_lateralligandWT, naoreslist)

dim_lateralligandWT$ligand_freq <- 0
for (i in llist){
  dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq <- dim_lateralligandWT[dim_lateralligandWT$`2nao_reslist` == i, ]$ligand_freq + 1
}

dim_lateralreceptorWT$receptor_freq <- 0
for (i in rlist){
  dim_lateralreceptorWT[dim_lateralreceptorWT$`6IRA_reslist` == i, ]$receptor_freq <- dim_lateralreceptorWT[dim_lateralreceptorWT$`6IRA_reslist` == i, ]$receptor_freq + 1
}

dim_lateralligandWT$`2nao_reslist` <- factor(dim_lateralligandWT$`2nao_reslist`, levels = dim_lateralligandWT$`2nao_reslist`)

dim_lateralligandWT$rel_freq <- dim_lateralligandWT$ligand_freq / max(dim_lateralligandWT$ligand_freq)
dim_lateralreceptorWT$rel_freq <- dim_lateralreceptorWT$receptor_freq / max(dim_lateralreceptorWT$receptor_freq)

dim_lateralligandWT$color <- dim_lateralligandWT$`2nao_reslist` == "ALA2U" | dim_lateralligandWT$`2nao_reslist` == "ALA2X"

#save file for b-factor coloring

write.csv(dim_lateralligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_lateralligandWT.csv")

write.csv(dim_lateralreceptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_lateralreceptorWT.csv")

##dim_stacked
#A2T
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "dimstack_AT/20201004_213339_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedA2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

dim_stackedreceptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorAT) <- c("6IRA_reslist", "receptor_freq")
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
  dim_stackedreceptorAT[dim_stackedreceptorAT$`6IRA_reslist` == i, ]$receptor_freq <- dim_stackedreceptorAT[dim_stackedreceptorAT$`6IRA_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandAT$`2nao_reslist` <- factor(dim_stackedligandAT$`2nao_reslist`, levels = dim_stackedligandAT$`2nao_reslist`)

dim_stackedligandAT$rel_freq <- dim_stackedligandAT$ligand_freq / max(dim_stackedligandAT$ligand_freq)
dim_stackedreceptorAT$rel_freq <- dim_stackedreceptorAT$receptor_freq / max(dim_stackedreceptorAT$receptor_freq)

dim_stackedligandAT$color <- dim_stackedligandAT$`2nao_reslist` == "THR2U" | dim_stackedligandAT$`2nao_reslist` == "THR2V"

#save file for b-factor coloring

write.csv(dim_stackedligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_stackedligandAT.csv")

write.csv(dim_stackedreceptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_stackedreceptorAT.csv")

#A2V
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "dimstack_AV/20201004_213527_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedA2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

dim_stackedreceptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorAV) <- c("6IRA_reslist", "receptor_freq")
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
  dim_stackedreceptorAV[dim_stackedreceptorAV$`6IRA_reslist` == i, ]$receptor_freq <- dim_stackedreceptorAV[dim_stackedreceptorAV$`6IRA_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandAV$`2nao_reslist` <- factor(dim_stackedligandAV$`2nao_reslist`, levels = dim_stackedligandAV$`2nao_reslist`)

dim_stackedligandAV$rel_freq <- dim_stackedligandAV$ligand_freq / max(dim_stackedligandAV$ligand_freq)
dim_stackedreceptorAV$rel_freq <- dim_stackedreceptorAV$receptor_freq / max(dim_stackedreceptorAV$receptor_freq)

dim_stackedligandAV$color <- dim_stackedligandAV$`2nao_reslist` == "VAL2U" | dim_stackedligandAV$`2nao_reslist` == "VAL2V"

#save file for b-factor coloring

write.csv(dim_stackedligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_stackedligandAV.csv")

write.csv(dim_stackedreceptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_stackedreceptorAV.csv")

#WT
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "dimstack_WT/20201004_213717_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/dim_stackedWT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

dim_stackedreceptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dim_stackedreceptorWT) <- c("6IRA_reslist", "receptor_freq")
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
  dim_stackedreceptorWT[dim_stackedreceptorWT$`6IRA_reslist` == i, ]$receptor_freq <- dim_stackedreceptorWT[dim_stackedreceptorWT$`6IRA_reslist` == i, ]$receptor_freq + 1
}

dim_stackedligandWT$`2nao_reslist` <- factor(dim_stackedligandWT$`2nao_reslist`, levels = dim_stackedligandWT$`2nao_reslist`)

dim_stackedligandWT$rel_freq <- dim_stackedligandWT$ligand_freq / max(dim_stackedligandWT$ligand_freq)
dim_stackedreceptorWT$rel_freq <- dim_stackedreceptorWT$receptor_freq / max(dim_stackedreceptorWT$receptor_freq)

dim_stackedligandWT$color <- dim_stackedligandWT$`2nao_reslist` == "ALA2U" | dim_stackedligandWT$`2nao_reslist` == "ALA2V"

#save file for b-factor coloring

write.csv(dim_stackedligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_stackedligandWT.csv")

write.csv(dim_stackedreceptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_dim_stackedreceptorWT.csv")

##mono
#A2T
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "mono_AT/20201004_214022_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_A2T_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

mono_receptorAT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorAT) <- c("6IRA_reslist", "receptor_freq")
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
  mono_receptorAT[mono_receptorAT$`6IRA_reslist` == i, ]$receptor_freq <- mono_receptorAT[mono_receptorAT$`6IRA_reslist` == i, ]$receptor_freq + 1
}

mono_ligandAT$`2nao_reslist` <- factor(mono_ligandAT$`2nao_reslist`, levels = mono_ligandAT$`2nao_reslist`)

mono_ligandAT$rel_freq <- mono_ligandAT$ligand_freq / max(mono_ligandAT$ligand_freq)
mono_receptorAT$rel_freq <- mono_receptorAT$receptor_freq / max(mono_receptorAT$receptor_freq)

mono_ligandAT$color <- mono_ligandAT$`2nao_reslist` == "THR2U"

#save file for b-factor coloring

write.csv(mono_ligandAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_mono_ligandAT.csv")

write.csv(mono_receptorAT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_mono_receptorAT.csv")

#A2V
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "mono_AV/20201004_214324_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_A2V_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

mono_receptorAV <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorAV) <- c("6IRA_reslist", "receptor_freq")
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
  mono_receptorAV[mono_receptorAV$`6IRA_reslist` == i, ]$receptor_freq <- mono_receptorAV[mono_receptorAV$`6IRA_reslist` == i, ]$receptor_freq + 1
}

mono_ligandAV$`2nao_reslist` <- factor(mono_ligandAV$`2nao_reslist`, levels = mono_ligandAV$`2nao_reslist`)

mono_ligandAV$rel_freq <- mono_ligandAV$ligand_freq / max(mono_ligandAV$ligand_freq)
mono_receptorAV$rel_freq <- mono_receptorAV$receptor_freq / max(mono_receptorAV$receptor_freq)

mono_ligandAV$color <- mono_ligandAV$`2nao_reslist` == "VAL2U"

#save file for b-factor coloring

write.csv(mono_ligandAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_mono_ligandAV.csv")

write.csv(mono_receptorAV, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_mono_receptorAV.csv")

#WT
rm(list = ls())

dir <- "c:/users/wessec/documents/research/rpi gb/clusprooutput/6IRA_blocked/"

basename <- "mono_WT/20201004_214513_output.txt"
path <- paste(dir, basename, collapse = "", sep = "")

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

naoreslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/mono_WT_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(naoreslist) <- "2nao_reslist"
pereslist <- read.csv("c:/users/wessec/documents/research/rpi gb/misc/pdb extract/6IRA_reslist.txt", skip = 0, header = F, sep = "", stringsAsFactors = F)
colnames(pereslist) <- "6IRA_reslist"

mono_receptorWT <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(mono_receptorWT) <- c("6IRA_reslist", "receptor_freq")
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
  mono_receptorWT[mono_receptorWT$`6IRA_reslist` == i, ]$receptor_freq <- mono_receptorWT[mono_receptorWT$`6IRA_reslist` == i, ]$receptor_freq + 1
}

mono_ligandWT$`2nao_reslist` <- factor(mono_ligandWT$`2nao_reslist`, levels = mono_ligandWT$`2nao_reslist`)

mono_ligandWT$rel_freq <- mono_ligandWT$ligand_freq / max(mono_ligandWT$ligand_freq)
mono_receptorWT$rel_freq <- mono_receptorWT$receptor_freq / max(mono_receptorWT$receptor_freq)

mono_ligandWT$color <- mono_ligandWT$`2nao_reslist` == "ALA2U"

#save file for b-factor coloring

write.csv(mono_ligandWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_mono_ligandWT.csv")

write.csv(mono_receptorWT, "c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus_mono_receptorWT.csv")
