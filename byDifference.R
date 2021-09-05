#Uses data formatted and calculated in Organizer scripts
#reformats for making by difference graphs between WT and mutants
#### Package Setup #####

rm(list=ls())

#required packages
library(ggplot2)
library(plyr)
library(tidyverse)

##### Zdock data #####
###lateral dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 200) - (dataWT$ligand_freq / 200)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 200) - (dataWT$ligand_freq / 200)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Lateral Dimer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_fill_manual(values = (color_tibble$color))


#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Lateral Dimer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_fill_manual(values = (color_tibble$color))


###stacked dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 200) - (dataWT$ligand_freq / 200)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 200) - (dataWT$ligand_freq / 200)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Stacked Dimer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_fill_manual(values = (color_tibble$color))


#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Stacked Dimer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_fill_manual(values = (color_tibble$color))


###monomer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/mono_ligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 200) - (dataWT$ligand_freq / 200)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 200) - (dataWT$ligand_freq / 200)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Monomer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_fill_manual(values = (color_tibble$color))


#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Monomer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_fill_manual(values = (color_tibble$color))

##### Zdock blocked 4pe5 #####
###lateral dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 172) - (dataWT$ligand_freq / 172)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 170) - (dataWT$ligand_freq / 172)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Lateral Dimer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.01, ymax = 0.05, alpha = 0.7, fill = "coral2") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.01, ymax = 0.03, alpha = 0.7, fill = "coral2") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Lateral Dimer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.01, ymax = 0.10, alpha = 0.7, fill = "coral2") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.01, ymax = 0.06, alpha = 0.7, fill = "coral2") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


###stacked dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 196) - (dataWT$ligand_freq / 196)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 196) - (dataWT$ligand_freq / 196)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Stacked Dimer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.02, ymax = 0.01, alpha = 0.7, fill = "coral2") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.02, ymax = 0.01, alpha = 0.7, fill = "coral2") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)


#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Stacked Dimer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.01, ymax = 0.12, alpha = 0.7, fill = "coral2") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.01, ymax = 0.14, alpha = 0.7, fill = "coral2") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.062, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


###monomer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/mono_ligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 103) - (dataWT$ligand_freq / 101)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 106) - (dataWT$ligand_freq / 101)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Monomer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.01, ymax = 0.065, alpha = 0.7, fill = "coral2") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Monomer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.01, ymax = 0.12, alpha = 0.7, fill = "coral2") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

##### Zdock blocked 4pe5 2#####
#redo images in new format
###lateral dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 172) - (dataWT$ligand_freq / 172)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 170) - (dataWT$ligand_freq / 172)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif

dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Lateral Dimer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2X", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Lateral Dimer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


###stacked dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 196) - (dataWT$ligand_freq / 196)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 196) - (dataWT$ligand_freq / 196)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2V", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Stacked Dimer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2V", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Stacked Dimer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.062, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


###monomer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/mono_ligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 103) - (dataWT$ligand_freq / 101)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 106) - (dataWT$ligand_freq / 101)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Monomer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Monomer Ligand by Freq. Difference (Zdock,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


##### Zdock blocked 6IRA #####
###lateral dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 189) - (dataWT$ligand_freq / 190)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 188) - (dataWT$ligand_freq / 190)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif

dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Lateral Dimer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2X", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Lateral Dimer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

#extra

ggplot(mapping = aes(label = as.character(dataAV$`2nao_reslist`))) +
  geom_point(aes(x = dataAT$bydif, y = dataAV$bydif)) +
  geom_text(aes(x = dataAT$bydif, y = dataAV$bydif), hjust=-.5,vjust=.5) +
  geom_abline(slope = 1, intercept = 0)

###stacked dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 177) - (dataWT$ligand_freq / 176)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 172) - (dataWT$ligand_freq / 176)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2V", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Stacked Dimer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.065, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2V", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Stacked Dimer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.062, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


###monomer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 143) - (dataWT$ligand_freq / 144)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 145) - (dataWT$ligand_freq / 144)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Monomer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Monomer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

##### Zdock blocked 6IRA_2#####
#alternate coloring method
###lateral dimer
rm(list = ls())

#load in data

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 189) - (dataWT$ligand_freq / 190)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAT[c(1:16, 43:58), ]$color <- T

dataAV$bydif <- (dataAV$ligand_freq / 188) - (dataWT$ligand_freq / 190)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAV[c(1:16, 43:58), ]$color <- T

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif

dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  labs(title = "Lateral Dimer A2T - WT", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  scale_fill_manual(values = c("steelblue4", "tomato"))

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2X", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif), show.legend = F) +
  labs(title = "Lateral Dimer A2V - WT", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  scale_fill_manual(values = c("steelblue4", "tomato"))

summary(dataAV$bydif)

#extra

ggplot(mapping = aes(label = as.character(dataAV$`2nao_reslist`))) +
  geom_point(aes(x = dataAT$bydif, y = dataAV$bydif)) +
  geom_text(aes(x = dataAT$bydif, y = dataAV$bydif), hjust=-.5,vjust=.5) +
  geom_abline(slope = 1, intercept = 0)

###stacked dimer
rm(list = ls())

#load in data

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 177) - (dataWT$ligand_freq / 176)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAT[c(1:16, 43:58), ]$color <- T

dataAV$bydif <- (dataAV$ligand_freq / 172) - (dataWT$ligand_freq / 176)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAV[c(1:16, 43:58), ]$color <- T

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2V", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif), show.legend = F) +
  labs(title = "Stacked Dimer A2T - WT", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.065, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  scale_fill_manual(values = c("steelblue4", "tomato"))

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2V", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif), show.legend = F) +
  labs(title = "Stacked Dimer A2V - WT", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.062, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  scale_fill_manual(values = c("steelblue4", "tomato"))

summary(dataAV$bydif)


###monomer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_ligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 143) - (dataWT$ligand_freq / 144)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 145) - (dataWT$ligand_freq / 144)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Monomer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Monomer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)



##### Zdock blocked 3SQ9 #####
###lateral dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 145) - (dataWT$ligand_freq / 146)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 160) - (dataWT$ligand_freq / 146)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif

dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Lateral Dimer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2X", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Lateral Dimer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

#extra

ggplot(mapping = aes(label = as.character(dataAV$`2nao_reslist`))) +
  geom_point(aes(x = dataAT$bydif, y = dataAV$bydif)) +
  geom_text(aes(x = dataAT$bydif, y = dataAV$bydif), hjust=-.5,vjust=.5) +
  geom_abline(slope = 1, intercept = 0)

###stacked dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 112) - (dataWT$ligand_freq / 114)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 111) - (dataWT$ligand_freq / 114)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2V", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Stacked Dimer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.065, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2V", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Stacked Dimer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.062, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


###monomer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3SQ9_blocked/zdock/mono_ligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 36) - (dataWT$ligand_freq / 36)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 36) - (dataWT$ligand_freq / 36)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)


#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Monomer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Monomer Ligand by Freq. Difference (Zdock,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

##### Zdock blocked 3SQ9_2#####
#alternate coloring method
###lateral dimer
rm(list = ls())

#load in data

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3sq9_blocked/zdock/dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3sq9_blocked/zdock/dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/3sq9_blocked/zdock/dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 145) - (dataWT$ligand_freq / 146)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAT[c(1:16, 43:58), ]$color <- T

dataAV$bydif <- (dataAV$ligand_freq / 160) - (dataWT$ligand_freq / 146)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAV[c(1:16, 43:58), ]$color <- T

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif

dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  labs(title = "Lateral Dimer A2T - WT", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  scale_fill_manual(values = c("steelblue4", "tomato"))

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2X", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif), show.legend = F) +
  labs(title = "Lateral Dimer A2V - WT", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.08, 0.14), breaks = seq(-.08, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  scale_fill_manual(values = c("steelblue4", "tomato"))

summary(dataAV$bydif)

#extra

ggplot(mapping = aes(label = as.character(dataAV$`2nao_reslist`))) +
  geom_point(aes(x = dataAT$bydif, y = dataAV$bydif)) +
  geom_text(aes(x = dataAT$bydif, y = dataAV$bydif), hjust=-.5,vjust=.5) +
  geom_abline(slope = 1, intercept = 0)

###stacked dimer
rm(list = ls())

#load in data

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 112) - (dataWT$ligand_freq / 114)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAT[c(1:16, 43:58), ]$color <- T

dataAV$bydif <- (dataAV$ligand_freq / 111) - (dataWT$ligand_freq / 114)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAV[c(1:16, 43:58), ]$color <- T

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2V", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif), show.legend = F) +
  labs(title = "Stacked Dimer A2T - WT", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.065, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  scale_fill_manual(values = c("steelblue4", "tomato"))

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2V", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif), show.legend = F) +
  labs(title = "Stacked Dimer A2V - WT", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.062, 0.14), breaks = seq(-.06, 0.14, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  scale_fill_manual(values = c("steelblue4", "tomato"))

summary(dataAV$bydif)



##### ClusPro data #####
###lateral dimer
rm(list = ls())

#load in data
dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_unblocked/clus_dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_unblocked/clus_dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_unblocked/clus_dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 25) - (dataWT$ligand_freq / 24)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 23) - (dataWT$ligand_freq / 24)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  labs(title = "A2T Lateral Dimer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(-0.25, 0.25), breaks = seq(-.25, 0.25, 0.05))

#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  labs(title = "A2V Lateral Dimer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(-0.25, 0.25), breaks = seq(-.25, 0.25, 0.05))

###stacked dimer
rm(list = ls())

#load in data
dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 30) - (dataWT$ligand_freq / 29)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 29) - (dataWT$ligand_freq / 29)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  labs(title = "A2T Stacked Dimer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(-0.25, 0.25), breaks = seq(-.25, 0.25, 0.05))

#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  labs(title = "A2V Stacked Dimer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(-0.25, 0.25), breaks = seq(-.25, 0.25, 0.05))

###monomer
rm(list = ls())

#load in data
dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_mono_ligandWT.csv", stringsAsFactors = F, header = T)

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 30) - (dataWT$ligand_freq / 30)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 30) - (dataWT$ligand_freq / 30)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

#plot AT
ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  labs(title = "A2T Monomer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(-0.25, 0.25), breaks = seq(-.25, 0.25, 0.05))

#plot AV
ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = color), show.legend = F) +
  labs(title = "A2V Monomer Ligand Residue Frequency by Difference", x = "Residue", y = "Difference in Frequency") +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(-0.25, 0.25), breaks = seq(-.25, 0.25, 0.05))


##### blocked ClusPro data #####
###lateral dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                               "TYR", "CYS", "ALA",
                               "THR", "HIS", "GLY", "SER", "GLN",
                               "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                               rep("slight hydrophobic", 3),
                               rep("neutral", 5),
                               rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 29) - (dataWT$ligand_freq / 27)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 27) - (dataWT$ligand_freq / 27)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Lateral Dimer Ligand by Freq. Difference (ClusPro,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.22), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2X", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Lateral Dimer Ligand by Freq. Difference (ClusPro,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.20), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

###stacked dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 29) - (dataWT$ligand_freq / 29)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 28) - (dataWT$ligand_freq / 29)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2V", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Stacked Dimer Ligand by Freq. Difference (ClusPro,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.22), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2V", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Stacked Dimer Ligand by Freq. Difference (ClusPro,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.22), breaks = seq(-.26, 0.22, 0.02)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

###monomer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/clus_blocked/clus_mono_ligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 22) - (dataWT$ligand_freq / 20)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 22) - (dataWT$ligand_freq / 20)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Monomer Ligand by Freq. Difference (ClusPro,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.22), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Monomer Ligand by Freq. Difference (ClusPro,4PE5)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.22), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

##### blocked ClusPro data 6IRA #####
###lateral dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_dim_lateralligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_dim_lateralligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_dim_lateralligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 25) - (dataWT$ligand_freq / 24)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 23) - (dataWT$ligand_freq / 24)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2X", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Lateral Dimer Ligand by Freq. Difference (ClusPro,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.22), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2X", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Lateral Dimer Ligand by Freq. Difference (ClusPro,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1X", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1X", xmax = "GLU3X", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.20), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1X", "ALA42X")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

###stacked dimer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_dim_stackedligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 30) - (dataWT$ligand_freq / 29)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 29) - (dataWT$ligand_freq / 29)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)

dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif
b <- dataAT[dataAT$X2nao_reslist == "THR2V", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Stacked Dimer Ligand by Freq. Difference (ClusPro,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.34, 0.22), breaks = seq(-.34, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAT$bydif)

#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif
b <- dataAV[dataAV$X2nao_reslist == "VAL2V", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Stacked Dimer Ligand by Freq. Difference (ClusPro,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "ASP1V", ymin = -Inf, ymax = Inf, alpha=0.1, fill = "cyan") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  annotate("rect", xmin = "ASP1V", xmax = "GLU3V", ymin = -0.005 + (0.01 * is.na(sqrt(b))), ymax = b + (0.005 - 0.01 * is.na(sqrt(b))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.34, 0.22), breaks = seq(-.34, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ASP1V", "ALA42V")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)

#extra

ggplot(mapping = aes(label = as.character(dataAV$`2nao_reslist`))) +
  geom_point(aes(x = dataAT$bydif, y = dataAV$bydif)) +
  geom_text(aes(x = dataAT$bydif, y = dataAV$bydif), hjust=-.5,vjust=.5)
  
  geom_label(data = dataAT$`2nao_reslist`)


geom_label(aes(mapping = dataAT$`2nao_reslist`))

ggplot(data=df,aes(x=A,y=B,label=genes)) +
  geom_point(aes(color=group)) +
  geom_text(hjust=-1,vjust=1)

###monomer
rm(list = ls())

#load in data
#hydrophobicity data from https://www.peptide2.com/N_peptide_hydrophobicity_hydrophilicity.php at pH 7
hydrophobicity <- data.frame(residue = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET",
                                         "TYR", "CYS", "ALA",
                                         "THR", "HIS", "GLY", "SER", "GLN",
                                         "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"),
                             hydrophobicity = c(rep("hydrophobic", 6),
                                                rep("slight hydrophobic", 3),
                                                rep("neutral", 5),
                                                rep("hydrophilic", 6)))

color_tibble <- tibble(hydrophobicity = c("hydrophobic", "slight hydrophobic", "neutral", "hydrophilic"),
                       color = c("red", "coral3", "chartreuse3", "deepskyblue4"))

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_mono_ligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_mono_ligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6IRA_blocked/clus/clus_mono_ligandWT.csv", stringsAsFactors = F, header = T)

#assign correct hydrophobicity factor
dataAT$residue <- substr(dataAT$X2nao_reslist, 1, 3)
dataAT <- join(dataAT, hydrophobicity, by = "residue")

dataAV$residue <- substr(dataAV$X2nao_reslist, 1, 3)
dataAV <- join(dataAV, hydrophobicity, by = "residue")

dataWT$residue <- substr(dataWT$X2nao_reslist, 1, 3)
dataWT <- join(dataWT, hydrophobicity, by = "residue")

#manipulate data
dataAT$bydif <- (dataAT$ligand_freq / 30) - (dataWT$ligand_freq / 30)
dataAT$'2nao_reslist' <- factor(dataAT$X2nao_reslist, levels = dataAT$X2nao_reslist)

dataAV$bydif <- (dataAV$ligand_freq / 30) - (dataWT$ligand_freq / 30)
dataAV$'2nao_reslist' <- factor(dataAV$X2nao_reslist, levels = dataAV$X2nao_reslist)


dataAT$hydrophobicity <- factor(dataAT$hydrophobicity, levels = color_tibble$hydrophobicity)
dataAV$hydrophobicity <- factor(dataAV$hydrophobicity, levels = color_tibble$hydrophobicity)

#plot AT
a <- dataAT[dataAT$X2nao_reslist == "THR2U", ]$bydif

ggplot(data = dataAT) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2T Monomer Ligand by Freq. Difference (ClusPro,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.22), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


#plot AV
a <- dataAV[dataAV$X2nao_reslist == "VAL2U", ]$bydif

ggplot(data = dataAV) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T) +
  labs(title = "A2V Monomer Ligand by Freq. Difference (ClusPro,6IRA)", x = "Residue", y = "Difference in Frequency") +
  annotate("rect", xmin = "ASP1U", xmax = "GLU3U", ymin = -0.005 + (0.01 * is.na(sqrt(a))), ymax = a + (0.005 - 0.01 * is.na(sqrt(a))), alpha = 0.7, fill = "darkslategrey") +
  theme(panel.background = element_rect(fill = "grey"),
        panel.grid.major.x = element_line(colour = "white", size = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20)) +
  scale_y_continuous(limits = c(-0.26, 0.22), breaks = seq(-.26, 0.22, 0.04)) +
  scale_x_discrete(breaks = c("ASP1U", "ALA42U")) +
  scale_fill_manual(values = (color_tibble$color)) +
  geom_col(aes(x = `2nao_reslist`, y = bydif, fill = hydrophobicity), show.legend = T)

summary(dataAV$bydif)


##### Zdock interface plots 6IRA #####
#modified Sweta program
library(dplyr)
library(ggplot2)
library(reshape)
library(hash)
#Clean up everything 
rm(list= ls())

#Set ggplot theme settings
t=  theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 10, face = "bold"),
          #   legend.title = element_blank(),
          plot.title=element_text(hjust=0.5,face="bold", size=14),
          axis.title=element_text(size=10, face="bold"),
          axis.text=element_text(size=8,face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))



setwd('c:/users/wessec/documents/research/rpi gb/data/interface')

#Colors for different chains in the receptor
colch=hash()
keys=c("Al", "Am", "Ah", "Bl", "Bm", "Bh","Cl", "Cm", "Ch","Dl", "Dm", "Dh")
keyval=c("cornflowerblue" , "cornsilk","red", "cornflowerblue", "honeydew", "red", "cornflowerblue", "lavenderblush","red", "darkslateblue", "lightcyan","red")
for (i in seq(1,length(keys))){
  colch[[keys[i]]]=keyval[i]
}

#mono, lateral dimer and stacked dimer
cases=c("mono_A2T","mono_A2V",
        "dim_stackedA2V","dim_stackedA2T",
        "dim_lateralA2V","dim_lateralA2T")
#cases=c("mono_WT")
chains=c("A","B","C","D")
#chains=c("D")
for (c in cases){
  datall=read.table(paste(c,'/res_contact_freq.dat',sep=''))
  datall[which(datall$V5 >=1.0),] #artifact
  
  wtname <- gsub("A2T|A2V", "WT", c) #replace current case with WT form
  datawtall <- read.table(paste(wtname, '/res_contact_freq.dat', sep = '')) #load in WT data
  
  for (ch in chains){
    data=datall[which(datall$V3==ch),]
    data$V1 = factor(data$V1, levels=unique(data$V1))
    data$V2 = factor(data$V2, levels=unique(data$V2))
    
    datawt <- datawtall[which(datawtall$V3 == ch), ] #only use specific chain
    
    bydif <- data$V5 - datawt$V5 #caculate by differnce values
    
    p=ggplot(data, aes(V2,V1, fill= bydif)) + 
      geom_tile() + 
      #coord_equal() + 
      coord_fixed(ratio=0.25) + 
      scale_y_discrete(breaks = levels(data$V1)[c(T, rep(F, 39))]) + 
      scale_x_discrete(breaks = levels(data$V2)[c(T, rep(F, 19))]) + 
      scale_fill_gradient2(low=colch[[paste(ch,"l",sep='')]], mid = colch[[paste(ch, "m", sep = '')]], high=colch[[paste(ch,"h",sep='')]], midpoint = 0, limits = c(min(bydif), max(bydif))) + 
      xlab("Ligand\n (L:Chain:ResidueID)") + 
      ylab("Receptor\n (R:Chain:ResidueID)") + 
      t  
    #theme(legend.position="bottom",
    #legend.box="horizontal") 
    p 
    ggsave(paste('bydif/',c,'_',ch,'_contactmap_bydif.png',sep=''), plot=p, width = 6, height = 8,  units= "in" ) 
  }
}

summary(data$V5)
summary(data$V3)
summary(datall$V3)

typeof(datawtall$V5)


##### misc #####
rm(list = ls())

dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralreceptorAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralreceptorAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/zdock_blocked/dim_lateralreceptorWT.csv", stringsAsFactors = F, header = T)

summary(dataAT$receptor_freq)
summary(dataAV$receptor_freq)
summary(dataWT$receptor_freq)

#which(dataWT[dataWT$receptor_freq == 39, ]$receptor_freq)

dataAT[which(dataAT$receptor_freq > 36.00),]
dataAV[which(dataAV$receptor_freq > 36.00),]
dataWT[which(dataWT$receptor_freq > 36.00),]



dataAT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/clus/clus_dim_stackedligandAT.csv", stringsAsFactors = F, header = T)
dataAV <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/clus/clus_dim_stackedligandAV.csv", stringsAsFactors = F, header = T)
dataWT <- read.csv("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/clus/clus_dim_stackedligandWT.csv", stringsAsFactors = F, header = T)


ggplot() +
  geom_point(aes(x = dataAT$ligand_freq, y = dataAV$ligand_freq))
