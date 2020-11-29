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
keys=c("Al", "Ah", "Bl", "Bh","Cl","Ch","Dl", "Dh")
keyval=c("cornsilk","red","honeydew","red","lavenderblush","red","lightcyan","red")
for (i in seq(1,length(keys))){
  colch[[keys[i]]]=keyval[i]
}

#mono, lateral dimer and stacked dimer
cases=c("mono_WT","mono_A2T","mono_A2V",
        "dim_stackedWT","dim_stackedA2V","dim_stackedA2T",
        "dim_lateralWT","dim_lateralA2V","dim_lateralA2T")
#cases=c("mono_WT")
chains=c("A","B","C","D")
#chains=c("D")
for (c in cases){
  datall=read.table(paste(c,'/res_contact_freq.dat',sep=''))
  datall[which(datall$V5 >=1.0),]
  for (ch in chains){
    data=datall[which(datall$V3==ch),]
    data$V1 = factor(data$V1, levels=unique(data$V1))
    data$V2 = factor(data$V2, levels=unique(data$V2))
    p=ggplot(data, aes(V2,V1, fill= V5)) + 
     geom_tile() + 
      #coord_equal() + 
      coord_fixed(ratio=0.25) + 
      scale_y_discrete(breaks = levels(data$V1)[c(T, rep(F, 39))]) + 
      scale_x_discrete(breaks = levels(data$V2)[c(T, rep(F, 19))]) + 
      scale_fill_gradient(low=colch[[paste(ch,"l",sep='')]], high=colch[[paste(ch,"h",sep='')]], limits=c(0, max(datall$V5))) + 
      xlab("Ligand\n (L:Chain:ResidueID)") + 
      ylab("Receptor\n (R:Chain:ResidueID)") + 
      t  
      #theme(legend.position="bottom",
    #legend.box="horizontal") 
   p 
    ggsave(paste('Plots/',c,'_',ch,'_contactmap.png',sep=''), plot=p, width = 6, height = 8,  units= "in" ) 
  }
}


