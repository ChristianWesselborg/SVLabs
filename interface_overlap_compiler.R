library(ggplot2)
#Clean up everything 
rm(list= ls())

#Set ggplot theme settings
t = theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          #legend.text = element_text(size = 15, face = "bold"),
          plot.title=element_text(hjust=0.5,face="bold", size=14),
          axis.title=element_text(size=15, face="bold"),
          axis.text=element_text(size=8,face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

path = 'c:/users/wessec/documents/research/rpi gb/data/interface overlap'
#path = 'c:/users/wessec/documents/research/rpi gb/test2'

dirpath = dir(path)


for(folder in dirpath){

  files = list.files(file.path(path, folder))

  for(file in files){

    data = read.csv(file = file.path(path, folder, file), skip = 2)
    data = data[data[ , 2] > 0, ]

    data[ , 1] = as.factor(data[ , 1])

    p = ggplot(data = data) +
      geom_col(aes(x = data[ , 1], y = data[ , 2]), orientation = "x", fill = "darkcyan") +
      t +
      labs(title = paste(folder, gsub(".txt", "", file), "Interface Overlap Plot", sep = "_"), x = "Model Number", y = "Residue Overlap")

    p
    ggsave(paste(folder, gsub(".txt", "", file), ".png"), path = "c:/users/wessec/documents/research/rpi gb/misc/graphs/interface overlap", plot=p, width = 6, height = 8,  units= "in" )
  }
}


#copy for testing purposes
# 
#   for(file in dirpath){
#     
#     data = read.csv(file = file.path(path, file), skip = 2)
#     data = data[data[ , 2] > 0, ]
#     
#     data[ , 1] = as.factor(data[ , 1])
#     
#     p = ggplot(data = data) +
#       geom_col(aes(x = data[ , 1], y = data[ , 2]), orientation = "x", fill = "darkcyan") +
#       t +
#       labs(title = paste(path, gsub(".txt", "", file), "Interface Overlap Plot", sep = "_"), x = "Model Number", y = "Residue Overlap")
#     
#     p
#     #ggsave(paste(folder, gsub(".txt", "", file), ".png"), path = "c:/users/wessec/documents/research/rpi gb/misc/graphs/interface overlap", plot=p, width = 6, height = 8,  units= "in" ) 
#     ggsave(paste(gsub(".txt", "", file), ".png"), path = "c:/users/wessec/documents/research/rpi gb/test2", plot=p, width = 6, height = 8,  units= "in" )
#   }
