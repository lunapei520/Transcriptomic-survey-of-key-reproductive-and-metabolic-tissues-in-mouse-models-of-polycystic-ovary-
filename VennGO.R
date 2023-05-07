library(devtools)
library(ggvenn)
library(tidyr)
library(dplyr)
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn") # install via GitHub (for latest version)
df1 = read.table("Hypo17NFgenes.txt", header= F, sep = "\t")
df2 = read.table("HypoHVgenes.txt", header = F, sep = "\t")
df3 = read.table("HypoPNAgenes.txt", header = F, sep = "\t")
df4 = read.table("HypoPPgenes.txt", header = F, sep = "\t")
x = list(Prepubertal = df4$V1, PNA = df3$V1, pHFHS = df2$V1, NF = df1$V1)

ggvenn(x, c("PNA", "pHFHS","Prepubertal", "NF"), show_percentage = FALSE, fill_color = c("orange" , "purple4","seagreen", "steelblue4"),
       fill_alpha = 0.5,stroke_size = 0.01,text_size = 6, set_name_size=6)

mydat = read.table("GOhypopvalue.txt", header = T, sep = "\t")

mydat$Type <- factor(mydat$Type, levels = c('PNA','HFHS','Prepubertal','NF')) 
 mydat$Class <- factor(mydat$Class, levels = c('PNA','Prepubertal','17NF','HV'))  
 ggplot(mydat, aes(x=LogP, reorder(GO.term, as.numeric(Type)), size= Number)) + geom_point(aes(color = Type),alpha = 1) + 
 scale_color_viridis(discrete = TRUE, option = "D")+
 geom_tile(aes(width = Inf, fill = Class),alpha = 0.2) +  
 scale_fill_manual(values = c("purple4","steelblue4", "seagreen","orange"  ))+
 theme_classic() +
 theme(axis.text = element_text(size = 15))+
 ggtitle("Hypothalamus GO")