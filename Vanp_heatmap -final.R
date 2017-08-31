require(ggplot2)
require(reshape2)
require(RColorBrewer) # for customise colors
require(ggtree) # for ploting phylogenetic tree and other things
require(gplots) # for plotting heatmap with "heatmap.2" function 
require(ape)
require(gridGraphics)
require(gridExtra)



#set working directory to the folder containing current working R file
# install.packages("rstudioapi") # run this if it's your first time using it to install
require(rstudioapi) # load it
# the following line is for getting the path of your current open file
current_path <- getActiveDocumentContext()$path 
# The next line set the working directory to the relevant one:
setwd(dirname(current_path ))
# you can make sure you are in the right directory
print( getwd() )



#input the BI tree, which is in nexus format
t1 <- read.beast("./Data/Vanp.BI.nex.con.tre")
get.fields(t1)
# converting the BI tree into dataframe
a1 <- fortify(t1) %>% dplyr::as_data_frame()
# plotting the BI tree
t1_plot <- ggtree(t1, layout="rectangular", ladderize = T, right = F)+
  geom_tiplab(align = T, hjust=-0.1)+
  # geom_tippoint()+
  geom_treescale(x=0, y=0, offset = -1, linesize = 1)+
  geom_point2(aes(subset= !is.na(prob_percent) & prob_percent > 90 & isTip == FALSE), color='red', size =4, alpha=1)  # highlight nodes with bootstrap values greater than 
t1_plot



# because the tip label is "OUT" instead of "OTU", below is the code to fix this mistake. solution from https://guangchuangyu.github.io/ggtree/faq/#modify-tip-labels 
lb = get.tree(t1)$tip.label
lb
d = data.frame(label=lb)
d = data.frame(label=lb, label2 = chartr("OUT", "OTU", d$label))
# replotting the tree with the correct tip labels
t1_plot <- ggtree(t1, layout="rectangular", ladderize = T, right = F) %<+% d + geom_tiplab(aes(label=label2), align = T, hjust=-0.1)+
  geom_treescale(x=0, y=0, offset = -1, linesize = 1)+
  geom_point2(aes(subset= !is.na(prob_percent) & prob_percent > 90 & isTip == FALSE), color='red3', size =4, alpha=1)  # highlight nodes with bootstrap values greater than 
t1_plot




#load abundance data
d1 <- read.csv("./Data/Vanp_count_table_heatmap.csv", header = T, sep=",", row.names = 1, check.names = F)


#converting numeric to factor in the data
d1$`4.5-1` <- cut(d1$`4.5-1`, breaks = c(0, 0.001, 0.1, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),right = FALSE)
d1$`4.5-2` <- cut(d1$`4.5-2`, breaks = c(0, 0.001, 0.1, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),right = FALSE)
d1$`4.5-3` <- cut(d1$`4.5-3`, breaks = c(0, 0.001, 0.1, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),right = FALSE)
d1$`7.5-1` <- cut(d1$`7.5-1`, breaks = c(0, 0.001, 0.1, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),right = FALSE)
d1$`7.5-2` <- cut(d1$`7.5-2`, breaks = c(0, 0.001, 0.1, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),right = FALSE)
d1$`7.5-3` <- cut(d1$`7.5-3`, breaks = c(0, 0.001, 0.1, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),right = FALSE)



#gheatmap with default setting
gheatmap(t1_plot, d1)



#!!
#because the number of colours in the above colour scheme is smaller than the total number of catergories, we need to do the following to fix that.
colfunc<-colorRampPalette(c("yellow","purple","darkblue"))
plot(rep(1,13),col=(colfunc(13)), pch=19,cex=2)

cols=(colfunc(13))
names(cols) = c("[0,0.001)", "[0.001,0.1)", "[0.1,1)", "[1,10)", "[10,20)", "[20,30)", "[30,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,100)")
label_heatmap = c("0-0.001", "0.001-0.1", "0.1-1", "1-10", "10-20", "20-30", "50-60", "60-70",  "90-100")
legend_title <- "relative abundance (%)"

gheatmap(t1_plot, d1, color="grey70", width = 1, offset = 0.05) + scale_fill_manual(values = cols)
gheatmap(t1_plot, d1, color="grey70", width = 1, offset = 0.05) + scale_fill_manual(legend_title, values = cols, labels=label_heatmap)



#to change the position of the colour key
gheatmap(t1_plot, d1, color="grey70", width = 0.5, offset = 0.01, font.size = 4) + 
  scale_fill_manual(values = cols, labels=label_heatmap) +
  theme(legend.position="left")

d1 <- read.csv("ddd.csv", header = T, sep=",", row.names = 1, check.names = F)


#to change the position of the colour key more procisely
gheatmap(t1_plot, d1, color="white", width = 0.8, offset = 0.03, font.size = 6) + 
  scale_fill_manual(legend_title, values = cols, labels=label_heatmap) +
  theme(legend.position=c(0.1,0.85))+
  geom_text2(aes(subset= !is.na(prob_percent) & isTip == FALSE, label=prob_percent), nudge_x = -0.01 , nudge_y = 0.4, size=3, colours="red3")