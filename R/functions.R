require(ggplot2)
require(reshape2)
require(RColorBrewer) # for customise colors
require(ggtree) # for ploting phylogenetic tree and other things
require(gplots) # for plotting heatmap with "heatmap.2" function 
require(ape)
require(gridGraphics)
require(gridExtra)


# install.packages("rstudioapi") # run this if it's your first time using it to install
require(rstudioapi) # load it
# the following line is for getting the path of your current open file
current_path <- getActiveDocumentContext()$path 
# The next line set the working directory to the relevant one:
setwd(dirname(current_path ))
# you can make sure you are in the right directory
print( getwd() )