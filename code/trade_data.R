######################################
#' this script provides the analysis for the trade dataset in the main paper.  
#' This script also contains the code to produce Figure 1 and figure 7.
######################################

#Load required packages:
source("misc.R")
source("HOOI.R")
source("membership_estimation.R")
source("estimate_rank.R")
require(maps)
library(viridis)
library(RColorBrewer)
library(maptools)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(igraph)

###########################################
# Process data
############################################'

nodes <- read.delim("../data/fao_trade_nodes.txt",sep = "")
layers <- read.delim("../data/fao_trade_layers.txt",sep="")
edges <- read.table("../data/fao_trade_multiplex.edges", sep="")

# we need to process our data to be a tensor object.  we will first
# store a list of adjacency matrices that we will then combine into a tensor.
all_Graphs <- list()
edges$weight <- edges$V4

# we now collect adjacency matrices
for ( i in c(1:nrow(layers))) {
  # first we create a dataframe of weighted edges between nodes
  layer_dat <- edges[which(edges$V1 == i),c("V2","V3","weight")]
  
  # we then create the corresponding graph, convert to undirected
  G <- graph_from_data_frame(layer_dat,vertices= nodes)
  G <-  as.undirected(G,mode="collapse")
  all_Graphs[[i]] <- G
}

nLayers <- nrow(layers)
nLayers
length(V(G))

# we will now process the data
toKeep <- list() # the list of networks we keep
layers_kept <- c() # the list of indices of networks we keep
j <- 1 

for (i in c(1:nLayers)) {
  
  # first we find the largest connected component in each network
  comps <- igraph::clusters(all_Graphs[[i]], mode="weak")
  
  # if the largest connected component is smaller than 150, we do not keep it
  if (max(comps$csize) >= 150) {
    toKeep[[j]] <- all_Graphs[[i]]
    layers_kept <- c(layers_kept,i)
    j <- j+1
  }
}

# we now update the number of layers
nLayers <- length(toKeep)
nLayers

# we now re-determine the dimensions of the tensor after processing
nV <- length(V(toKeep[[1]]))
nLayers <- length(toKeep)

# we now instantiate the tensor
tens <- array(0,dim=c(nLayers,nV,nV))

# within each layer we will store the weighted adjacency matrix
for ( i in c(1:nLayers)) {
  A <- as_adjacency_matrix(toKeep[[i]],attr="weight")
  rownames(A) <- V(toKeep[[i]])$nodeLabel
  colnames(A) <- V(toKeep[[i]])$nodeLabel
  tens[i,,] <- as.matrix(A)
}

#finally, we obtain the tensor object
tens <- as.tensor(tens)
dim(tens)

###########################################
#start analysis
###########################################

#first try to embed by deleting the diagonal:
mat1 <- k_unfold(tens,1)@data
toembed <- mat1 %*% t(mat1)
diag(toembed) <- 0
embed1 <- eigen(toembed)
togetelbows <- embed1$values[embed1$values >= 0]
togetelbows #only a few positive values
r1 <- getElbows(sqrt(togetelbows) )
r1
r1 <- r1[2] #results in only 8 positive eigenvalues.  But we use r1 = 5

#get second embedding dimension
mat2 <- k_unfold(tens,2)@data
toembed <- mat2 %*% t(mat2)
diag(toembed) <- 0
embed1 <- eigen(toembed)
togetelbows <- embed1$values[embed1$values >= 0]
togetelbows <- sqrt(togetelbows)#only look at postitive
r2 <- getElbows(togetelbows)
r2 #4 11 23, so choose 4 for interpretability
r2 <- r2[1]

# #get third embedding dimension
 mat3 <- k_unfold(tens,3)@data
 toembed <- mat3 %*% t(mat3)
 diag(toembed) <- 0
 embed1 <- eigen(toembed)
 togetelbows <- embed1$values[embed1$values >= 0]
 togetelbows <- sqrt(togetelbows)#only look at postitive
 r3 <- getElbows(togetelbows)
 r3 #4 11 23, so choose 4 for interpretability
 r3 <- r3[1]


rs <- c(r1,r2,r3)

#find embedding and estimate memberships
Uhats <- HOOI_dd(tens,r= rs, niter=20)
clusts <- SPAMM(uhats = Uhats,threshold=T,tval=10^(-20)) #outputs mixed memberships

###################
# process output
########################

#first determine the labels of the goods kept:
#(note clusts[[2]] is a list of three different indices of labels)
layers[layers_kept[clusts[[2]][[1]]],"layerLabel"]

# store the mixed memberships as a data frame
dat_layers <- round(clusts[[1]][[1]],3)
dat_layers <- as.data.frame(dat_layers)
names(dat_layers) <- layers[layers_kept[clusts[[2]][[1]]],"layerLabel"]
dat_layers$good <- layers[layers_kept,"layerLabel"]
View(dat_layers)

# now we combine processed food into a single component
dat_layers$processed <- 
  dat_layers$Food_prep_nes + 
  dat_layers$`Cheese,_whole_cow_milk` + dat_layers$`Beverages,_distilled_alcoholic`

#same with unprocessed:
dat_layers$unprocessed <- dat_layers$Crude_materials + dat_layers$Maize

View(dat_layers[,c("processed","unprocessed","Beverages,_distilled_alcoholic","good")])

#below outputs primarily processed, primarily unprocessed, and primarily neither
dat_layers[dat_layers$processed > .7,"good"]
dat_layers[dat_layers$unprocessed > .7,"good"]
dat_layers[dat_layers$unprocessed <= .7 &dat_layers$processed <= .7 ,"good"]

View(dat_layers[dat_layers$processed <= .7,c("unprocessed","Crude_materials","Maize","good")])

# we now analyze the country mode: first we figure out the labels and edit for
# readability
countries <- nodes[V(toKeep[[1]]),"nodeLabel"]
order(as.character(nodes[V(toKeep[[1]]),"nodeLabel"]))
countries <- gsub("_", " ", countries)
countries[countries == "Viet Nam"] <- "Vietnam"
countries[countries == "United States of America"] <- "USA"
countries[countries == "Venezuela (Bolivarian Republic of)"] <- "Venezuela"
countries[countries == "Bolivia (Plurinational State of)"] <- "Bolivia"
countries[countries == "Republic of Korea"] <- "South Korea"
countries[countries == "Russian Federation"] <- "Russia"
countries[countries == "Iran (Islamic Republic of)"] <- "Iran"
countries[countries == "The former Yugoslav Republic of Macedonia"] <- "Macedonia"
countries[countries == "United Kingdom"] <- "UK"
countries[countries == "Republic of Moldova"] <- "Moldova"
countries[countries=="Micronesia (Federated States of)" ] <- "Micronesia"
countries[countries=="Democratic People's Republic of Korea"  ] <- "North Korea"
countries[countries=="Netherlands Antilles"] <- "Netherlands"
countries[countries == "Syrian Arab Republic" ] <- "Syria"
countries

# we determine the pure nodes for each mode
pure_nodes1 <-countries[clusts[[2]][[2]]]
pure_nodes2 <- countries[clusts[[2]][[3]]]  
#they are the same as we only used undirected relationships

##########################
# Create Figure 1
##########################


world <- map_data("world")
world <- fortify(world, region="id")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(0, 1))

#function below takes in a pure node and plots the intensity of membership 
# for that pure node for each country based on the palette above.
# note: if multiple=TRUE then we can combine pure nodes.  this is needed for Figure
# 7 in the paper, where we plot both USA and Canada together.
plot_pure_node <- function(node,multiple=FALSE,legend=FALSE) {
  if(multiple) {
    tit <- "Pure Node: "
    ddf$value <- 0
    for ( i in c(1:length(node))) {
      eval(parse(text=paste0("value <- ddf$",node[i])))
      ddf$value <- value + ddf$value
      if (i > 1) {
        tit <- paste(tit, "and")
      }
      tit <- paste(tit,node[i])
    }
  } else {
    eval(parse(text=paste0("ddf$value <- ddf$",node)))
    tit <- paste0("Pure Node: ",node)
  }
  
  
  world %>%
    merge(ddf, by.x = "region", by.y = "country", all.x = T) %>%
    arrange(group, order) %>%
    ggplot(aes(x = long, y = lat, group = group, fill = value)) +
    geom_polygon(color = "white", size = 0.2) +
    sc +
    scale_y_continuous(limits=c(-60,90))+
      theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank()) +
    ggtitle(tit)  

}

# we now plot the different pure node memberships.  First we format the data:
dat <- round(as.data.frame(clusts[[1]][[2]]),4)
var <-pure_nodes1
var
names(dat) <- var #label the communities by pure node
dat$country <- countries
ddf <- dat

#preliminary plots:
plot_pure_node(var[1]) 

#USA and Canada combined:
plot_pure_node(c(var[1],var[3]),TRUE)
plot_pure_node(var[2])
plot_pure_node(var[3])
plot_pure_node(var[4])

# We now produce a plot for each pure node, which we combine manually in the paper
png("../output/Pure_node_USA.png",units='in',width=8,height=5,res=300)
plot_pure_node(var[1])
dev.off()

png("../output/Pure_node_Japan.png",units='in',width=8,height=5,res=300)
plot_pure_node(var[2])
dev.off()

png("../output/Pure_node_Canada.png",units='in',width=8,height=5,res=300)
plot_pure_node(var[3])
dev.off()

png("../output/Pure_node_Germany.png",units='in',width=8,height=5,res=300)
plot_pure_node(var[4])
dev.off()

# we also produce figure 7:
png("../output/Pure_node_USACanada.png",units='in',width=8,height=5,res=300)
plot_pure_node(c(var[1],var[3]),TRUE)
dev.off()

