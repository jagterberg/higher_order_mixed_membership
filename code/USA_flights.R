###########################
#' This script contains the code to analyze the USA Flight network from the supplementary 
#' material.  In particular, this script contains the code to produce Figure 4, 
#' # Figure 5, and Figure 6.
###########################

source("misc.R")
source("HOOI.R")
source("membership_estimation.R")
source("estimate_rank.R")
source("pltfunctions.R")
library(airportr)
require(geonames)
require(maps)
library(ggplot2)

#################################
#process data
##################################

#first produce the network as a tensor
load("../data/US_airport_networks-only48states.RData")
air_tensor <- array(0,c(69,343,343))
for (i in c(1:length(Adj_list_new))) {
  air_tensor[i,,] <- as.matrix(Adj_list_new[[i]])
}
air_tensor <- as.tensor(air_tensor)

# we now edit the names slightly for formatting using airport_lookup
airport_names <- rep("",length(lcc))
for (i in 1:length(lcc)) {
  if (i == 285) {
    airport_names[i] <- "Grant County"
  } else {
    airport_names[i] <- airport_lookup(lcc[i],output_type="city")
  }
  
}

################################
#analyze data
#################################

# select 1st dimension
#first try to embed by deleting the diagonal:
mat1 <- k_unfold(air_tensor,1)@data
toembed <- mat1 %*% t(mat1)
diag(toembed) <- 0
embed1 <- eigen(toembed)
togetelbows <- embed1$values[embed1$values >= 0]
togetelbows #only 2 positive values
r1 <- getElbows(togetelbows) 
r1 #results in only two positive eigenvalues.

#try again using vanilla svd scree plot
mat1 <- k_unfold(air_tensor,1)@data
embed1 <-svd(mat1)
togetelbows <- embed1$d# embed1$values[embed1$values >= 0]
r1 <- getElbows(togetelbows)
r1 #1 4 18
#don't want to choose 1, so choose 4 for interpretability
r1 <- r1[2]


#get second embedding dimension
mat2 <- k_unfold(air_tensor,2)@data
toembed <- mat2 %*% t(mat2)
diag(toembed) <- 0
embed1 <- eigen(toembed)
togetelbows <- embed1$values[embed1$values >= 0]
togetelbows <- sqrt(togetelbows)#only look at postitive
r2 <- getElbows(togetelbows)
r2 #1 3 9 
r2 <- r2[2] #keep r2 = 3

# by symmetry we set the third mode to have the same rank as the second
rs <- c(r1,r2,r2)

#find embedding and estimate memberships
Uhats <- HOOI_dd(air_tensor,r= rs, niter=20)
clusts <- SPAMM(uhats = Uhats,threshold=T,tval=10^(-20))

##########################
# Analyze results and produce figures
############################

#process data for plotting
rownames(clusts[[1]][[2]]) <- lcc
dat <- as.data.frame(clusts[[1]][[2]])
names(dat) <- paste("Pure_node",lcc[clusts[[2]][[2]]],sep="_")
dat$airport <- airport_names
dat$airport2 <- lcc

# we will load the latitude and longitude data for plotting
cities <- read.csv("../data/uscities.csv")
latitude <- rep(0,nrow(dat))
longitude <- rep(0,nrow(dat))
for (i in c(1:length(latitude))) {
  if(any(cities$city == dat$airport[i])) {
    latitude[i] <- cities[which(cities$city == dat$airport[i])[1],"lat"]
    longitude[i] <- cities[which(cities$city == dat$airport[i])[1],"lng"]
    
  # a few don't work so we do it manually:
  } else if (dat$airport[i] == "Dallas-Fort Worth") {
    latitude[i] <- cities[which(cities$city == "Dallas")[1],"lat"]
    longitude[i] <- cities[which(cities$city == "Dallas")[1],"lng"]
    
  } else if (dat$airport[i] == "Sault Ste Marie" ) {
    latitude[i] <- 46.5136
    longitude[i] <- 84.3476
  } else if (dat$airport[i] == "Wassau") {
    latitude[i] <- 44.9591
    longitude[i] <- 89.6301
  } else if (dat$airport[i] == "Montrose CO") {
    latitude[i] <- 38.4783
    longitude[i] <-  107.8762
  } else if (dat$airport[i] ==   "Windsor Locks" ) {
    latitude[i] <-41.9243
    longitude[i] <-  72.6454
  } else if (dat$airport[i] ==  "Manchester NH" ) {
    latitude[i] <-42.9956
    longitude[i] <-  71.4548
  } else if (dat$airport[i] ==  "Raleigh-durham"  ) {
    latitude[i] <- 35.8992
    longitude[i] <-  78.8636
  } else if (dat$airport[i] ==  "Charlottesville VA" ) {
    latitude[i] <-38.0293
    longitude[i] <- 78.4767
  } else if (dat$airport[i] ==  "Lexington KY") {
    latitude[i] <-38.0406
    longitude[i] <-  84.5037
  } else if (dat$airport[i] ==  "Roanoke VA" ) {
    latitude[i] <-37.2710
    longitude[i] <-  79.9414
  } else if (dat$airport[i] ==  "BRISTOL") {
    latitude[i] <-36.5951
    longitude[i] <-  82.1887
  } else if (dat$airport[i] ==  "Jacksn Hole" ) {
    latitude[i] <-43.4799
    longitude[i] <-  110.7624
  } else if (dat$airport[i] ==  "Bush Field"  ) {
    latitude[i] <-33.3700
    longitude[i] <-  81.9652
  } else if (dat$airport[i] ==  "Islip"  ) {
    latitude[i] <-40.7298
    longitude[i] <-  73.2104
  } else if (dat$airport[i] ==  "Mcallen") {
    latitude[i] <-26.2034
    longitude[i] <-  98.2300
  } else if (dat$airport[i] ==  "Jacksonville NC" ) {
    latitude[i] <-30.3322
    longitude[i] <-  81.6557
  } else if (dat$airport[i] ==  "State College Pennsylvania") {
    latitude[i] <-40.7934
    longitude[i] <-  77.8600
  } else if (dat$airport[i] ==  "Redmond-Bend" ) {
    latitude[i] <-44.2726
    longitude[i] <- 121.1739
  } else if (dat$airport[i] ==  "MONTGOMERY" ) {
    latitude[i] <-32.3792
    longitude[i] <-  86.3077
  } else if (dat$airport[i] ==  "Hattiesburg/Laurel"  ) {
    latitude[i] <-31.4682
    longitude[i] <-  89.3354
  } else if (dat$airport[i] ==  "Arcata CA" ) {
    latitude[i] <-40.8665
    longitude[i] <- 124.0828
  } else if (dat$airport[i] ==   "Vineyard Haven MA" ) {
    latitude[i] <-41.4543
    longitude[i] <-  70.6036
  } else if (dat$airport[i] ==  "Barnstable") {
    latitude[i] <-41.7003
    longitude[i] <-  70.3002
  } else if (dat$airport[i] ==  "PARKERSBURG"  ) {
    latitude[i] <-39.2667
    longitude[i] <-  81.5615
  } else if (dat$airport[i] ==  "Devils Lake" ) {
    latitude[i] <-48.1128
    longitude[i] <-  98.8651
  } else if (dat$airport[i] ==  "Grant County") {
    latitude[i] <-40.4731
    longitude[i] <-  85.6846
  } else if (dat$airport[i] ==  "Riverton WY" ) {
    latitude[i] <-43.0247
    longitude[i] <-  108.3806
  } else if (dat$airport[i] ==  "Nantucket" ) {
    latitude[i] <-41.2835
    longitude[i] <-  70.0995
  } else if (dat$airport[i] ==  "PADUCAH"  ) {
    latitude[i] <-37.0834
    longitude[i] <-  88.6000
  } else if (dat$airport[i] ==  "Columbus Mississippi" ) {
    latitude[i] <-33.4957
    longitude[i] <-  88.4273
  } else if (dat$airport[i] ==  "Bar Harbor"  ) {
    latitude[i] <-44.3876
    longitude[i] <-  68.2039
  } else {
    latitude[i] <- NA
    longitude[i] <- NA
  }
  
}
dat$LAT <- latitude
dat$LON <- -abs(longitude)
dat[which(is.na(dat$lat)),"airport"]
purenode_airports <- lcc[clusts[[2]][[2]]]
purenode_airports2 <- airport_names[clusts[[2]][[2]]]
dat$purenode1 <- ifelse(dat$airport2 %in% purenode_airports,40,20)
dat$purenode2 <- as.factor(ifelse(dat$airport2 %in% purenode_airports,"100","0"))
purenode_airports2 # list of pure nodes

# Draw maps of airports (figure 5)
gs <- list() # list of plots

# for each pure node, we plot the memberships and manually format for the paper
for (i in c(1:length(purenode_airports))) {
  gs[[i]] <- make_US_plot(dat,purenode_airports[i]) #calls from the pltfunctions.R script
  flnme <- paste0("../output/pure_node_",purenode_airports[i],".png")
  png(flnme,units='in',width=5,height=5,res=300)
  print(gs[[i]])
  dev.off()
}
gs[[1]]
gs[[2]]
gs[[3]]

# Draw time plots (figure 4)
start_v <- as.Date("2016-01-01")
dtvals <- seq.Date(from = start_v, length.out = 69, by = "month")
plts <- list()
# for each pure node we make a smoothed plot
for ( i in c(1:rs[1])) {
  dt <- dtvals[clusts[[2]][[1]][i]]
  plts[[i]] <- make_time_plot(i,tt= paste0("Pure Node:",dt)) #uses pltfunctions.R
  flnme <- paste0("../output/pure_node_",dt,".png")
  png(flnme,units='in',width=5,height=5,res=300)
  print(plts[[i]])
  dev.off()
}

# we also make figure 6 by combining pure nodes
png("../output/janmarch.png",units='in',width=5,height=5,res=300)
make_time_plot(c(1,3),combine=T,tt="January 2021 and March 2020")
dev.off()
png("../output/janaug.png",units='in',width=5,height=5,res=300)
make_time_plot(c(1,2),combine=T,tt="January 2021 and August 2021")
dev.off()


