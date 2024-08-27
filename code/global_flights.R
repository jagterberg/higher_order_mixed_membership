############################
#' This script contains the code to analyze the global flight dataset from the 
#' supplementary material.  In particular, it contains the code to produce figure 3.
##########################
source("misc.R")
source("HOOI.R")
source("membership_estimation.R")
library(reshape)
library(dplyr)

##############################
# Process data and obtain memberships
################################
load('../data/flight_route.RData')
air_tensor <- as.tensor(air_tensor)
rs <- c(5,5,5)
Uhats <- HOOI_dd(air_tensor,r=rs,niter=20)
clusts <- SPAMM(uhats = Uhats,threshold=T,tval=10^(-10))
round(clusts[[1]][[2]],3)
round(clusts[[1]][[1]],3)

# obtain the pure nodes
purelines <- pick_linename[clusts[[2]][[1]]]
pureports <- pick_portname[clusts[[2]][[2]]]

##############################
# Examine the memberships for the airlines
##############################

# look at pure nodes
select_airl_info[select_airl_info$V4 %in% purelines,]$V2

# associate the row indices to their names
rownames(clusts[[1]][[1]]) <- select_airl_info[order(select_airl_info$V4),]$V2

#save the memberships for the airlines according to the pure nodes
# (note: hardcoded because of ordering)
purelines
cols <- c("united","US airways","british airways","delta","air china")
colnames(clusts[[1]][[1]]) <- cols

# final airline memberships:
round(clusts[[1]][[1]],3)

#look at just chinese airlines and american airlines separately
america_airport_membs <- round(clusts[[1]][[2]][which(rownames(clusts[[1]][[2]]) %in% 
                                                        select_airport_info[select_airport_info$V4 == 'United States',]$V2), ],3)
china_airport_membs <-  round(clusts[[1]][[2]][which(rownames(clusts[[1]][[2]]) %in% 
                                                       select_airport_info[select_airport_info$V4 == 'China',]$V2), ],3)

###########################
#Examine the memberhsips for the airports
###########################

# first find the row names for the airports  
rownames(clusts[[1]][[2]]) <- select_airport_info[order(select_airport_info$V5),]$V2

#examine pure nodes:
select_airport_info[select_airport_info$V5 %in% pureports,]$V2

# we now save the column names according to the pure node:
#note: hardcoded due to ordering
pureports
rowz <- c("london","atlanta","chicago","Beijing","Newark")
colnames(clusts[[1]][[2]]) <- rowz

# view the modified data
round(clusts[[1]][[2]],3)

# look at US and China memberships
china_airline_membs <- round(clusts[[1]][[1]][which(rownames(clusts[[1]][[1]]) %in% 
                                                      select_airl_info[select_airl_info$V7 == 'China',]$V2), ],3)
america_airline_membs <- round(clusts[[1]][[1]][which(rownames(clusts[[1]][[1]]) %in% 
                                                        select_airl_info[select_airl_info$V7 == 'United States',]$V2), ],3)
###########################
# plot airport memberships in each community
##########################

# make a data frame from the memberships
dat <- as.data.frame(clusts[[1]][[2]])
dat$country <- "Other"
dat$country[which(rownames(clusts[[1]][[2]]) %in% 
                    select_airport_info[select_airport_info$V4 
                                        == 'United States',]$V2)] <- "USA"
dat$country[which(rownames(clusts[[1]][[2]]) %in% 
                    select_airport_info[select_airport_info$V4
                                        == "China",]$V2)]  <- "China"
dat$airport <- rownames(clusts[[1]][[2]])

dat_new <- melt(data = dat,ids = c("country","airport"),variable_name = "community")
dat_new2 <- dat_new %>% group_by(country,`community`) %>%
  summarise(`Mean Membership` = mean(value))

# produce Figure 3 for the airports:
g1 <- ggplot(dat_new2, aes( x=country,fill=`community`,y=`Mean Membership`)) + 
  geom_bar(position="stack", stat="identity") + 
  ggtitle("Mean Community Membership: Airports")

png("../output/worldwide_flights_1.png",units='in',width=5,height=5,res=300)
g1
dev.off()

###################################
# plot airline memberships
###################################

dt <- as.data.frame(clusts[[1]][[1]])
dt$country <- "Other"
dt$country[which(rownames(clusts[[1]][[1]]) %in% 
        select_airl_info[select_airl_info$V7 == 
         'United States',]$V2)] <- "USA"
dt$country[which(rownames(clusts[[1]][[1]]) %in% 
        select_airl_info[select_airl_info$V7 == 
         'China',]$V2)] <- "China"
dt$airline <- rownames(clusts[[1]][[1]])

dt <- melt(dt,ids=c("country","airline"),variable_name = 'community')  
dt2 <- dt %>% group_by(country,`community`) %>%
  summarise(`Mean Membership` = mean(value))

# produce figure 3 for the airlines (right)
g2 <- ggplot(dt2, aes( x=country,fill=`community`,y=`Mean Membership`)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Mean Community Membership: Airlines")

png("../output/worldwide_flights_2.png",units='in',width=5,height=5,res=300)
g2
dev.off()

################################
# Additional Analysis for the airlines (not in paper)
################################
#united airlines community (whose membership is >.5)
clus1 <- pick_linename[which(clusts[[1]][[1]][,1] > .5)]
select_airl_info[select_airl_info$V4 %in% clus1,]$V2

#second airline is US airways
clus2 <- pick_linename[which(clusts[[1]][[1]][,2] > .5)]
select_airl_info[select_airl_info$V4 %in% clus2,]$V2

#third is british airways
clus3 <- pick_linename[which(clusts[[1]][[1]][,3] > .5)]
select_airl_info[select_airl_info$V4 %in% clus3,]$V2

#fourth is delta
clus4 <- pick_linename[which(clusts[[1]][[1]][,4] > .5)]
select_airl_info[select_airl_info$V4 %in% clus4,]$V2

#fifth is air china
clus5 <- pick_linename[which(clusts[[1]][[1]][,5] > .5)]
select_airl_info[select_airl_info$V4 %in% clus5,]$V2
clusts[[1]][[1]][which(rownames(clusts[[1]][[1]]) %in% select_airl_info[select_airl_info$V4 %in% clus5,]$V2),]

#####################################
# Additional  Analysis for the airports (not in paper)
#####################################

#first cluster mostly corresponds to london
ports1 <- pick_portname[which(clusts[[1]][[2]][,1] > .5)]
select_airport_info[select_airport_info$V5 %in% ports1,]$V2

#second cluster corresponds primarily to high-degree or Atlanta
ports2 <- pick_portname[which(clusts[[1]][[2]][,2] > .5)]
select_airport_info[select_airport_info$V5 %in% ports2,]$V2

#third cluster corresponds primarily to USA, but contains many airports
ports3 <- pick_portname[which(clusts[[1]][[2]][,3] > .5)]
select_airport_info[select_airport_info$V5 %in% ports3,]$V2

#fourth airport cluster corresponds to PEK which is Beijing.  
ports4 <- pick_portname[which(clusts[[1]][[2]][,4] > .5)]
select_airport_info[select_airport_info$V5 %in% ports4,]$V2

#fifth cluster is basically just miscellaneous
ports5 <- pick_portname[which(clusts[[1]][[2]][,5] > .9)]
select_airport_info[select_airport_info$V5 %in% ports5,]$V2

