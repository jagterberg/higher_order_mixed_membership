####################################
#' This script contains the auxiliary plotting functions required to analyze
#' the USA flight data from the supplementary material. 
###################################

library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggplot2)

# this function plots the memberships in a given pure node by overlaying onto
# a map of the USA
make_US_plot <- function(dat,pure) {
  var <- paste("Pure_node",pure,sep="_")
  dat <- dat[order(dat[,var],decreasing = FALSE),]
  dat$pure <- dat[,var]
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))
  states <- map_data("state") 
  us <- map_data('world','usa')
  g1 <- ggplot() + geom_polygon(data=us, aes(x=long, y=lat,group=group), colour="grey20", fill="grey80")+
    geom_polygon(data=states,aes(x=long,y=lat,group=group),color="black",alpha=0) +
    geom_point(data = dat,aes(x=LON, y = LAT, 
                              color = pure,
                              shape=purenode2,
                              size=purenode1
    ),alpha=.85)+
    coord_map(projection = "mercator", xlim=c(-128, -65),ylim=c(25,50))+
    theme_bw() + sc +theme(axis.title = element_blank(),axis.text = element_blank()) +
    labs(color = paste(pure,"Membership"),purenode=NA)  +
    guides(size=FALSE,shape=FALSE) + scale_size(range=c(2,5))
  g1
}

# this function plots the membership in a pure node as a function of time, and
# it adds a smoothing line.  
make_time_plot <- function(wv,tt=NULL,togroup=TRUE,combine=F) {
  if(!combine) {
    start_v <- as.Date("2016-01-01")
    
    
    df <- data.frame(nothing=c(1:69),
                     membership=clusts[[1]][[1]][,wv],
                     year=c(rep("2016",12),rep("2017",12),rep("2018",12),
                            rep("2019",12),rep("2020",12),rep("2021",9)))
    df$time <- seq.Date(from = start_v, length.out = dim(df)[1], by = "month")
  } else {
    start_v <- as.Date("2016-01-01")
    df <- data.frame(nothing=c(1:69),
                     membership=rowSums(clusts[[1]][[1]][,wv]),
                     year=c(rep("2016",12),rep("2017",12),rep("2018",12),
                            rep("2019",12),rep("2020",12),rep("2021",9)))
    df$time <- seq.Date(from = start_v, length.out = dim(df)[1], by = "month")
    
  }
   if(is.null(tt)) {
      tt <-paste("Community",wv,"by month")
    }
    if(!togroup) {
      ggplot(df,aes(x=time,y=membership)) +
        geom_smooth(method="loess") + 
        geom_point() +
        ggtitle(tt)+ 
        scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
        ylim(0,1)
    } else {
      ggplot(df,aes(x=time,y=membership)) +
        geom_smooth(aes(group=year),method="loess") + 
        geom_point() +
        ggtitle(tt)+ 
        scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
        ylim(0,1)
    }
  
}
