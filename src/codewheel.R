data=read.table(file="data.txt",header=T,sep="\t")

library(reshape)
library(ggplot2)
library(plyr)
library(grDevices)

nba=data 
jpeg("wheel.jpg",width=1200,height=1200) 

nba.m <- melt(nba)
nba.m <- ddply(nba.m, .(variable), transform, value = value)
nba.m$var2 = as.numeric(nba.m$variable) + 35
y_labels = levels(nba.m$variable)
y_breaks = seq_along(y_labels) + 35

nba.labs <- subset(nba.m, variable==levels(nba.m$variable)[nlevels(nba.m$variable)])
nba.labs <- nba.labs[order(nba.labs$Name),]
nba.labs$ang <- seq(from=(360/nrow(nba.labs))/1.5, to=(1.5*(360/nrow(nba.labs)))-360, length.out=nrow(nba.labs))+80
nba.labs$hjust <- 0
nba.labs$hjust[which(nba.labs$ang < -90)] <- 1
nba.labs$ang[which(nba.labs$ang < -90)] <- (180+nba.labs$ang)[which(nba.labs$ang < -90)]
p2 = ggplot(nba.m, aes(x=Name, y=var2, fill=value)) +
     geom_tile(colour="white") +
     geom_text(data=nba.labs, aes(x=Name, y=var2+0.6,
        label=Name, angle=ang, hjust=hjust), size=4) +
     scale_fill_gradient(low = "white", high = "darkred") +
     ylim(c(0, max(nba.m$var2) + 6.6)) +
     scale_y_discrete(breaks=y_breaks, labels=y_labels) +
     coord_polar(theta="x") +
     theme(panel.background=element_blank(),
           axis.title=element_blank(),
           panel.grid=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks=element_blank(),
           axis.text.y=element_text(size=3))
print(p2)

dev.off() 
