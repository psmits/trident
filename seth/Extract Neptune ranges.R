


# read in taxon csv files from download

#rads <- read.csv("~/Dropbox/Neptune_database.Feb.15.2016/Radiolaria_2016-02-15_17-25-56.csv",header=TRUE)
#nannos <- read.csv("~/Dropbox/Neptune_database.Feb.15.2016/Nannos_2016-02-15_18-49-09.csv",header=TRUE)
#forams <- read.csv("~/Dropbox/Neptune_database.Feb.15.2016/Foraminifera_2016-02-15_18-04-51.csv",header=TRUE)
#forams$Fossil.Group <- rep("F",nrow(forams))
#dias <- read.csv("~/Dropbox/Neptune_database.Feb.15.2016/Diatoms_2016-02-15_19-31-50.csv",header=TRUE)

#neptune <- rbind(rads,nannos,forams,dias)

#saveRDS(neptune,"~/Dropbox/Neptune_database.Feb.15.2016/Neptune.Raw.rds")

nep <- readRDS("~/Dropbox/Neptune_database.Feb.15.2016/Neptune.Raw.rds")

## define bins
nep$age.top <- floor(nep$Age..Ma..Gradstein.et.al..2012)
nep$age.bottom <- ceiling(nep$Age..Ma..Gradstein.et.al..2012)
nep$age.bin <- nep$age.bottom - .5
nep$resolved_fossil_name <- paste(nep$Resolved.Genus,nep$Resolved.Species)

# load packages
library(fields)
library(plyr)
library(PBSmapping)
library(gdata)

# Load the equal-area grid:
load("~/Documents/Databases/PBDB marine inverts May 2014/global_45x14.rda")
LatSeq <- global_45x14$latitude
LongSeq <- global_45x14$longitude
polys<- PBSmapping::makeGrid (x=round(LongSeq,0), y=round(LatSeq,0), projection="UTM")


### to be safe, remove occurrences that fall off of the equal-area grid
nep <- nep[(nep$Latitude  < 82) & (nep$Latitude  > -75) ,]


### now compile the first and last occurrence of each genus to merge back in at the end
spec.FA.bottom <- function (df) max(df$Age..Ma..Gradstein.et.al..2012)
spec.LA.top <- function(df) min(df$Age..Ma..Gradstein.et.al..2012)
spec.range <- ddply(nep,.(resolved_fossil_name),each(spec.FA.bottom,spec.LA.top))


ranges <- function(df) {
  #df <- nep[(nep$age.bin  == 15.25) ,]
  int <- mean(df$age.bin)
  #df <- nep[(nep$age.bin >= int - .5 & nep$age.bin <= int + .5),]
  
  if (length(na.omit(df$Latitude)) > 1) {
    events <- data.frame(EID = 1:length(df$Latitude),
                         X = as.numeric(as.vector(df$Longitude)),
                         Y = as.numeric(as.vector(df$Latitude)))
    events <- PBSmapping::as.EventData(na.omit(events), projection = "UTM")
    fc <- findCells(events, polys)
    eac <- length(unique(subset(fc, select = c(PID, SID)))[,
                                                           1])
    gcd <- max(rdist.earth(events[, 2:3], mile = FALSE))
    MaxLat <- max(df$Latitude)
    MinLat <- min(df$Latitude)
    MaxAbsLat <- max(abs(df$Latitude))
    MinAbsLat <- min(abs(df$Latitude))
    MeanAbsLat <- mean(abs(df$Latitude))
    AbsLatRange <-  MaxAbsLat - MinAbsLat
    occurrences <- length(df$Latitude)
  }
  else {
    eac = 1
    gcd = 1
    MaxLat <- max(df$Latitude)
    MinLat <- min(df$Latitude)
    MaxAbsLat <- max(abs(df$Latitude))
    MinAbsLat <- min(abs(df$Latitude))
    MeanAbsLat <- mean(abs(df$Latitude))
    AbsLatRange <-  MaxAbsLat - MinAbsLat
    occurrences <- length(df$Latitude)
  }
  return(data.frame(occurrences, eac, gcd, MaxAbsLat,
                    MinAbsLat, MeanAbsLat, AbsLatRange))
}

#####  now apply it to the PBDB data
bySpec.Neptune <- ddply(nep, .(Fossil.Group,resolved_fossil_name,age.top,age.bottom,age.bin),ranges,.progress="text")


colnames(bySpec.Neptune) <- c("group","species","top","bottom","bin.med","occurrences","eac","gcd","MaxAbsLat","MinAbsLat","MeanAbsLat","AbsLatRange")

Neptune.out <- merge(bySpec.Neptune,spec.range,by.x="species",by.y="resolved_fossil_name")

Neptune.out$species.age <- Neptune.out$spec.FA.bottom - Neptune.out$top 
Neptune.out$Ex <- ifelse(Neptune.out$spec.LA.top >= Neptune.out$top & Neptune.out$spec.LA.top < Neptune.out$bottom ,1,0)
Neptune.out$Or <- ifelse(Neptune.out$spec.FA.bottom <= Neptune.out$bottom & Neptune.out$spec.FA.bottom > Neptune.out$top,1,0)

saveRDS(Neptune.out,"~/Dropbox/Neptune_database.Feb.15.2016/Neptune.processed.1myr.bins.rds")

num.spec <- function(df) length(df$Ex)
num.ex <- function(df) sum(df$Ex)
exes <- ddply(Neptune.out,.(top),each(num.spec,num.ex))
exes$prop.ex <- exes$num.ex/exes$num.spec
