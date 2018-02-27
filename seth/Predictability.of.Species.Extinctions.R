
library(plyr)
library(ggplot2)
library(gbm)
library(pROC)
library(dismo)
setwd(dir="~/Dropbox/Neptune_database.Feb.15.2016")

#neptune.data <- readRDS("~/Dropbox/Rise and Fall IPC 2014/Neptune.processed.stages.rds")
#neptune.data <- readRDS("~/Dropbox/Rise and Fall IPC 2014/Neptune.processed.rds")
neptune.data <- readRDS("~/Dropbox/Neptune_database.Feb.15.2016/Neptune.processed.1myr.bins.rds")
#neptune.data <- readRDS("~/Dropbox/Rise and Fall IPC 2014/Neptune.processed.1myr.bins.rds")
#resolution <- .1


nep.dat <- unique(neptune.data[c("species","spec.FA.bottom","spec.LA.top","group")])
nep.dat$duration <- nep.dat$spec.FA.bottom - nep.dat$spec.LA.top
quartz()
p <- ggplot(nep.dat,aes(duration))
p + geom_histogram() + facet_wrap(~group)


extinctions.per.bin <- function(df) sum(df$Ex)
ext.prop <- function(df) sum(df$Ex)/length(df$Ex)
exnums <- ddply(neptune.data,.(top,group),each(extinctions.per.bin,ext.prop ))


### calculate gapiness
prop.comp <-function(df){
  #df <- neptune.data[(neptune.data$species == "Abathomphalus intermedius"),]
  obs <-  nrow(df)*2
  exp <- max(df$bottom) - min(df$top)
  prop <- obs/exp 
  return(prop)
}
gaps <- ddply(neptune.data,.(species),each(prop.comp),.progress="text")

neptune.data <- merge(neptune.data,gaps)

# cull based on completeness
neptune.data <- neptune.data[(neptune.data$prop.comp >= .5),]

## standardize variables

max.occs <- function(df) max(df$occurrences)
max.eac <- function(df) max(df$eac)
maxoccs <- ddply(neptune.data,.(top),each(max.occs,max.eac))
neptune.data <- merge(neptune.data,maxoccs)
neptune.data$prop.occs <- neptune.data$occurrences/neptune.data$max.occs
neptune.data$prop.eac <- neptune.data$eac/neptune.data$max.eac

# limit data range
neptune.data <- neptune.data[(neptune.data$bottom < 75),]
#subset desired columns
mic1 <- neptune.data
mic1$dataset <- rep("Pelagic Microfossil Species",nrow(mic1))
mic1$resolution <- rep("Sp.",nrow(mic1))
mic1 <- mic1[(mic1$occurrences > 1),]


#rename columns and bind datasets together
#colnames(mic1) <- c("group","taxon","bottom","top","occurrences","eac","gcd","richness","MaxAbsLat","MinAbsLat","MeanAbsLat", "AbsLatRange","taxon.age","Ex","Or","dataset","resolution")


top <- data.frame(table(neptune.data$top))
Bin_Num <- seq(1,nrow(top),1)
bins <- data.frame(top,Bin_Num)
colnames(bins) <- c("top","N.Spec.","Bin_Num")

mic2 <- merge(mic1,bins)
mic2 <- mic2[(mic2$top != 0 & mic2$Bin_Num != 1),]
#mic2$comp.interval <- 

make.rocs <- function(df){
   # df <- mic2[(mic2$Bin_Num == 4 & mic2$group == "R"),]
    mod.int <- mic2[(mic2$Bin_Num == df$Bin_Num[1] & mic2$group == df$group[1]),]
    #mod.int$Ex <- sample(mod.int$Ex,nrow(mod.int),replace=FALSE)
    prev.int <- mic2[(mic2$Bin_Num > df$Bin_Num[1] & mic2$Bin_Num < df$Bin_Num[1]+6 & mic2$group == df$group[1]),]
    
    # make vector of weights
    obs <- length(prev.int$Ex)
    exes <- sum(prev.int$Ex)
    ExFreq <- (exes/obs)
    SurFreq <- (1-(exes/obs))
    MaxFreq <- max(ExFreq,SurFreq)
    ExWeight <- 1/(ExFreq/MaxFreq)
    SurWeight <- 1/(SurFreq/MaxFreq)
    prev.int.weights <- ifelse(prev.int$Ex==1,ExWeight,SurWeight)

   
prev.gbm <- glm(Ex ~ prop.eac + gcd + prop.occs + AbsLatRange + species.age, data = prev.int, family = binomial(link=logit),weights=prev.int.weights)
      
      # prev.pred <- predict(prev.gbm,type="response",n.trees = prev.gbm$n.trees,newdata=mod.int)

      # randomize extinctions if desired for comparison
      #mod.int$Ex <- sample(mod.int$Ex,nrow(mod.int),replace=FALSE)

      prev.pred <- predict(prev.gbm,type="response",newdata=mod.int)
      #prev.pred <- predict(prev.gbm,type="response",newdata=prev.int)
      predictions <- cbind(mod.int,prev.pred)
      
      ROC <- roc(predictions$Ex ~ predictions$prev.pred)
      sensitivity <- ROC$sensitivities
      specificity <- ROC$specificities
      taxa <- rep(length(predictions$Ex),length(sensitivity))
      ex.prop <- rep(mean(predictions$Ex),length(sensitivity))
      ex.num <- rep(sum(predictions$Ex),length(sensitivity))
      AUC <- rep(ROC$auc,length(sensitivity))
      sequence <- seq(1,length(sensitivity))
      top <- rep(predictions$top[1],length(sensitivity))
      output <- data.frame(top,sensitivity,specificity,taxa,ex.prop,ex.num,AUC,sequence)
      return(output)}

output <- ddply(mic2,.(group,Bin_Num),failwith(NULL,make.rocs))
output <- output[(output$ex.num > 4 & output$top < 20),]

quartz()
p <- ggplot(data=output ,aes(1-specificity,sensitivity,colour=group))
p + geom_abline(colour="gray") + geom_path(size=.5) + facet_wrap(~top,ncol=3) + coord_equal() + theme_bw() + xlab("False extinction prediction rate")+ ylab("True extinction prediction rate") + theme(axis.text.x  = element_text(angle=-90, vjust=0.5))