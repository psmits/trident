library(plyr)
library(ggplot2)
library(gbm)
library(pROC)
library(irr)
setwd(dir="~/Dropbox/Rise and Fall IPC 2014")

#neptune.data <- readRDS("~/Dropbox/Rise and Fall IPC 2014/Neptune.processed.stages.rds")
#neptune.data <- readRDS("~/Dropbox/Rise and Fall IPC 2014/Neptune.processed.rds")
neptune.data <- readRDS("~/Dropbox/Rise and Fall IPC 2014/Neptune.processed.2myr.bins.rds")
#neptune.data <- readRDS("~/Dropbox/Rise and Fall IPC 2014/Neptune.processed.1myr.bins.rds")
#resolution <- .1

extinctions.per.bin <- function(df) sum(df$Ex)
ext.prop <- function(df) sum(df$Ex)/length(df$Ex)
exnums <- ddply(neptune.data,.(top,group),each(extinctions.per.bin,ext.prop ))


### calculate gapiness
prop.comps <-function(df){
  # df <- neptune.data[(neptune.data$species == "Abathomphalus intermedius"),]
  obs <-  nrow(df)
  exp <- max(df$bottom) - min(df$top)
  prop.comp <- obs/exp 
  total.occs <- sum(df$occurrences)
  out <- data.frame(prop.comp,total.occs)
  return(out)
}
gaps <- ddply(neptune.data,.(species),each(prop.comps),.progress="text")

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
neptune.data <- neptune.data[(neptune.data$bottom < 66),]
#subset desired columns
mic1 <- neptune.data
mic1$dataset <- rep("Pelagic Microfossil Species",nrow(mic1))
mic1$resolution <- rep("Sp.",nrow(mic1))
mic1 <- mic1[(mic1$total.occs > 5),]


#rename columns and bind datasets together
#colnames(mic1) <- c("group","taxon","bottom","top","occurrences","eac","gcd","richness","MaxAbsLat","MinAbsLat","MeanAbsLat", "AbsLatRange","taxon.age","Ex","Or","dataset","resolution")


top <- data.frame(table(neptune.data$top))
Bin_Num <- seq(1,nrow(top),1)
bins <- data.frame(top,Bin_Num)
colnames(bins) <- c("top","N.Spec.","Bin_Num")

mic2 <- merge(mic1,bins)
mic2 <- mic2[(mic2$top != 0),]

### order and then attach past and future ranges
mic3 <- mic2[order(mic2$species, -mic2$top),]
mic4 <- mic3[c("prop.occs","prop.eac","gcd","AbsLatRange")]

past <- mic4[1:nrow(mic4)-1,]
colnames(past) <- c("prop.occs.p","prop.eac.p","gcd.p","AbsLatRange.p")
first <- data.frame(t(c(NA,NA,NA,NA)))
colnames(first) <- c("prop.occs.p","prop.eac.p","gcd.p","AbsLatRange.p")
past <- rbind(first,past)

mic5 <- cbind(mic3,past)
mic5$prop.occs.p <- ifelse(mic5$Or == 1,0,mic5$prop.occs.p)
mic5$prop.eac.p <- ifelse(mic5$Or == 1,0,mic5$prop.eac.p)
mic5$gcd.p <- ifelse(mic5$Or == 1,0,mic5$gcd.p)
mic5$AbsLatRange.p <- ifelse(mic5$Or == 1,0,mic5$AbsLatRange.p)

mic5$prop.occs.diff <- mic5$prop.occs - mic5$prop.occs.p
mic5$prop.eac.diff <- mic5$prop.eac - mic5$prop.eac.p
mic5$gcd.diff <- mic5$gcd - mic5$gcd.p
mic5$AbsLatRange.diff <- mic5$AbsLatRange - mic5$AbsLatRange.p

mic6 <- mic5[c("Ex","prop.occs.diff","prop.eac.diff","gcd.diff","AbsLatRange.diff","prop.occs","prop.eac","gcd","AbsLatRange","species.age")]

model <- glm(Ex ~., data = mic6 , family = binomial) %>%
  stepAIC(trace = FALSE)
summary(model)

mic5 <- mic5[order(-mic5$top),]

randEx <- mic5[c("Ex","top")]
randEx$rand <- runif(nrow(randEx))
randEx <- randEx[order(-randEx$top,randEx$rand),]
#mic5$Ex <- randEx$Ex



comps <- function(df){
  df <- mic5[(mic5$top==6 & mic5$group == "Coccolithophora"),]
  #df <- mic5[(mic5$Bin_Num==10),]
  #mod.int <- mic5[(mic5$Bin_Num == df$Bin_Num[1] & mic5$group == df$group[1]),]
  mod.int <- mic5[(mic5$top == df$top[1]),]
  #mod.int$Ex <- sample(mod.int$Ex,nrow(mod.int),replace=FALSE)
  #next.int <- mic5[(mic5$Bin_Num == df$Bin_Num[1]+1 & mic5$group == df$group[1]),]
  next.int <- mic5[(mic5$top == df$top[1]-2),]

  num.ex <- sum(mod.int$Ex)
  num.pred.ex <- sum(next.int$Ex)

  int.past=gbm(Ex ~ prop.occs.diff + prop.eac.diff + gcd.diff + AbsLatRange.diff + prop.occs + prop.eac +  gcd + AbsLatRange + species.age +group  ,data = mod.int,distribution = "bernoulli",n.trees = 10000,shrinkage = 0.01, interaction.depth = 4,bag.fraction = 0.5,train.fraction = .67,cv.folds=3)

  int.nopast=gbm(Ex ~ prop.occs + prop.eac +  gcd + AbsLatRange + species.age +group  ,data = mod.int,distribution = "bernoulli",n.trees = 10000,shrinkage = 0.01, interaction.depth = 4,bag.fraction = 0.5,train.fraction = .67,cv.folds=3)

  past.best <- gbm.perf(int.past, method="OOB", plot.it=FALSE, oobag=TRUE, overlay=TRUE)
  nopast.best <- gbm.perf(int.nopast, method="OOB", plot.it=FALSE, oobag=TRUE, overlay=TRUE)

  int.past=gbm(Ex ~ prop.occs.diff + prop.eac.diff + gcd.diff + AbsLatRange.diff + prop.occs + prop.eac +  gcd + AbsLatRange + species.age +group  ,data = mod.int,distribution = "bernoulli",n.trees = past.best,shrinkage = 0.01, interaction.depth = 4,bag.fraction = 0.5,train.fraction = .67,cv.folds=3)

  int.nopast=gbm(Ex ~ prop.occs + prop.eac +  gcd + AbsLatRange + species.age +group  ,data = mod.int,distribution = "bernoulli",n.trees = nopast.best,shrinkage = 0.01, interaction.depth = 4,bag.fraction = 0.5,train.fraction = .67,cv.folds=3)

  past.auc <- auc(roc(mod.int$Ex,predict(int.past,n.trees = int.past$n.trees,type="response")))
  past.auc.ci <- ci.auc(roc(mod.int$Ex,predict(int.past,n.trees = int.past$n.trees,type="response")))
  past.auc.ci <- past.auc.ci[3] - past.auc.ci[2]
  nopast.auc <- auc(roc(mod.int$Ex,predict(int.nopast,n.trees = int.nopast$n.trees,type="response")))
  nopast.auc.ci <- ci.auc(roc(mod.int$Ex,predict(int.nopast,n.trees = int.nopast$n.trees,type="response")))
  nopast.auc.ci <- nopast.auc.ci[3] - nopast.auc.ci[2]
  past.pred.auc <- auc(roc(next.int$Ex,predict(int.past,newdata=next.int,n.trees = int.past$n.trees,type="response")))
  past.pred.auc.ci <- ci.auc(roc(next.int$Ex,predict(int.past,newdata=next.int,n.trees = int.past$n.trees,type="response")))
  past.pred.auc.ci <- past.pred.auc.ci[3] - past.pred.auc.ci[2]
  nopast.pred.auc <- auc(roc(next.int$Ex,predict(int.nopast,newdata=next.int,n.trees = int.nopast$n.trees,type="response")))
  nopast.pred.auc.ci <- ci.auc(roc(next.int$Ex,predict(int.nopast,newdata=next.int,n.trees = int.nopast$n.trees,type="response")))
  nopast.pred.auc.ci <- nopast.pred.auc.ci[3] - nopast.pred.auc.ci[2]


  out <- data.frame(past.auc,past.auc.ci,nopast.auc,nopast.auc.ci,past.pred.auc,past.pred.auc.ci,nopast.pred.auc,nopast.pred.auc.ci,num.ex,num.pred.ex)
  return(out)}

aucs <- ddply(mic5,.(bin.med,top,bottom),failwith(NULL,comps),.progress="text")
aucs

aucs2 <- aucs[c(aucs$num.ex > 1 & aucs$num.pred.ex > 1),]

quartz()
p <- ggplot(aucs2,aes(past.auc,past.pred.auc,fill=bin.med))
p + geom_abline() + geom_errorbarh(aes(xmax=past.auc + past.auc.ci,xmin=past.auc - past.auc.ci),colour="gray") +  geom_errorbar(aes(ymax=past.pred.auc + past.pred.auc.ci,ymin=past.pred.auc - past.pred.auc.ci),size=.25) + geom_point(size=4,pch=21)  + coord_equal(xlim=c(.5,1),ylim=c(.5,1)) + theme_bw()

quartz()
p <- ggplot(aucs2,aes(nopast.auc,nopast.pred.auc))
p + geom_abline() + geom_errorbarh(aes(xmax=nopast.auc + nopast.auc.ci,xmin=nopast.auc - nopast.auc.ci),colour="gray") +  geom_errorbar(aes(ymax=nopast.pred.auc + nopast.pred.auc.ci,ymin=nopast.pred.auc - nopast.pred.auc.ci),colour="black",size=.25) + geom_point(size=4,pch=21,colour="white",fill="black")  + coord_equal(xlim=c(.5,1),ylim=c(.5,1)) + theme_bw()
