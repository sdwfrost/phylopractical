library(ape)
library(phangorn)
library(adephylo)
library(magrittr)
library(ggplot2)
library(ggtree)

tr <- read.tree("data/ebola.fas.treefile")
annotations <- read.table("data/ebola.txt",header=T,row.names=NULL,sep="\t")
td <- as.double(annotations[match(tr$tip.label,annotations[,1]),2])
tr.rtt <- rtt(tr,td,objective="rsquared")
tr.rtt$edge.length[tr.rtt$edge.length<0] <- 0.0
rd <- distRoot(tr.rtt)
td <- tr.rtt$tip.label %>% strsplit(.,"-|_") %>% lapply(.,tail,1) %>% unlist %>% as.double
rtt.lm <- lm(rd~td)

root.time <- unname(-as.double(coef(rtt.lm)[1])/coef(rtt.lm)[2])
results <- data.frame(TMRCA=root.time,SubstRate=as.double(coef(rtt.lm)[2]))
results
plot(rd~td,xlab="Time",ylab="Root to tip distance",ylim=c(0,max(rd)),xlim=c(root.time,max(td)),pch=16,col="red")
abline(rtt.lm)
abline(h=0,lty=2)

calibrating.values <- makeChronosCalib(tr.rtt)
root.time <- 1975
max.time <- max(td)
calibrating.values$age.min <- max.time - root.time
calibrating.values$age.max <- max.time - root.time
# pins the tips to sampling years
calibrating.values <- rbind(calibrating.values,
                            data.frame(node=seq(1,length(td)),
                                       age.min=max.time - td,
                                       age.max=max.time - td,
                                       soft.bounds=FALSE))

dated.tree <- RLchronos(tr.rtt, 
                     lambda=1, 
                     model="discrete", 
                     calibration=calibrating.values,
                     control=chronos.control(nb.rate.cat=1)
                     )
dated.tree2 <- read.tree(text=write.tree(dated.tree))

g <- ggtree(dated.tree2)+theme_tree2()
g <- g %<+% annotations
g + geom_tippoint(aes(color=Country),size=2)+theme(legend.position="right")
