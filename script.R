library(ape)
library(phangorn)
library(adephylo)
library(magrittr)

tr <- read.tree("data/ebola.fas.treefile")
td <- tr$tip.label %>% strsplit(.,"-|_") %>% lapply(.,tail,1) %>% unlist %>% as.double
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
root.time <- results$TMRCA[1]
max.time <- max(td)
calibrating.values$age.min <- max.time - root.time
calibrating.values$age.max <- max.time - root.time
# pins the tips to sampling years
calibrating.values <- rbind(calibrating.values,
                            data.frame(node=seq(1,length(td)),
                                       age.min=max.time - td,
                                       age.max=max.time - td,
                                       soft.bounds=FALSE))

dated.tree <- chronos(tr.rtt, 
                     lambda=1, 
                     model="discrete", 
                     calibration=calibrating.values,
                     control=chronos.control(nb.rate.cat=1)
                     )

lasv.chronos <- tr
lasv.host <- rep("Human",length(lasv.chronos$tip.label))
lasv.host[grep("Josiah",lasv.chronos$tip.label)] <- "Lab"
lasv.host[grep("LM",lasv.chronos$tip.label)] <- "Mastomys"
lasv.host[grep("ZO",lasv.chronos$tip.label)] <- "Mastomys"
lasv.annotations <- data.frame(Taxa=lasv.chronos$tip.label,Time=td,Host=lasv.host)
write.table(lasv.annotations,"data/lasv.txt",col.names=TRUE,row.names=FALSE,sep="\t")

country <- rep("DRC",29)
country[13:19]="SierraLeone"
country[20:22]="Guinea"
country[c(3,26,27,28)]="Gabon"
ebola.annotations <- data.frame(Taxa=tr$tip.label,Time=td,Country=country)
write.table(ebola.annotations,"data/ebola.txt",col.names=TRUE,row.names=FALSE,sep="\t")
