## R Code for beta diversity workshop
## QCBS R symposium 2018
## by Vincent Fugere & Katrine Turgeon

#### prep work ####

rm(list=ls())

#loading data for temporal diversity example (it takes some time)
library(adklakedata)
zoops <- adk_data('crustacean')
pH <- adk_data('chem')

#### Example I: Ugandan streams (spatial beta div) ####
library(RCurl)
library(vegan)
library(gclus)
library(RColorBrewer)
library(ade4)
library(adespatial)

#getting the dataset from GitHub
kib <-read.csv(text=getURL("https://raw.githubusercontent.com/VFugere/QCBS_betadiv/master/kibale_inverts.csv"), header=T)

# #Alternative: getting the dataset from the Open Science Framework (OSF) website
# devtools::install_github('chartgerink/osfr')
# library(osfr)
# osfr::download_files('zc6xt', '/Users/vincentfugere/Desktop/') #change path to working directory
# kib <- read.csv('/Users/vincentfugere/Desktop/kibale_inverts.csv', header=T)
# 
# #if both fail, you can always try Dropbrox! (link might become obsolete eventually)
# https://www.dropbox.com/s/qj8vswqp322add2/kibale_inverts.csv?dl=0

#fixing some issues
colnames(kib)[1] <- 'site.code' #corrects excel generated junk
kib$land.use <- relevel(kib$land.use, 'park') #more intuitive to have park/forest as reference level
#have a look at the data
head(kib)
str(kib)
#will model community composition (columns aeshnidae:veliidae) ~ land.use

#some useful matrices
spe <- kib[,-c(1:3)]
row.names(spe) <- kib$site.code
site <- kib[,2:3]
rownames(site) <- kib[,1]

## 'Classic' community composition & alpha diversity analysis ##

#dissimilarity matrix (log to give more weight to rare species). Chose Bray-Curtis but could use something else
spe.dis <- vegdist(log1p(spe), method='bray')

# do land use categories differ in community composition? clustering of sites using dissimilarity matrix & UPGMA clustering
spe.dis.UPGMA <- hclust(spe.dis, method="average")
spe.dis.o <- reorder.hclust(spe.dis.UPGMA, spe.dis)
dend <- as.dendrogram(spe.dis.o)
com.tx <- apply(t(spe),1,sum) #common taxa for heat map plot
spe.com <- spe[,which(com.tx > 20)]
or <- vegemite(spe.com, spe.dis.o,scale='log')
heatmap(t(spe.com[or$species]), Rowv=NA, Colv=dend, cexRow=0.5, cexCol=0.5, col=c('white', brewer.pal(5,"Greens")), scale="column", revC = T, margin=c(4.5,8))
#farm sites (Kab & Kam) are very different from forest sites (For)

#confirm overall difference withpermutational multivariate analysis of variance
adonis(spe.dis ~ land.use, data=site, permutations=1000)
#yep

# do land use categories differ in alpha diversity?
div <- site[,1:2]
div$abund <- apply(spe,MARGIN=1,FUN=sum)
div$richness <- specnumber(spe)
div$r.richness <- rarefy(spe,min(div$abund)) #uses 77 random individuals, the lowest abundance among samples
div$shannon <- diversity(spe, index='shannon')

par(mfrow=c(2,2), mar=c(4,4,1,1), cex=1)
boxplot(richness ~ land.use, div, ylab="richness",whisklty=1,boxwex=0.5,staplelwd=-1)
boxplot(r.richness ~ land.use, div, ylab="rarefied richness",whisklty=1,boxwex=0.5,staplelwd=-1)
boxplot(shannon ~ land.use, div, ylab=expression("Shannon-Wiener "*italic("H")),whisklty=1,boxwex=0.5,staplelwd=-1)
plot(specaccum(spe), xlab = "streams sampled", ylab = "taxa sampled",ci.type='polygon',ci.col='light gray',ci.lty=0,lwd=2)
plot(specaccum(spe[which(div$land.use == 'park'),]),add=TRUE,ci=0,lwd=2,col='red')
plot(specaccum(spe[which(div$land.use == 'farm'),]),add=TRUE,ci=0,lwd=2,col='blue')
points(x=c(34,11,23),y=c(57,44,48),pch=16,col=c(1,2,4),cex=1.2)
legend('bottom',legend=c('park','farm','both'),bty='n',pt.cex=1.2,pch=16,col=c('red','blue','black'))
par(mfrow=c(1,1))

#beta diversity, at last
beta.BC <- betadisper(spe.dis, site$land.use, type = 'centroid')
plot(beta.BC) #ugly

mybetadivplot <- function (x, cex=1){
  g <- scores(x, choices = 1:2)
  plot(g$sites, type = "n", ylab='Dimension 2', xlab='Dimension 1')
  points(g$sites, pch = c(21,22,21)[site$region], bg=c('white','white','black')[site$region])
  points(g$centroids, pch = 16, cex = 1.5, col = "red")
  for (i in levels(x$group)) {
    ch <- chull(g$sites[x$group == i, ])
    ch <- c(ch, ch[1])
    lines(x$vectors[x$group == i, 1:2][ch, ], col = "black", lty = "dashed")
  }
  legend('topright',cex=1, pt.cex=1, inset = c(0.2,0),legend=c("Kab (farm)","Kam (farm)",'For (park)'),pch=c(21,22,21),pt.bg=c('white','white','black'),bty='n')
}
mybetadivplot(beta.BC)
legend('bottomright',bty='n',inset=c(0.2,0),legend = expression(italic(Bray-Curtis)))

#with incidence-based index instead
dist.jac <- vegdist(spe, method='jaccard', binary=T)
beta.jac <- betadisper(dist.jac, site$land.use, type='centroid')
mybetadivplot(beta.jac)
legend('bottomright',bty='n',inset=c(0.2,0),legend = expression(italic(Jaccard)))

#significant difference?
permutest(beta.jac)
boxplot(beta.jac$distances ~ site$land.use, ylab="Distance to centroid",whisklty=1,boxwex=0.5,staplelwd=-1,frame.plot=F,xaxt='n',ylim=c(.2,.7))
axis(1,at=c(1,2),labels=c('park','farm'),tick=F)

#decompose D into richness distance + replacement distance
park.jac <- beta.div.comp(spe[site$land.use == 'park',],coef='J',quant=F)
farm.jac <- beta.div.comp(spe[site$land.use == 'farm',],coef='J',quant=F)
park.jac

BD.mat <- cbind(park.jac$part[2:3],farm.jac$part[2:3],park.jac$part[4:5],farm.jac$part[4:5])
barplot(BD.mat, names.arg=rep(c('park','farm'),2), col=c(4,2), width=0.5, ylab='contribution to total dissimilarity')
axis(1, lwd=0, lwd.ticks=0, at=c(0.66,1.85), labels = c('absolute','relative'), line = 2, cex = 1.5)
legend('topleft', inset=c(0.1,0), legend=c(expression(paste(Delta,' richness',sep='')),'replacement'), pch=22, pt.bg=c(2,4), bty='n', y.intersp = 2)

par(mfrow=c(1,2))
triang1 = cbind(park.jac$rich,(1-park.jac$D), park.jac$repl)
triangle.plot(as.data.frame(triang1), show=F, addmean=T, labeltriangle=F, scale=F)
triang2 = cbind(farm.jac$rich,(1-farm.jac$D), farm.jac$repl)
triangle.plot(as.data.frame(triang2), show=F, labeltriangle=F, addmean = T, scale=F)

#### Example II: Adirondack lakes (temporal beta div) ####
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)

#getting TBI() function. I DL'd it from Pierre's site and UL'd to GitHub (it might be in a package too)
TBIfunc <- getURL("https://raw.githubusercontent.com/VFugere/QCBS_betadiv/master/TBI.R", ssl.verifypeer = FALSE)
eval(parse(text = TBIfunc))

#loading and formatting dataset. Large-ish dataset so moving to the tidyverse
zoops$Sp <- paste0(zoops$Genus,'.',zoops$Species) %>% as.factor
zoops %<>% select(-PERMANENT_ID, -date, -(month:ug_WWperind), -mgWW.l) %>%
  group_by(lake.name, year, Sp) %>%
  summarize(density = mean(org.l)) %>%
  spread(Sp, density) 

#only 50% of lakes have data up to 2012. Others stop at 2006
with(zoops, table(lake.name, year))
table <- with(zoops, table(lake.name, year))
names(table[table[,18] == 1,18]) -> lakes2kp
zoops %<>% filter(lake.name %in% lakes2kp)

#adding pH
pH %<>% select(lake.name, year, pH) %>%
  group_by(lake.name, year) %>%
  summarize(pH = mean(pH))
zoops <- left_join(zoops, pH, by = c('lake.name','year')) %>%
  select(lake.name, year, pH, everything()) %>%
  as.data.frame

#Pierre Legendres's TBI function

#1) subsetting data to use oldest and most recent time points
com.y1 <- zoops %>% filter(year == 1994) %>% select(Aglaodiaptomus.leptpus:unknown.unknown) %>% as.matrix
com.y2 <- zoops %>% filter(year == 2012) %>% select(Aglaodiaptomus.leptpus:unknown.unknown) %>% as.matrix

#2) setting up a results matrix and adding pH change as an explanatory variable
TBI.mat <- zoops %>% filter(year == 2012) %>% select(lake.name:pH) %>% as.data.frame
TBI.mat$ph.t1 <- zoops$pH[zoops$year == 1994]
TBI.mat$ph.diff <- zoops$pH[zoops$year == 2012] - zoops$pH[zoops$year == 1994]

#3) computing dissimilarity in time and partitioning into constituents
TBI.res <- TBI(com.y1, com.y2, method="jaccard", pa.tr=T, save.BC = T)
TBI.res$TBI #these are Jaccard dissimilarities among the 2 sampling occasions for each site
TBI.res$BC #number of taxa lost (B) and gained (C) over study period
TBI.res$BCD.mat #relative contributions of gains and losses to TBI

#4) adding results to dataframe and plotting beta div ~ pH and pH change
TBI.mat <- cbind(TBI.mat, TBI.res$TBI, TBI.res$BC, TBI.res$BCD.mat[,1:2]) %>% as.data.frame
TBI.mat %<>% rename('TBI' = `TBI.res$TBI`, 'loss' = B, 'gain' = C)
colnames(TBI.mat)[9:10] <- c('rel.loss','rel.gain')

cols <- brewer.pal(5,'Dark2')
par(mfrow=c(3,5),mar=c(4,4,0,0),cex=1)

plot(TBI ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[1],0.5),bty='l')
plot(loss ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[2],0.5),bty='l')
plot(gain ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[3],0.5),bty='l')
plot(rel.loss ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[4],0.5),bty='l')
plot(rel.gain ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[5],0.5),bty='l')
plot(TBI ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[1],0.5),bty='l')
plot(loss ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[2],0.5),bty='l')
plot(gain ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[3],0.5),bty='l')
plot(rel.loss ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[4],0.5),bty='l')
plot(rel.gain ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[5],0.5),bty='l')
plot(TBI ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[1],0.5),bty='l')
plot(loss ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[2],0.5),bty='l')
plot(gain ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[3],0.5),bty='l')
plot(rel.loss ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[4],0.5),bty='l')
plot(rel.gain ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[5],0.5),bty='l')

#5) repeating the whole thing with %diff (Bray-Curtis) instead of p/a
TBI.mat <- zoops %>% filter(year == 2012) %>% select(lake.name:pH) %>% as.data.frame
TBI.mat$ph.t1 <- zoops$pH[zoops$year == 1994]
TBI.mat$ph.diff <- zoops$pH[zoops$year == 2012] - zoops$pH[zoops$year == 1994]

### THIS IS THE ONLY LINE THAT CHANGES ###
TBI.res <- TBI(com.y1, com.y2, method="ruzicka", pa.tr=F, save.BC = T)
### THIS IS THE ONLY LINE THAT CHANGES ###

TBI.mat <- cbind(TBI.mat, TBI.res$TBI, TBI.res$BC, TBI.res$BCD.mat[,1:2]) %>% as.data.frame
TBI.mat %<>% rename('TBI' = `TBI.res$TBI`, 'loss' = B, 'gain' = C)
colnames(TBI.mat)[9:10] <- c('rel.loss','rel.gain')

plot(TBI ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[1],0.5),bty='l')
plot(loss ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[2],0.5),bty='l')
plot(gain ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[3],0.5),bty='l')
plot(rel.loss ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[4],0.5),bty='l')
plot(rel.gain ~ ph.t1, TBI.mat, xlab='initial pH', pch=17, col=alpha(cols[5],0.5),bty='l')
plot(TBI ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[1],0.5),bty='l')
plot(loss ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[2],0.5),bty='l')
plot(gain ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[3],0.5),bty='l')
plot(rel.loss ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[4],0.5),bty='l')
plot(rel.gain ~ pH, TBI.mat, xlab='final pH', pch=16, col=alpha(cols[5],0.5),bty='l')
plot(TBI ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[1],0.5),bty='l')
plot(loss ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[2],0.5),bty='l')
plot(gain ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[3],0.5),bty='l')
plot(rel.loss ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[4],0.5),bty='l')
plot(rel.gain ~ ph.diff, TBI.mat, pch=15, col=alpha(cols[5],0.5),bty='l')

par(mfrow=c(1,1))

#### Example III: Quebec hydropower reservoirs: ####
#### visualizing beta div & identifying important contributors

library(tibble)

#loading the dataset (emailed to you)
#PLEASE DELETE FILE AFTER WORKSHOP!!! UNPUBLISHED DATASET (BUT WILL BE ON FIGSHARE SOON)
LG_Y <- read.csv('/Users/vincentfugere/Desktop/LG_Y.csv', header=T)
str(LG_Y)

LG_Y_I<-LG_Y[LG_Y$I_R=="I",c(2,8:40)] # SUBSET OF IMPACTED STATIONS
LG_Y_R<-LG_Y[LG_Y$I_R=="R",c(2,8:40)] # SUBSET OF REFERENCE SITES

# CREATING RAREFIED RICHNESS METRICS FOR IMPACTED STATIONS
RRI<-specpool(LG_Y_I, smallsample=T, pool=LG_Y_I$TSI)
I_R <-c("I") # CREATING ONE VECTOR = IMPACTED 
RRI <-cbind(I_R, RRI) # MERGE THE RR METRICS AND THE VECTOR
rownames_to_column(RRI, var = "TSI")
RRI <- as_tibble(rownames_to_column(RRI, var = "TSI"))

# CREATING RAREFIED RICHNESS FOR REFERENCE SITES
RRR<-specpool(LG_Y_R, smallsample=T, pool=LG_Y_R$TSI)
I_R <-c("R")
RRR <-cbind(I_R, RRR)
rownames_to_column(RRR, var = "TSI")
RRR <- as_tibble(rownames_to_column(RRR, var = "TSI"))

RR_IR<-bind_rows(list(RRI, RRR)) # combine both matrices
RR_IR$TSI<-as.integer(RR_IR$TSI)

# LINEAR MODEL; RICHNESS IN RELATION TO TIME SINCE IMPOUNDMENT
RR_CBIR <- lm(jack2 ~ TSI*I_R, weight=RR_IR$n, data=RR_IR)
summary(RR_CBIR)

# GRAPH OF RICHNESS OVER TIME FOR IMPACTED AND REFERENCE SITES
ggplot(RR_IR, aes(TSI, jack2, color=I_R, size = n)) + geom_point() + 
  geom_smooth(method = "lm", se= T) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#BETA DIVERSITY

LG_com<-LG_Y[,-1:-7] #Keep only species columns

BDLG_dec<-beta.div.comp(LG_com)
BDLG_dec$part

BDLG<-beta.div(LG_com,method="hellinger")
?beta.div

#LCBD

LGLCBD = data.frame(BDLG$LCBD, BDLG$p.LCBD)
LGLCBD <-cbind(LGLCBD, LG_Y)
LGLCBD$pval <- ifelse(LGLCBD$BDLG.p.LCBD >= 0.05, c("non-sig"), c("sig")) 

LGLCBD$REFUD<- factor(LGLCBD$REFUD, levels = c("UR", "UL", "D", "R"))

ggplot(LGLCBD, aes(TSI, STATION, color=REFUD, size =BDLG.LCBD, shape=pval)) + geom_point() + 
  scale_shape_manual(values=c(1, 19)) +
  scale_color_manual(values=c("firebrick3", "darkorange", "blue", "chartreuse4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#SCBD

LGSCBD = data.frame(BDLG$SCBD)
LGSCBD <- as_tibble(rownames_to_column(LGSCBD, var = "SPECIES"))
LGSCBD_s <- LGSCBD[c(1:3,5,9,11,13,15,16,17,20,22),] # Select the most abundant species

# ROSE PLOT
ggplot(LGSCBD_s, aes(x=SPECIES, y=BDLG.SCBD)) + geom_bar(aes(fill=SPECIES), stat="identity") +
  scale_y_continuous(breaks = 0:10) +
  coord_polar() + labs(x = "", y = "")

# SUBSET for ROBERT-BOURASSA
RB_Y<-LG_Y[LG_Y$RES=="RB",-1:-7] 
BDRB<-beta.div(RB_Y,method="hellinger")
RBSCBD = data.frame(BDRB$SCBD)
RBSCBD <- as_tibble(rownames_to_column(RBSCBD, var = "SPECIES"))
RBSCBD_s <- RBSCBD[c(1:3,5,9,11,13,15,16,17,20,22),]
ggplot(RBSCBD_s, aes(x=SPECIES, y=BDRB.SCBD)) + geom_bar(aes(fill=SPECIES), stat="identity") +
  scale_y_continuous(breaks = 0:10) +
  coord_polar() + labs(x = "", y = "")

# SUBSET for OPINACA
OP_Y<-LG_Y[LG_Y$RES=="OP",-1:-7]
BDOP<-beta.div(OP_Y,method="hellinger")
OPSCBD = data.frame(BDOP$SCBD)
OPSCBD <- as_tibble(rownames_to_column(OPSCBD, var = "SPECIES"))
OPSCBD_s <- OPSCBD[c(1:3,5,9,11,13,15,16,17,20,22),]
ggplot(OPSCBD_s, aes(x=SPECIES, y=BDOP.SCBD)) + geom_bar(aes(fill=SPECIES), stat="identity") +
  scale_y_continuous(breaks = 0:10) +
  coord_polar() + labs(x = "", y = "")

# SUBSET for CANIAPISCAU
CA_Y<-LG_Y[LG_Y$RES=="CA",-1:-7]
BDCA<-beta.div(CA_Y,method="hellinger")
CASCBD = data.frame(BDCA$SCBD)
CASCBD <- as_tibble(rownames_to_column(CASCBD, var = "SPECIES"))
CASCBD_s <- CASCBD[c(1:3,5,9,11,13,15,16,17,20,22),]
ggplot(CASCBD_s, aes(x=SPECIES, y=BDCA.SCBD)) + geom_bar(aes(fill=SPECIES), stat="identity") +
  scale_y_continuous(breaks = 0:10) +
  coord_polar() + labs(x = "", y = "")
