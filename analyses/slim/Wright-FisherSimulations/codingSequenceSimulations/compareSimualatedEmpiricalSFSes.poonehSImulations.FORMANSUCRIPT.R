require(ggplot2)
require(dplyr)
require(reshape2)
# Lets start with one population then can figure out how to generalize

### solved a scary R thing: if you don't properly facet over all your vars, then ggplot can place the bar plots on top of each other, masking true values. Doesn't give a warning! Bad. A couple safeguards in place: lower the alpha value so you can see if there are multiple plots hiding, but also: giving an internalID to each unique grouping (all replicates within that grouping have the same internalID), and add group=internalID to your aes(). That way, if you forget to facet over something, it's okay, you'll get a bunch of plots next to each other that will make you notice. They won't stack #######
######## initial set up #######
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/SFSes/plots/"
#pops=c("CA","AK","AL","COM","KUR")
pops=c("AK","CA")
todaysdate=format(Sys.Date(),format="%Y%m%d")
genotypeDate=20181119
projection.date=20181221
hetFilter=0.75 # level that hets were filtered
dir.create(plot.dir,recursive = T)

############ functions ###################
# ####################### fold sfs function ###############
foldSFS <- function(sfs){
  foldedSFS <- data.frame()
  ss=length(sfs$frequency) - 1 # this is the ss in chromosomes
  foldedBin=ss/2  # half as many as ss ; do seq from 0 --> foldedBins to include monomorphics
  # note you have to start seq at 0 to include monomorphic bin
  for(i in seq(0,foldedBin)){
    # if it's not the center point (which doesn't get added together)
    # see wakeley coalescent eq 1.2
    if(i==ss-i){
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count)/2)) # if it's add midpoint (foldedLen, including 0 monomorphic bin), add together and divide by two (equivalent of not adding together)
    }
    else{
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count))) # if not at mid point, just add together like normal
    }
  }
  return(foldedSFS)
}

############## get empirical ###########
syn.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/cds/synonymous/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2/",sep="")
mis.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/cds/missense/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2/",sep="")
neut.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/neutral/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2-plusMonomorphic/",sep="")

all.empirical.sfs.folded.noMono <- data.frame()
for(pop in pops){
  
  ############## prepare the data ##############
  neut.input <- list.files(neut.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  syn.input <- list.files(syn.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  mis.input <- list.files(mis.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  neut.sfs <- melt(read.table(neut.input,skip = 1,header = T)) # skip first line: "13 folded "KUR""
  neut.sfs$mutLabel <- "neutral"
  syn.sfs <- melt(read.table(syn.input,skip = 1,header = T))
  syn.sfs$mutLabel <- "synonymous"
  mis.sfs <- melt(read.table(mis.input,skip = 1,header = T))
  mis.sfs$mutLabel <- "missense"
  # want to plot them separately and all together for the population
  # want to exclude d0_0 bin and 0-bins
  ########### put the syn, neut and mis sfses all together: ######
  
  all.sfs <- rbind.data.frame(syn.sfs,mis.sfs,neut.sfs)
  all.sfs$frequency <- sapply(strsplit(as.character(all.sfs$variable),"_"),"[",2)
  ######## exclude bins that are above the fold #########
  all.sfs_noMono_noZero <- all.sfs[all.sfs$value!=0 & all.sfs$frequency!="0",]
  ######## add proportional as well: ########
  all.sfs_noMono_noZero <- all.sfs_noMono_noZero %>%
    group_by(mutLabel) %>%
    mutate(proportion=value/sum(value))
  ####### order factors: #########
  all.sfs_noMono_noZero$mutLabel <- factor(all.sfs_noMono_noZero$mutLabel, levels=c("neutral","synonymous","missense"))
  ############ add identifiers to match simulated data #####
  all.sfs_noMono_noZero$population <- pop # population
  all.sfs_noMono_noZero$category <- "empirical"
  all.sfs_noMono_noZero$state <- "empirical"
  all.sfs_noMono_noZero$count <-   all.sfs_noMono_noZero$value
  
  ############ add to overall empirical set ########
  all.empirical.sfs.folded.noMono <- rbind(all.empirical.sfs.folded.noMono,data.frame(all.sfs_noMono_noZero))
  ########### plot count-sfses all together ########
  # sfsPlot <- ggplot(all.sfs_noMono_noZero,aes(x=as.numeric(frequency),y=proportion,fill=label))+
  #   geom_bar(stat="identity",position = "dodge")+
  #   theme_bw()+
  #   theme(legend.title=element_blank(),legend.background = element_rect(fill = "transparent"),legend.position = c(0.5,0.9),legend.direction = "horizontal",legend.spacing.x = unit(.15,"cm"),text=element_text(size=18))+
  #   xlab("frequency")+
  #   ggtitle(paste(pop," folded SFS (projected)\nhet filter: ",hetFilter,sep=""))+
  #   scale_x_continuous(breaks=c(seq(1,max(as.numeric(all.sfs_noMono_noZero$frequency)))))+
  #   scale_fill_manual(values=c("gray","dodgerblue","darkred"))
  # sfsPlot
  # ggsave(paste(plot.dir,pop,".cds.neutral.projected.prop.SFS.hetFilter.",hetFilter,".pdf",sep=""),sfsPlot,device = "pdf",height=5,width=7)
  # 
  ########## Plot each individually ############
  sfsPlot2 <- ggplot(all.sfs_noMono_noZero,aes(x=as.numeric(frequency),y=value,fill=mutLabel))+
    geom_bar(stat="identity",position = "dodge")+
    theme_bw()+
    #theme(legend.title=element_blank(),legend.background = element_rect(fill = "transparent"),legend.position = "none",legend.direction = "horizontal",legend.spacing.x = unit(.15,"cm"),text=element_text(size=18))+
    xlab("frequency")+
    ylab("count")+
    ggtitle(paste(pop," folded SFS (projected)\nhet filter: ",hetFilter,sep=""))+
    scale_x_continuous(breaks=c(seq(1,max(as.numeric(all.sfs_noMono_noZero$frequency)))))+
    scale_fill_manual(values=c("gray","dodgerblue","darkred")) #+
    #facet_wrap(~state)
  sfsPlot2
  ggsave(paste(plot.dir,pop,".cds.neutral.projected.counts.faceted.SFS.hetFilter.",hetFilter,".pdf",sep=""),sfsPlot2,device = "pdf",height=5,width=12)

}
####### empirical sfses are now in: all.empirical.sfs.folded.noMono are folded and have no monomorphic counts ###########
################################## get simulated ###############################
sim.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/SFSes/"
#rundate=20190424 # set of simulations you're interested in (if is across different rundates you can list popsModelsRundates explicitly)
hs=c("0.5","0") # set of hs you're interested in
#popMods=c("AL/1D.2Epoch.1.5Mb.cds", "AK/1D.2Epoch.1.5Mb.cds") # population and corresponding models you're interested in
######## update with Ponoeh's new simulations (CA only)
#popModsDate=c("AK/1D.3Epoch.LongerRecovery/20191202/","CA/1D.3Epoch.LongerRecovery/20191013/") 
popModsDate="CA/1D.3Epoch.35GenBottleneck.LongerRecovery.newUseThis/20200310/"
popmodels=list()
for(i in popModsDate){
  for(h in hs){
    pm=paste(i,"h_",h,"/",sep="")
    popmodels=c(popmodels,pm)
  }
}

#numreps=1 # number of reps 
# go through pre and post 
# empty df for count sfses
all.sim.sfs.folded.noMono <- data.frame()
#### add a unique identifier for each one to keep track of total 
internalID=0 # counting each distinct category (within each internal ID will be many replicates)
for(pm in popmodels){ # go through each population-model dir
  for(state in c("Post","Pre")){  # do separtely for pre and post contraction 
    #for(rep in seq(1,numreps)){ # go through each replicate
    for(mutType in c(1,2)){ # go through each mutation type
      # mutation type 1 
      # for all reps
      sim.files <- list.files(paste(sim.dir,pm,sep=""),pattern=paste("*rep.*.",state,"Contraction.slim.mutType",as.character(mutType),".unfolded.sfs.R.format",sep=""),full.names = T) # this is for one rep, but could potentially do all reps and average here by turning rep into *
      # Fill in info: 
      internalID = internalID+ 1
      for(sim.file in sim.files){
        ################## mutation type 1 ####################
        sfs <- read.table(sim.file,header=T)
        sfs.folded <- foldSFS(sfs) #### the input SFS from slim is UNFOLDED, fold here ***
        # get rid of monomorphic sites
        sfs.folded.noMono <- sfs.folded[sfs.folded$frequency!=0,]
        sfs.folded.noMono$mutType <- mutType
        sfs.folded.noMono$popModel <- as.character(pm)
        sfs.folded.noMono$population <- unlist(lapply(strsplit(sfs.folded.noMono$popModel,"/"),"[",1))
        sfs.folded.noMono$model <- unlist(lapply(strsplit(sfs.folded.noMono$popModel,"/"),"[",2))
        sfs.folded.noMono$modelDate <- unlist(lapply(strsplit(sfs.folded.noMono$popModel,"/"),"[",3))
        sfs.folded.noMono$dominance_h <- unlist(lapply(strsplit(sfs.folded.noMono$popModel,"/"),"[",4))
        sfs.folded.noMono$state <- paste("simulated: ",state,"-Contraction",sep="")
        sfs.folded.noMono$category <- "simulated"
        # get the replicate number: oh this doesn't work when it's two digits. 
        sfs.folded.noMono$replicate <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(sim.file,".replicate_"),"[",2)),"\\."),"[",1)))
        sfs.folded.noMono$internalID <- internalID # to make sure things don't get combined falsely
        # also get proportions: 
        sfs.folded.noMono <- sfs.folded.noMono %>%
          mutate(proportion=count/sum(count))
        
        all.sim.sfs.folded.noMono <- rbind(all.sim.sfs.folded.noMono, sfs.folded.noMono)
      }
    }
  }
}

####### simulated sfses are now in: all.sim.sfs.folded.noMono and are folded and have no monomorphic counts ###########
## Based on my slim script, mutation type 1 is "synonymous" and mutation type 2 is "missense" (not this is different from Chris' script labels)
# so assign those labels
all.sim.sfs.folded.noMono$mutLabel <- NA
all.sim.sfs.folded.noMono[all.sim.sfs.folded.noMono$mutType=="1",]$mutLabel <- "synonymous"
all.sim.sfs.folded.noMono[all.sim.sfs.folded.noMono$mutType=="2",]$mutLabel <- "missense"
# skipping neutral for now, but should find those simulated neutral sfses

################## get average of each bin of SFS #####################
require(dplyr)
# count simulations per category -- should be 20-21 or so per category. AK has 21 in some, but that's okay because each is an indepedent replicate
# if there are 40 then somethings gone weird with numbering 
ns <- all.sim.sfs.folded.noMono %>%
  group_by(popModel,model,modelDate,mutLabel,category,state,dominance_h,population,frequency,internalID) %>%
  summarise(numSims=n())
means <- all.sim.sfs.folded.noMono %>%
  group_by(popModel,model,modelDate,mutLabel,category,state,dominance_h,population,frequency,internalID) %>%
  summarise_at(c("count", "proportion"), list(MEAN=mean),na.rm=TRUE)
sds <- all.sim.sfs.folded.noMono %>%
  group_by(popModel,model,modelDate,mutLabel,category,state,dominance_h,population,frequency,internalID) %>%
  summarise_at(c("count", "proportion"), list(SD=sd),na.rm=TRUE)
# join ns, means and sds together: 
summarizedSFSes <-  left_join(ns,means) %>%
  left_join(., sds) 
# Joining, by = c("popModel", "modelDate", "mutLabel", "category", "state", "dominance_h", "population", "frequency")
summarizedSFSes$count_SE <- summarizedSFSes$count_SD / sqrt(summarizedSFSes$numSims)
summarizedSFSes$proportion_SE <- summarizedSFSes$proportion_SD / sqrt(summarizedSFSes$numSims)


########### Plot Average with +-1 SD ###############
################# BE CAREFUL; IF YOU FORGOT TO GROUP BY ONE YOUR VARIABLES THE PLOTS WILL ECLIPES EACH OTHER OR ADD TOGETHER ODDLY ############
######### choose the model you want so you don't have to factor over it:
desiredModelToPlot="1D.3Epoch.35GenBottleneck.LongerRecovery.newUseThis"
#summarizedSFSes$state <- factor(summarizedSFSes$state,levels=c("simulated: Pre-Contraction","simulated: Post-Contraction"))
summarizedSFSesToPlot <- summarizedSFSes[summarizedSFSes$model==desiredModelToPlot,]
summarizedSFSesToPlot$state <- factor(summarizedSFSesToPlot$state,levels=c("simulated: Pre-Contraction","simulated: Post-Contraction"))
countMeanPlot1 <- ggplot(summarizedSFSesToPlot,aes(x=frequency,y=count_MEAN,group=internalID,fill=state))+
  geom_bar(stat="identity",position="dodge",alpha=0.75)+
  geom_errorbar(aes(ymin=count_MEAN-count_SD, ymax=count_MEAN+count_SD), width=.2,
                position=position_dodge(0.9)) +
  scale_fill_manual(values=c("blue","darkred"))+
  #facet_grid(~population~model~dominance_h~mutLabel)+
  facet_grid(~population~dominance_h~mutLabel)+
  theme_bw()+
  scale_x_continuous(breaks=c(seq(1,max(as.numeric(summarizedSFSes$frequency)))))+
  ggtitle("Simulated SFS averaged over simluation replicates\nError Bars denote one standard deviation")
countMeanPlot1
#ggsave(paste(plot.dir,"AK.CA.SimulatedSFSes.3EpochLongerContraction.forMS.",todaysdate,".pdf",sep=""),countMeanPlot1,device="pdf",width=7,height=5)


################## combine empirical and simulated cds #################
# things in common that I need to plot are:
# "mutLabel"   "frequency"  "proportion" "pop"        "category"   "state"      "count" "popModel"
# make sure there aren't multiple models or dates per population here though
all.empirical.sfs.folded.noMono$popModel <- "empirical"
all.empirical.sfs.folded.noMono$model <- "empirical"
all.empirical.sfs.folded.noMono$internalID <- 0
all.empirical.sfs.folded.noMono$popLabel2 <- all.empirical.sfs.folded.noMono$population
### create some dummy variables:
all.empirical.sfs.folded.noMono[,c("count_SD","proportion_SD","count_SE","proportion_SE","dominance_h")] <- NA
## add columns to sims:
summarizedSFSes$count <- summarizedSFSes$count_MEAN
summarizedSFSes$proportion <- summarizedSFSes$proportion_MEAN
#### add a second population label to encompass the two models:
summarizedSFSes$popLabel2 <- summarizedSFSes$population
# empirical is just the same: 

########## need to rename columns to get them all to match ###########

colsIWant <- intersect(names(all.empirical.sfs.folded.noMono),names(summarizedSFSes))
# and want to exclude some pops:
#popsIWant <- c("AK","AL") # skip COM for now -- it's weird
#popsIWant <- c("AK","genericPop.LongerContract")# renamed genericPop to AK
EMPpopsIWant <- c("CA")
SIMpopsIWant <- c("CA")
SIMmodelsIWant <- "1D.3Epoch.35GenBottleneck.LongerRecovery.newUseThis"
# okay since AK and longer contraction don't make a big diff (dadi MLE official (AK) vs. elsewhere on the MLE curve from the grid search, going to use the genericPop.LongerContract only)
# SEPARATE by h = 0 , h = 0.5:
# h = 0
combo.cds.h_0 <- rbind(all.empirical.sfs.folded.noMono[all.empirical.sfs.folded.noMono$population %in% EMPpopsIWant,colsIWant],summarizedSFSes[summarizedSFSes$population %in% SIMpopsIWant & summarizedSFSes$dominance_h=="h_0" & summarizedSFSes$model %in% SIMmodelsIWant,colsIWant])
# h = 0.5
combo.cds.h_0.5 <- rbind(all.empirical.sfs.folded.noMono[all.empirical.sfs.folded.noMono$population %in% EMPpopsIWant,colsIWant],summarizedSFSes[summarizedSFSes$population %in% SIMpopsIWant & summarizedSFSes$dominance_h=="h_0.5"  & summarizedSFSes$model %in% SIMmodelsIWant,colsIWant])
# exclude neutral for now --> 
# order factors
####### order factors: #########
combo.cds.h_0$state <- factor(combo.cds.h_0$state, levels=c("empirical","simulated: Pre-Contraction","simulated: Post-Contraction"))
combo.cds.h_0.5$state <- factor(combo.cds.h_0.5$state, levels=c("empirical","simulated: Pre-Contraction","simulated: Post-Contraction"))
########## plot recessive h=0 only with empirical ##########
#### just using AK genericPop.LongerContraction instead of "AK" dadi MLE
# yield same SFS 
# using popModel to group which is a combo of pop,date,model and dom coefficient -- so same as internal ID, and interacting with state (pre or post contraction)
p1a <- ggplot(combo.cds.h_0[combo.cds.h_0$mutLabel!="neutral",],aes(x=as.numeric(frequency),y=as.numeric(proportion),fill=state,group=interaction(popModel,state)))+
  geom_bar(position='dodge',stat="identity",alpha=0.75)+
  facet_wrap(as.factor(popLabel2)~mutLabel,scales="free_x")+
  geom_errorbar(aes(ymin=proportion-proportion_SE, ymax=proportion+proportion_SE), width=.2,
                position=position_dodge(0.9)) +
  theme_bw()+
  scale_fill_manual(values=c("darkgray","blue","darkred"))+
  xlab("Frequency")+
  ylab("Proportion")+
  scale_x_continuous(breaks=c(seq(1,max(as.numeric(combo.cds.h_0$frequency)))))+
  ggtitle(paste("Dominance Coefficient (h): 0\nError Bars denote one standard error",sep=""))+
  theme(legend.position="right",legend.background = element_rect("transparent"),legend.title=element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size=14))
p1a
# you'll get this warning:
#Warning message:
#  Removed 14 rows containing missing values (geom_errorbar). 
# because empirical doesn't have error bars' == that's the expected behavior.

ggsave(paste(plot.dir,"ComparingSimulatedEmpiricalCDS.SFSes.h_0.",todaysdate,".pdf",sep=""),p1a,device="pdf",width=11,height=6)

######## plot additive  #############
p1b <- ggplot(combo.cds.h_0.5[combo.cds.h_0.5$mutLabel!="neutral",],aes(x=as.numeric(frequency),y=as.numeric(proportion),fill=state,group=interaction(popModel,state)))+
  geom_bar(position='dodge',stat="identity",alpha=0.75)+
  facet_wrap(as.factor(popLabel2)~mutLabel,scales="free_x")+
  geom_errorbar(aes(ymin=proportion-proportion_SD, ymax=proportion+proportion_SD), width=.2,
                position=position_dodge(0.9)) +
  theme_bw()+
  scale_fill_manual(values=c("darkgray","blue","darkred"))+
  xlab("Frequency")+
  ylab("Proportion")+
  scale_x_continuous(breaks=c(seq(1,max(as.numeric(combo.cds.h_0$frequency)))))+
  ggtitle(paste("Dominance Coefficient (h): 0.5\nError Bars denote one standard deviation",sep=""))+
  theme(legend.position="right",legend.background = element_rect("transparent"),legend.title=element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size=14))
p1b
# you'll get this warning:
#Warning message:
#  Removed 14 rows containing missing values (geom_errorbar). 
# because empirical doesn't have error bars' == that's the expected behavior.
ggsave(paste(plot.dir,"ComparingSimulatedEmpiricalCDS.SFSes.h_0.5.",todaysdate,".pdf",sep=""),p1b,device="pdf",width=11,height=6)
###################### PLOT for manuscript, just plotting CA and plotting separately so you can stack them and just show one synonymous #######
### Kirk request: show empirical next to post contraction
combo.cds.h_0$state <- factor(combo.cds.h_0$state, levels=c("simulated: Pre-Contraction","simulated: Post-Contraction","empirical"))
combo.cds.h_0.5$state <- factor(combo.cds.h_0$state, levels=c("simulated: Pre-Contraction","simulated: Post-Contraction","empirical"))

pops=c("CA")
for(pop in pops){
  for(mutLabel in c("synonymous","missense")){
    ##### recessive: #####
    p2a <- ggplot(combo.cds.h_0[combo.cds.h_0$mutLabel==mutLabel & combo.cds.h_0$population==pop,],aes(x=as.numeric(frequency),y=as.numeric(proportion),fill=state,group=interaction(popModel,state)))+
      geom_bar(position='dodge',stat="identity",alpha=0.75)+
      facet_wrap(~mutLabel,scales="free_x")+
      geom_errorbar(aes(ymin=proportion-proportion_SE, ymax=proportion+proportion_SE), width=.2,
                    position=position_dodge(0.9)) +
      theme_bw()+
      scale_fill_manual(values=c("blue","darkred","darkgray"))+
      xlab("Frequency")+
      ylab("Proportion")+
      scale_x_continuous(breaks=c(seq(1,max(as.numeric(combo.cds.h_0$frequency)))))+
      ggtitle("simulated mutations are recessive")+
      #ggtitle(paste("Dominance Coefficient (h): 0\nError Bars denote one standard error",sep=""))+
      theme(legend.position="none",legend.background = element_rect("transparent"),legend.title=element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size=14))
    p2a
    # you'll get this warning:
    #Warning message:
    #  Removed 14 rows containing missing values (geom_errorbar). 
    # because empirical doesn't have error bars' == that's the expected behavior.
    
    ggsave(paste(plot.dir,pop,".ComparingSimulatedEmpiricalCDS.35GenContraction.SFSes.h_0.",mutLabel,".",todaysdate,".pdf",sep=""),p2a,device="pdf",width=3.5,height=2.5)
    
    ######## plot additive  #############
    p2b <- ggplot(combo.cds.h_0.5[combo.cds.h_0.5$mutLabel==mutLabel & combo.cds.h_0.5$population==pop,],aes(x=as.numeric(frequency),y=as.numeric(proportion),fill=state,group=interaction(popModel,state)))+
      geom_bar(position='dodge',stat="identity",alpha=0.75)+
      facet_wrap(~mutLabel,scales="free_x")+
      geom_errorbar(aes(ymin=proportion-proportion_SE, ymax=proportion+proportion_SE), width=.2,
                    position=position_dodge(0.9)) +
      theme_bw()+
      scale_fill_manual(values=c("blue","darkred","darkgray"))+
      xlab("Frequency")+
      ylab("Proportion")+
      scale_x_continuous(breaks=c(seq(1,max(as.numeric(combo.cds.h_0$frequency)))))+
      ggtitle("simulated mutations are additive")+
      #ggtitle(paste("Dominance Coefficient (h): 0\nError Bars denote one standard error",sep=""))+
      theme(legend.position="none",legend.background = element_rect("transparent"),legend.title=element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size=14))
    p2b
    # you'll get this warning:
    #Warning message:
    #  Removed 14 rows containing missing values (geom_errorbar). 
    # because empirical doesn't have error bars' == that's the expected behavior.
    ggsave(paste(plot.dir,pop,".ComparingSimulatedEmpiricalCDS.35GenContraction.SFSes.h_0.5.",mutLabel,".",todaysdate,".pdf",sep=""),p2b,device="pdf",width=3.5,height=2.5)
  }}



############################ so clearly post-contraction fits better; focus just on that + missense with both h's ##########################




# don't need: 
additiveRecessiveMissense <- rbind(combo.cds.h_0.5[!combo.cds.h_0.5$mutLabel %in% c("synonymous","neutral") & (combo.cds.h_0.5$state %in% c("empirical","simulated: Post-Contraction")),],combo.cds.h_0[!combo.cds.h_0$mutLabel %in% c("synonymous","neutral") & (combo.cds.h_0$state %in% c("empirical","simulated: Post-Contraction")),])
# label empirical h_0 as unknown
additiveRecessiveMissense[additiveRecessiveMissense$category=="empirical",]$dominance_h <- "empirical"
p1c <- ggplot(additiveRecessiveMissense,aes(x=as.numeric(frequency),y=as.numeric(proportion),fill=dominance_h))+
  geom_bar(position='dodge',stat="identity",alpha=0.75)+
  facet_wrap(as.factor(popLabel2)~mutLabel,scales="free_x")+
  geom_errorbar(aes(ymin=proportion-proportion_SE, ymax=proportion+proportion_SE), width=.2,
                position=position_dodge(0.9)) +
  theme_bw()+
  scale_fill_manual(values=c("darkgray","firebrick3","orange4"))+
  xlab("Frequency")+
  ylab("Proportion")+
  scale_x_continuous(breaks=c(seq(1,max(as.numeric(combo.cds.h_0$frequency)))))+
  ggtitle(paste("Dominance Coefficient (h): 0.5\nError Bars denote one standard error",sep=""))+
  theme(legend.position="right",legend.background = element_rect("transparent"),legend.title=element_blank(),legend.key.size = unit(1,"cm"),legend.text = element_text(size=14))
p1c
# you'll get this warning:
#Warning message:
#  Removed 14 rows containing missing values (geom_errorbar). 
# because empirical doesn't have error bars' == that's the expected behavior.
ggsave(paste(plot.dir,"ComparingFitOfAdditive.Recessive.MissenseToEmpirical.",todaysdate,".pdf",sep=""),p1c,device="pdf",width=11,height=5)
#################### make sure this isn't aggregating over something incorrectly; so far I think t's okay ########
# plot ratio of NS and S 

# NS_S_Ratio.h_0.5 <- combo.cds.h_0.5 %>%
#   group_by(state,population,mutLabel,model) %>%
#   tally(sum(count)) %>%
#   summarise(Ratio = n[mutLabel == "missense" ] / n[mutLabel == "synonymous"])
# 
# NS_S_Ratio.h_0 <- combo.cds.h_0 %>%
#   group_by(state,population,mutLabel) %>%
#   tally(sum(count)) %>%
#   summarise(Ratio = n[mutLabel == "missense"] / n[mutLabel == "synonymous"])
# 
# NS_S_Ratio.h_0.5$dominance_h <- "h_0.5"
# NS_S_Ratio.h_0$dominance_h <- "h_0"
# 
# NS_S_Ratio = rbind(NS_S_Ratio.h_0.5,NS_S_Ratio.h_0)
# p2 <- ggplot(NS_S_Ratio,aes(x=population,y=Ratio,fill=state))+
#   geom_bar(stat="identity", position="dodge")+
#   ylab("Ratio of Missense / Synonymous SNPs")+
#   theme_bw()+
#   facet_wrap(~dominance_h)+
#   scale_fill_manual(values=c("darkgray","blue","darkred"))+
#   ggtitle(paste("Simulation info: ", unique(combo.cds.h_0.5[combo.cds.h_0.5$category!="empirical",]$model," ",sep="")))
# p2
# ggsave(paste(plot.dir,"Ratio_NS_S.",todaysdate,".pdf",sep=""),p2,device="pdf",width=11,height=5)
