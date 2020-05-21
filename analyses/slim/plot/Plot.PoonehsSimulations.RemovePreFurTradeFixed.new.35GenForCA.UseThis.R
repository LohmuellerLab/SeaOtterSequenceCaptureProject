################### Plotting Pooneh's simulations ################
# this is a new version of the script 20200312 -- going to plot things for manuscript
# the CA 3-Epoch and AK 5-Epoch and the 2D models
# now CA has 35 gen bneck just like AK
# and have burn-in fixed sites removed (but fixed sites during bneck are kept)
# 
# want to plot generations after bottleneck
# want to facet by additive/recessive and by population 
require(ggplot2)
require(dplyr)
require(RColorBrewer)
require(viridis) # nice colors
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
#migrationColors=list(NoGeneFlow="grey22","1perGen"="cyan4","5perGen"="#8C96C6","10perGen"="#88419D","25for2Gen"="#6E016B","25perGen"="darkgoldenrod4")
migrationColors=list(NoGeneFlow="#000000","1perGen"="#0072B2","5perGen"="#8C96C6","10perGen"="#88419D","25for2Gen"="#E69F00","25perGen"="#CC79A7")
textsize=14
################ I. California 3 Epoch (plot with and without recovery) ##################
# new CA sims with 35 gen bneck and removing fixed vars:
CA_3Epoch <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/concattedSummaries/RemovingFixedVariants_useThis/CA/1D.3Epoch.35Generations.LongerRecovery/2020311/20200311_CA_LoadPerGeneration.ThroughTime.AllReps.RemovedBurninFixedVar.txt",sep="\t",header=T)
# note in these new simulations, CA bneck lasts 35 generations
CA_3Epoch_noRec <- CA_3Epoch[CA_3Epoch$generation<36,]
# this now has no recovery periods
# rest should work as before (check it)
CA_3Epoch_noRec$label <- ""
CA_3Epoch_noRec[CA_3Epoch_noRec$h==0,]$label <- "recessive"
CA_3Epoch_noRec[CA_3Epoch_noRec$h==0.5,]$label <- "additive"
# gather means:
CA_3Epoch_noRec_MeansOverSimReplicates = CA_3Epoch_noRec %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanLoad=mean(L))

######################### p1: Plot the bottleneck load only ( no recovery ) ################
for(dominance in c("recessive","additive")){
  p1 <- ggplot(CA_3Epoch_noRec_MeansOverSimReplicates[CA_3Epoch_noRec_MeansOverSimReplicates$label==dominance,],aes(x=generation,y=meanLoad,color=population))+
    geom_line(data=CA_3Epoch_noRec[CA_3Epoch_noRec$label==dominance,],aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since contraction")+
    ylab("mean genetic load")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
  p1
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/NewPlots_RemoveFixedVars_ForManuscript/CA.geneticLoad.ContractionOnly.RemovedBurnInFixedVars.",dominance,".pdf",sep=""),p1,height=2,width=3)
}
######## plot with dominance faceted  #########
CA_3Epoch_noRec_MeansOverSimReplicates$label2 <- paste(CA_3Epoch_noRec_MeansOverSimReplicates$label," mutations",sep="")
CA_3Epoch_noRec_MeansOverSimReplicates$label2 <- factor(CA_3Epoch_noRec_MeansOverSimReplicates$label2, levels=c("recessive mutations","additive mutations"))
CA_3Epoch_noRec$label2 <- paste(CA_3Epoch_noRec$label," mutations",sep="")
CA_3Epoch_noRec$label2 <- factor(CA_3Epoch_noRec$label2, levels=c("recessive mutations","additive mutations"))

p1c <- ggplot(CA_3Epoch_noRec_MeansOverSimReplicates,aes(x=generation,y=meanLoad,color=population))+
  geom_line(data=CA_3Epoch_noRec,aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.5)+
  geom_line(size=1.5)+
  theme_bw()+
  facet_wrap(~label2,scales="free_y",ncol=1)+
  theme(text = element_text(size=12),legend.position = "none",strip.background = element_rect("transparent"))+
  xlab("generations since initial contraction")+
  ylab("mean genetic load")+
  scale_color_manual(values=unlist(colors))+
  theme(text = element_text(size=11),strip.background = element_blank(),legend.title = element_blank(),strip.text = element_text(hjust=0,size=textsize))
p1c
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/NewPlots_RemoveFixedVars_ForManuscript/CA.geneticLoad.ContractionOnly.RemovedBurnInFixedVars.FacetedDominance.pdf",sep=""),p1c,height=4,width=3)

########## II. Contraction + Recovery (same model but going beyond 36 gen) #################
CA_3Epoch$label <- ""
CA_3Epoch[CA_3Epoch$h==0,]$label <- "recessive"
CA_3Epoch[CA_3Epoch$h==0.5,]$label <- "additive"
# gather means:
CA_3Epoch_meansOverSimReplicates = CA_3Epoch %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanLoad=mean(L))

######################### p3: Plot with recovery load ################
hLines=data.frame(population=c("CA"),intercept=c(36))
for(dominance in c("recessive","additive")){
  
  p3 <- ggplot(CA_3Epoch_meansOverSimReplicates[CA_3Epoch_meansOverSimReplicates$label==dominance,],aes(x=generation,y=meanLoad,color=population))+
    geom_line(data=CA_3Epoch[CA_3Epoch$label==dominance,],aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since contraction")+
    ylab("mean genetic load")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
    
  p3
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/NewPlots_RemoveFixedVars_ForManuscript/CA.geneticLoad.Recovery.RemovedBurnInFixedVars.",dominance,".pdf",sep=""),p3,height=4,width=7)
  
}
##### plot faceted by dominance: #####
CA_3Epoch_meansOverSimReplicates$label2 <- paste(CA_3Epoch_meansOverSimReplicates$label," mutations",sep="")
CA_3Epoch_meansOverSimReplicates$label2 <- factor(CA_3Epoch_meansOverSimReplicates$label2, levels=c("recessive mutations","additive mutations"))
CA_3Epoch$label2 <- paste(CA_3Epoch$label," mutations",sep="")
CA_3Epoch$label2 <- factor(CA_3Epoch$label2, levels=c("recessive mutations","additive mutations"))
p3b <- ggplot(CA_3Epoch_meansOverSimReplicates,aes(x=generation,y=meanLoad,color=population))+
  geom_line(data=CA_3Epoch,aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.5)+
  geom_line(size=1.5)+
  geom_vline(data=hLines,aes(xintercept = intercept),linetype="dashed",color="black")+
  theme_bw()+
  facet_wrap(~label2,scales="free_y",ncol=1)+
  theme(text = element_text(size=12),legend.position = "none",strip.background = element_rect("transparent"))+
  xlab("generations since initial contraction")+
  ylab("mean genetic load")+
  scale_color_manual(values=unlist(colors))+
  theme(text = element_text(size=11),strip.background = element_blank(),legend.title = element_blank(),strip.text = element_text(hjust=0,size=textsize))
p3b
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/NewPlots_RemoveFixedVars_ForManuscript/CA.geneticLoad.Recovery.RemovedBurnInFixedVars.FacetedDominance.pdf",sep=""),p3b,height=4,width=3)

################### III. AK 5 epoch double bottleneck ##########
hLines2=data.frame(population=c("AK"),intercept=c(36,36+14,36+14+6))

doubleBottleneckAK <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/concattedSummaries/RemovingFixedVariants_useThis/AK/1D.5Epoch/20200311/20200313_LoadPerGeneration.ThroughTime.AllReps.RemovedBurninFixedVar.txt",header=T)
doubleBottleneckAK$label <- NA
doubleBottleneckAK[doubleBottleneckAK$h==0,]$label <- "recessive mutations"
doubleBottleneckAK[doubleBottleneckAK$h==0.5,]$label <- "additive mutations"
doubleBottleneckAK$label <- factor(doubleBottleneckAK$label,levels=c("recessive mutations","additive mutations"))
# gather means:
doubleBottleneckAK_means = doubleBottleneckAK %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanFitness=mean(W),meanLoad=mean(L))

p4 <- ggplot(doubleBottleneckAK_means,aes(x=generation,y=meanLoad,color=population))+
  geom_line(data=doubleBottleneckAK,aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.5)+
  geom_line(size=1.5)+
  geom_vline(data=hLines2,aes(xintercept = intercept),linetype="dashed",color="black")+
  theme_bw()+
  facet_wrap(~label,scales="free_y",ncol=1)+
  theme(text = element_text(size=12),legend.position = "none",strip.background = element_rect("transparent"))+
  xlab("generations since contraction")+
  ylab("mean genetic load")+
  scale_color_manual(values=unlist(colors))+
  theme(text = element_text(size=11),strip.background = element_blank(),legend.title = element_blank(),strip.text = element_text(hjust=0,size=textsize))
p4
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/NewPlots_RemoveFixedVars_ForManuscript/AK.geneticLoad.5Epoch.RemovedBurnInFixedVars.FacetedDominance.pdf",sep=""),p4,height=4,width=3)
################# IV. Model with Migration ###########
migModels <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/concattedSummaries/RemovingFixedVariants_useThis/CA_AK/2D.3Epoch.VaryingLevelsTranslocation/20200310/20200310_CA_AK_LoadPerGeneration.ThroughTime.AllReps.RemovedBurninFixedVar.txt",header=T)
bneckGen=4002 # after pops separated for 4K 
migModels_recent=migModels[(migModels$subpop==1 & migModels$generation>=bneckGen)|(migModels$subpop==2 & migModels$generation>=bneckGen),]
migModels_recent$genRescaled <- 0
migModels_recent[migModels_recent$subpop==2,]$genRescaled <- migModels_recent[migModels_recent$subpop==2,]$generation-bneckGen
migModels_recent[migModels_recent$subpop==1,]$genRescaled <- migModels_recent[migModels_recent$subpop==1,]$generation-bneckGen
##### get averages:
migModels_recent$label <- NA
migModels_recent[migModels_recent$h==0,]$label <- "recessive"
migModels_recent[migModels_recent$h==0.5,]$label <- "additive"
migModels_recent$pop <- NA
migModels_recent[migModels_recent$subpop==1,]$pop <- "AK"
migModels_recent[migModels_recent$subpop==2,]$pop <- "CA"
### nicer model names:
migModels_recent$modelID <- unlist(lapply(strsplit(as.character(migModels_recent$model),"\\."),tail,n=1)) # choose last element -- is element 4 for most models, but 3 for 'notranslocation'
# rename "notranslocation" to no migration
migModels_recent[migModels_recent$modelID=="NoTranslocation",]$modelID <- "NoGeneFlow"
#order factors:
migModels_recent$modelID <- factor(migModels_recent$modelID,levels=c("NoGeneFlow","1perGen","5perGen","10perGen","25for2Gen","25perGen"))
## the no translocation ones are bad 
# gather means:
migModels_recent_means = migModels_recent %>%
  group_by(model,pop,label,generation,modelID,genRescaled,subpop)  %>%
  summarise(meanFitness=mean(W),meanLoad=mean(L))
## when does migration start for each population?
hLines3=data.frame(pop=c("AK","AK","AK","CA","CA","CA"),intercept=c(36,54))
dominance="recessive"
  p5a <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanLoad,color=as.factor(modelID)))+
    #geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=L_mutationLoad,color=as.factor(modelID)),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines3,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~pop,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
    xlab("generations since contraction")+
    ylab("mean genetic load")+
    ggtitle(paste(dominance," mutations",sep=""))+
    scale_color_manual(values=unlist(migrationColors))
  
  p5a
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/NewPlots_RemoveFixedVars_ForManuscript/MIGRATION.",dominance,".no5-10.geneticLoad.AK.CA.NoFineLines.pdf",sep=""),p5a,height=3,width=8)

  
  dominance="additive"
  p5b <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanLoad,color=as.factor(modelID)))+
    #geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(replicate,modelID),x=genRescaled,y=L,color=as.factor(modelID)),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines3,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~pop,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
    scale_y_continuous(limits=c(.15,.275))+
    xlab("generations since contraction")+
    ylab("mean genetic load")+
    ggtitle(paste(dominance," mutations",sep=""))+
    scale_color_manual(values=unlist(migrationColors))
  
  p5b
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/NewPlots_RemoveFixedVars_ForManuscript/MIGRATION.",dominance,".no5-10.geneticLoad.AK.CA.NoFineLines.pdf",sep=""),p5b,height=3,width=8)
  
ggplot(migModels_recent,aes(x=generation,y=popsizeDIP))+
  geom_point()+
  facet_wrap(pop~h)
############ Plot faceted by dominance: THIS IS UGLY; just combine the top two plots which are plotted by dominance seprately in ppt ########
migModels_recent_means$label2 <- paste(migModels_recent_means$label," mutations",sep="")
migModels_recent_means$label2 <- factor(migModels_recent_means$label2, levels=c("recessive mutations","additive mutations"))
p5d <- ggplot(migModels_recent_means[!(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanLoad,color=as.factor(modelID)))+
  geom_line(size=1)+
  theme_bw()+
  #geom_vline(data=hLines3[hLines3$pop=="AK",],aes(xintercept = intercept),linetype="dashed",color="black")+
  facet_grid(~label2~pop)+
  theme(text = element_text(size=11),strip.background = element_blank(),legend.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(0.9)),strip.text = element_text(hjust=0,size=textsize))+
  xlab("generations since contraction")+
  ylab("mean genetic load")+
  scale_color_manual(values=unlist(migrationColors))+
  guides(colour=guide_legend(ncol=2,nrow=2,byrow=TRUE))

p5d
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.FacetDominance.no5-10.geneticLoad.AK.Only.NoFineLines.pdf",sep=""),p5d,height=4.8,width=3)

