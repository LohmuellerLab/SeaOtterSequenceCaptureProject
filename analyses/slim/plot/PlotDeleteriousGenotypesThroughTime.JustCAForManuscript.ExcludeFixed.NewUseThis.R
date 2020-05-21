colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
textsize=14
require(ggplot2)
require(dplyr)
todaysdate=format(Sys.Date(),format="%Y%m%d")

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/concattedSummaries/RemovingFixedVariants_useThis/"
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/loadPlots/"
dir.create(plot.dir)
#popModDates=c("AK/1D.3Epoch.LongerRecovery/20191202/","CA/1D.3Epoch.LongerRecovery/20191013/","AK/1D.5Epoch/20190802/")

popModDates=c("CA/1D.3Epoch.35Generations.LongerRecovery/2020311/")
reps=c(seq(1,25)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)
# note these simulations start at gen 0, when pop is at 4500 (taht's the only gen with that pop size)
# input: has already had avg homozygous derived calculated on hoffman, and fixed vars during burn in removed. it's ready to plot!
allAvgdInputs <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/concattedSummaries/RemovingFixedVariants_useThis/CA/1D.3Epoch.35Generations.LongerRecovery/2020311/20200311_CA_AvgHomozygousDerivedGTs.PerInd.ThroughTime.AllReps.RemovedBurninFixedVar.txt",sep="\t",header=T)
# label h:
allAvgdInputs$hLabel <- ""
allAvgdInputs[allAvgdInputs$h==0,]$hLabel <- "recessive mutations"
allAvgdInputs[allAvgdInputs$h==0.5,]$hLabel <- "additive mutations"

### need to get average across number of reps:
input_AvgdAcrossReps <- allAvgdInputs %>%
  group_by(generation,h,hLabel,population,sCat,model,popModDate,popsizeDIP) %>%
  summarise(avgHomPerInd_overAllReps=mean(avgHomPerInd))
# order factors:
input_AvgdAcrossReps$sCat <- factor(input_AvgdAcrossReps$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))
input_AvgdAcrossReps$hLabel <- ""
input_AvgdAcrossReps[input_AvgdAcrossReps$h==0,]$hLabel <- "recessive mutations"
input_AvgdAcrossReps[input_AvgdAcrossReps$h==0.5,]$hLabel <- "additive mutations"


allAvgdInputs$sCat <- factor(allAvgdInputs$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))

#modelsUsed = paste(unique(avgHomPerIndPersCat$popModDate),collapse="; ")
# aha, the "collapse" option in paste works for this to collapse list elements and sep by something:
pops=c("CA")
########## do recessive and additive separately #######
hlines=data.frame(intercept=36)
for(pop in pops){
  for(h in c(0,0.5)){
    p1 <- ggplot(allAvgdInputs[allAvgdInputs$h==h & allAvgdInputs$population==pop,],aes(y=avgHomPerInd,x=generation))+
      geom_line(size=0.3,alpha=0.5,aes(group=replicate,color=population))+
      geom_line(data=input_AvgdAcrossReps[input_AvgdAcrossReps$h==h & input_AvgdAcrossReps$population==pop,],aes(y=avgHomPerInd_overAllReps,x=generation,color=population),size=1.5)+
      theme_bw()+
      geom_vline(data=hlines,aes(xintercept=intercept),linetype="dashed")+
      facet_wrap(~sCat,nrow=1,scales="free")+
      #facet_grid(~sCat~population,scales="free_y")+
      ggtitle(unique(allAvgdInputs[allAvgdInputs$h==h,]$hLabel))+
      #ggtitle(paste("Average number of homozygous derived genotypes per individual, per category\nneutral [0]; weak (0 to -0.001]; moderate (-0.001 to -0.01]; strong (-0.01 to -1]\n",h,sep=""))+
      ylab("Avg. Hom. Derived\nGenotypes Per Individual") +
      xlab("generations since contraction")+
      scale_color_manual(values=unlist(colors))+
      theme(legend.position = 'none')
    p1
    ggsave(paste(plot.dir,"avgNumHomDerivedPerInd.ThroughTime.h.",h,".",pop,".ONLY.pdf",sep = ""),p1,width=8,height=2.5)
    
  }}

########### Plot just the contraction period ########
for(pop in pops){
  for(h in c(0,0.5)){
    p1b <- ggplot(allAvgdInputs[allAvgdInputs$h==h & allAvgdInputs$population==pop & allAvgdInputs$generation<=36,],aes(y=avgHomPerInd,x=generation))+
      geom_line(size=0.3,alpha=0.5,aes(group=replicate,color=population))+
      geom_line(data=input_AvgdAcrossReps[input_AvgdAcrossReps$h==h & input_AvgdAcrossReps$population==pop & input_AvgdAcrossReps$generation<36,],aes(y=avgHomPerInd_overAllReps,x=generation,color=population),size=1.5)+
      theme_bw()+
      facet_wrap(~sCat,nrow=2,scales="free")+
      #facet_grid(~sCat~population,scales="free_y")+
      ggtitle(unique(allAvgdInputs[allAvgdInputs$h==h,]$hLabel))+
      #ggtitle(paste("Average number of homozygous derived genotypes per individual, per category\nneutral [0]; weak (0 to -0.001]; moderate (-0.001 to -0.01]; strong (-0.01 to -1]\n",h,sep=""))+
      ylab("Avg. Hom. Derived\nGenotypes Per Individual") +
      xlab("generations since contraction")+
      scale_color_manual(values=unlist(colors))+
      theme(legend.position = 'none')
    p1b
    ggsave(paste(plot.dir,"avgNumHomDerivedPerInd.ContractionONLY.ThroughTime.h.",h,".",pop,".ONLY.pdf",sep = ""),p1b,width=4,height=2.5)
    
  }}

############## plot avg het per ind ############

avgHetPerIndPersCat_AvgAcrossReps <- allAvgdInputs %>%
  group_by(generation,h,hLabel,population,sCat,model,popModDate,popsizeDIP) %>%
  summarise(avgHetPerInd_overAllReps=mean(avgHetPerInd))
# order factors:
avgHetPerIndPersCat_AvgAcrossReps$sCat <- factor(avgHetPerIndPersCat_AvgAcrossReps$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))


pops=c("CA")
for(pop in pops){
  for(h in c(0,0.5)){
    p2 <- ggplot(allAvgdInputs[allAvgdInputs$h==h & allAvgdInputs$population==pop,],aes(y=avgHetPerInd,x=generation))+
      geom_line(size=0.3,alpha=0.5,aes(group=replicate,color=population))+
      geom_line(data=avgHetPerIndPersCat_AvgAcrossReps[avgHetPerIndPersCat_AvgAcrossReps$h==h & avgHetPerIndPersCat_AvgAcrossReps$population==pop,],aes(y=avgHetPerInd_overAllReps,x=generation,color=population),size=1.5)+
      theme_bw()+
      ggtitle(unique(allAvgdInputs[allAvgdInputs$h==h,]$hLabel))+
        #facet_grid(~sCat~population,scales="free")+
      facet_wrap(~sCat,nrow=1,scales="free")+
      #ggtitle(paste("Average number of heterozygous genotypes per individual, per category\nneutral [0]; weak (0 to -0.001]; moderate (-0.001 to -0.01]; strong (-0.01 to -1]\nh=",h,sep=""))+
      ylab("Average Heterozygous\nGenotypes Per Individual") +
      xlab("generations since contraction")+
      scale_color_manual(values=unlist(colors))+
      theme(legend.position = 'none')
    p2
    ggsave(paste(plot.dir,"avgNumHetDerivedPerInd.ThroughTime.h.",h,".",pop,".ONLY.pdf",sep = ""),p2,width=8,height=2.5)
    
  }}
##################### what else would be good? total number of alleles? #######
totalDerivedAlleles_AcrossReps <- allAvgdInputs %>%
  group_by(generation,h,hLabel,population,sCat,popModDate) %>%
  summarise(avgDerivedAcrossReps=mean(avgDerivedAllelesPerInd))

totalDerivedAlleles_AcrossReps$sCat <- factor(totalDerivedAlleles_AcrossReps$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))
########### Plot just the contraction period ########

for(pop in pops){
  for(h in c(0,0.5)){
    p3 <- ggplot(allAvgdInputs[allAvgdInputs$h==h & allAvgdInputs$population==pop  & allAvgdInputs$generation<=36,],aes(y=avgDerivedAllelesPerInd,x=generation,color=population))+
      geom_line(data=totalDerivedAlleles_AcrossReps[totalDerivedAlleles_AcrossReps$h==h & totalDerivedAlleles_AcrossReps$population==pop & totalDerivedAlleles_AcrossReps$generation<=36,],aes(y=avgDerivedAcrossReps,x=generation,color=population),size=1.5)+
      geom_line(size=0.2,alpha=0.5,aes(group=replicate,color=population))+
      theme_bw()+
      #facet_grid(~sCat~population,scales="free")+
      facet_wrap(~sCat,nrow=2,scales="free")+
      xlab("generations since contraction")+
      scale_color_manual(values=unlist(colors))+
      theme(legend.position = 'none')+
      ggtitle(unique(allAvgdInputs[allAvgdInputs$h==h,]$hLabel))+
      ylab("Average Derived\nAlleles Per Individual") 
    p3
    ggsave(paste(plot.dir,"avgDerivedAllelesPerInd.ThroughTime.h.",h,".",pop,".ONLY.ContractionOnly.pdf",sep = ""),p3,width=4,height=2.5)
    
  }}

