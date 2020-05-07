# Grid search
# also want to input fsc parameters
require(ggplot2)
populations=c("CA","AK","KUR","AL")
genotypeDate="20181119"
model="1D.2Epoch"
mu=8.64e-09
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/",genotypeDate,"/grid.search/",sep = "")
# modification from Kirk: plot relative to MLE LL rather than to data:data LL
############################### set up vectors of inferred paramaters as reference points ###########
### from dadi and fsc results, set reference points (figure out good way to make this standardized)
# draw from: /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/excelResultsComparisons/20190110/excelResultsSummary.dadi.FSCSummaryAdded.AllPops.20190110.xlsx
# set up a dictionary-like list of with the MLE results of dadi and fsc inference from above file
# contraction sizes in diploids
######## eventually write out this out as a table in a different script to automate #############
# total sites assessed: (these values are also in excel file)
Ls = vector(mode="list",length=length(populations))
Ls["CA"] <- 5989967
Ls["AK"] <- 6335196
Ls["AL"] <- 6379260
Ls["KUR"] <- 6416481

dadiNus = vector(mode="list",length=length(populations))
names(dadiNus) <- populations
dadiNus["CA"] <- 0.8
dadiNus["AK"] <- 29.5
dadiNus["AL"] <- 114.2
dadiNus["KUR"] <- 1.7

# sizes in DIPLOIDS
fscNus = vector(mode="list",length=length(populations))
names(fscNus) <- populations
fscNus["CA"] <- 309.5
fscNus["AK"] <- 389.5
fscNus["AL"] <- 1811
fscNus["KUR"] <- 1054.5

# times in generations 
dadiTs = vector(mode="list",length=length(populations))
names(dadiTs) <- populations
dadiTs["CA"] <- 0.2
dadiTs["AK"] <- 3.7
dadiTs["AL"] <- 10.2
dadiTs["KUR"] <- 0.1

fscTs = vector(mode="list",length=length(populations))
names(fscTs) <- populations
fscTs["CA"] <- 73
fscTs["AK"] <- 56
fscTs["AL"] <- 306
fscTs["KUR"] <- 80

dadiNancs = vector(mode="list",length=length(populations))
names(dadiNancs) <- populations
dadiNancs["CA"] <- 3585.3
dadiNancs["AK"] <- 4381.4
dadiNancs["AL"] <- 4834.3
dadiNancs["KUR"] <- 4348.6


fscNancs = vector(mode="list",length=length(populations))
names(fscNancs) <- populations
fscNancs["CA"] <- 3576.5
fscNancs["AK"] <- 4411.5
fscNancs["AL"] <- 4869.5
fscNancs["KUR"] <- 4289.5


############################### read in results and plot #####################

allTops <- data.frame()
for(pop in populations){
  indir=paste(data.dir,pop,"/",sep="")
  input <- read.table(paste(indir,"dadi.grid.search.",pop,".",model,".LL.output.txt.gz",sep=""),header=T,sep="\t")
  input$population <- pop
  # INSTEAD want to do it relative to MLE
  input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
  # scale results by Nanc that is inferred by dadi, rather than the given Nanc 
  # get the Nanc from theta by dividing by 4*mu*L
  # and then rescale based on that (nu * Nanc_from_theta) and (T*2*Nanc_from_theta)
  input$Nanc_from_theta <- input$theta / (4*mu*Ls[[pop]])
  input$T_scaledByNanc_from_Theta_gen <- input$T * 2* input$Nanc_from_theta
  input$nu_scaledByNanc_from_Theta_dip <- input$nu * input$Nanc_from_theta
  # get grid.search MLE coordinates: 
  MLE_nu=input[input$LL_model==max(input$LL_model),]$nu
  MLE_T=input[input$LL_model==max(input$LL_model),]$T
  MLE_Nanc=input[input$LL_model==max(input$LL_model),]$Nanc_from_theta
  # get grid.search coords that are within 1pt of MLE
  closeToMLE <- input[input$deltaLL <= 1,]
  # pull out range of parameters that are within 1 pt of MLE
  allTops = rbind(allTops,closeToMLE)
 
}


write.table(allTops,paste(data.dir,"AllGridSearchResultsWithin1PtofMle.2Epoch.AK.CA.AL.KUR.txt",sep=""),quote=F,row.names=F,sep="\t")


# get min and max ranges: 
require(dplyr)

## summarize min max based on the theta inferred by each run:
minMax <- allTops %>%
  group_by(population) %>%
  summarise(minNanc=min(Nanc_from_theta),maxNanc=max(Nanc_from_theta),minTgen=min(T_scaledByNanc_from_Theta_gen),maxTgen=max(T_scaledByNanc_from_Theta_gen),minNdip=min(nu_scaledByNanc_from_Theta_dip),maxNdip=max(nu_scaledByNanc_from_Theta_dip),minNu=min(nu),maxNu=max(nu),minT=min(T),maxT=max(T))

write.table(minMax,paste(data.dir,"AllGridSearchResultsWithin1PtofMle.2Epoch.AK.CA.AL.KUR.MinMaxRanges.txt",sep=""),quote=F,row.names=F,sep="\t")

# need the Nanc in there somehow
ggplot(allTops,aes(x=nu_scaledByNanc_from_Theta_dip,y=T_scaledByNanc_from_Theta_gen,color=population))+
  geom_point()+
  geom_point(data=allTops,aes(x=Nanc_from_theta,y=T_scaledByNanc_from_Theta_gen),color="black")+
  facet_wrap(~population)+
  scale_x_log10()+
  scale_y_log10()

allTopsInfo <-allTops %>%
  group_by(population) %>%
  summarise(Nanc=unique(Nanc_given))



mlePlot <- ggplot(allTops,aes(x=nu_scaledByNanc_from_Theta_dip,y=T_scaledByNanc_from_Theta_gen,color=population,fill=population))+
  #geom_bar(stat="identity",size=3)+
  geom_point()+
  facet_wrap(~population)+
  geom_vline(data=allTopsInfo,mapping=aes(xintercept = Nanc),color="black")+
  theme_bw()+
  ylab("Time of contraction (Generations ago)")+
  xlab("Population effective size during contraction")+
  ggtitle("Comparing the maximum likelihood estimates of a contraction model\nBlack line is the pre-contraction ancestral size\nThe contraction lasts for the number of generations on the y axis at the size on the x axis\nOnly points that are within 1 log-likelihood point of the MLE are shown")
mlePlot
ggsave(paste(data.dir,"comparingPopulationsMLEs.2Epoch.pdf",sep=""),mlePlot,height=5,width=9)



# Figure out some way to plot the more recent stuff? maybe just zoom in?

# plausible fur trade range: 1980 - 1740 = 240 years
# 240 / 6 = 40 gen; 240 / 7 = 34 gen; 240/10 = 24 gen; 240/4 = 60 gen
# so let's zoom in on 20-60 generations
mlePlotZoom <- ggplot(allTops,aes(x=nu_scaledByNanc_from_Theta_dip,y=T_scaledByNanc_from_Theta_gen,color=population,fill=population))+
  geom_point()+
  facet_wrap(~population,scales="free")+
  theme_bw()+
  ylab("Time of contraction (Generations ago)")+
  xlab("Population effective size during contraction")+
  coord_cartesian(ylim = c(20,60), xlim=c(0,1200),expand = TRUE,
                  default = FALSE, clip = "on")+ # instead of setting limits which strips data, you can just zoom in using coord_cartesian
  ggtitle("Zooming in on 20-60 Generations (reasonable time frame for fur trade)\nComparing the maximum likelihood estimates of a contraction model\nThe contraction lasts for the number of generations on the y axis at the size on the x axis\nOnly points that are within 1 log-likelihood point of the MLE are shown")
mlePlotZoom
ggsave(paste(data.dir,"comparingPopulationsMLEs.2Epoch.ZoomedIn20-60gen.pdf",sep=""),mlePlotZoom,height=5,width=9)


########### want min and max per population ###########
