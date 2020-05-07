# Grid search
# also want to input fsc parameters
require(ggplot2)
require(reshape2)
populations=c("CA","AK","KUR","AL")
#populations="CA"
genotypeDate="20181119"
model="1D.2Epoch"
mu=8.64e-09
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/plotsForMS/"
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
Ls["COM"] <- 6424414
######################## Choose a time #####################
# 35 generations: this is because 1986( around samplign time) - 1741 (start of trade) = 245 years / 7 yr/gen = 35 generations. Pick this for all.

timeChoice=35 # gen
############################### read in results and plot #####################
allMLEResults <- data.frame()
allSFSes <- data.frame()
for(pop in populations){
  indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/",genotypeDate,"/grid.search/",pop,"/",sep="")
  input <- read.table(paste(indir,"dadi.grid.search.",pop,".",model,".LL.output.txt.gz",sep=""),header=T,sep="\t")
  #input$deltaLL = input$LL_model - input$LL_data # this is relative to data:data
  # INSTEAD want to do it relative to MLE 
  input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
    # scale results by Nanc that is inferred by dadi, rather than the given Nanc 
  # get the Nanc from theta by dividing by 4*mu*L
  # and then rescale based on that (nu * Nanc_from_theta) and (T*2*Nanc_from_theta)
  input$Nanc_from_theta <- input$theta / (4*mu*Ls[[pop]])
  input$T_scaledByNanc_from_Theta_gen <- input$T * 2* input$Nanc_from_theta
  input$nu_scaledByNanc_from_Theta_dip <- input$nu * input$Nanc_from_theta
  
  #### restrict to within 1 pt of MLE ####
  inputMLE <- input[input$deltaLL<=1,]
  
  #### choose a point that is closest to the chosen time (and output how close it is)
  # do that by minimizing the difference between T and the Time choice:
  inputMLE_T <- inputMLE[which.min(abs(inputMLE$T_scaledByNanc_from_Theta_gen - timeChoice)),] 
  # make sure diff is less than 1
  if(inputMLE_T$T_scaledByNanc_from_Theta_gen - timeChoice >1){
    print("You are more than 1 generation away from time choice!!")
  }
  inputMLE_T$pop <- pop

  # want to get range from grid search that is within one point of MLE
  MLE_nu_max <- max(inputMLE$nu)
  MLE_nu_min <- min(inputMLE$nu)
  inputMLE_T$maxNu <- MLE_nu_max
  inputMLE_T$minNu <- MLE_nu_min
  MLE_Nanc_max <- max(inputMLE$Nanc_from_theta)
  MLE_Nanc_min <- min(inputMLE$Nanc_from_theta)
  inputMLE_T$maxNanc <- MLE_Nanc_max
  inputMLE_T$minNanc <- MLE_Nanc_min
  # pull out exp sfs
  expSFS=unlist(strsplit(as.character(inputMLE_T$expectedSFS_fold_Theta1),", "))
  obsSFS = unlist(strsplit(as.character(inputMLE_T$observedSFS_folded),", "))
  # get rid of "None" values so that only folded SFS is left:
  sfsDF=data.frame(countScaledByTheta1=as.numeric(expSFS[grep("None",expSFS,invert = T)]),obsCounts=as.numeric(obsSFS[grep("None",obsSFS,invert=T)]))
  sfsDF$theta <- inputMLE_T$theta
  # expSFS multiply by theta
  sfsDF$countScaledByThetaFull <- sfsDF$countScaledByTheta1*sfsDF$theta
  sfsDF$frequency <- row.names(sfsDF)
  sfsDF$pop <- pop
  allMLEResults <- rbind(inputMLE_T,allMLEResults)
  allSFSes <- rbind(allSFSes,sfsDF)
}

####### p1: Nanc only ######################
# order the factor levels
allMLEResults$pop <- factor(allMLEResults$pop,levels=c("CA","AK","AL","KUR"))

p1 <- ggplot(allMLEResults,aes(x=pop,y=Nanc_from_theta,color=pop))+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=minNanc,ymax=maxNanc,width=0))+
  ggtitle("Values of Nancestral within 1 LL point of MLE\nDot corresponds to T = 35 gen\nBars are min-max from grid search for range of T")+
  theme_bw()+
  ylab("Nanc")+
  xlab("")+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none")+
  geom_text(aes(label=round(Nanc_from_theta)),hjust=-0.2,vjust=-0.3)
p1
ggsave(paste(plot.dir,"Nanc.PerPop.1PtOfMle.pdf",sep=""),p1,height=4,width=4)


####### p2: Nu only ######################

p2 <- ggplot(allMLEResults,aes(x=pop,y=nu,color=pop))+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=minNu,ymax=maxNu,width=0))+
  ggtitle("Values of nu within 1 LL point of MLE\nDot corresponds to T = 35 gen\nBars are min-max from grid search for range of T")+
  theme_bw()+
  ylab("Fraction of pop. remaining\n(Ncur/Nanc)")+
  xlab("")+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none")+
  geom_text(aes(label=paste(round(nu,2)*100,"%",sep="")),hjust=-0.3,vjust=-0.3)

p2
ggsave(paste(plot.dir,"Nu.PerPop.1PtOfMle.pdf",sep=""),p2,height=4,width=4)

##### p3: frac of pop lost (1-nu) #######
p3 <- ggplot(allMLEResults,aes(x=pop,y=1-nu,color=pop))+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=1-minNu,ymax=1-maxNu,width=0))+
  ggtitle("Values of nu within 1 LL point of MLE\nDot corresponds to T = 35 gen\nBars are min-max from grid search for range of T")+
  theme_bw()+
  ylab("Fraction of population lost\n(1- Ncur/Nanc)")+
  xlab("")+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none")+
  geom_text(aes(label=paste(round(1-nu,2)*100,"%",sep="")),vjust=-0.3,hjust=-0.4)
p3
ggsave(paste(plot.dir,"FracLost.1-nu.PerPop.1PtOfMle.pdf",sep=""),p3,height=4,width=4)

################### p4: show cascade from Nanc to Ncur ##########
#### try to show the cascade
allMLEResults_melt <- melt(allMLEResults)
# rename varialbes:
toPlot <- allMLEResults_melt[allMLEResults_melt$variable %in% c("nu_scaledByNanc_from_Theta_dip","Nanc_from_theta"),]
toPlot$label <- NA
toPlot[toPlot$variable=="nu_scaledByNanc_from_Theta_dip",]$label <- "Contraction"
toPlot[toPlot$variable=="Nanc_from_theta",]$label <- "Ancestral"

p4a <- ggplot(toPlot,aes(x=label,y=value,color=pop))+
  geom_point(size=4)+
  geom_step(group=pop,direction = "vh",size=1)+
  theme_bw()+
  xlab("")+
  ylab("Effective pop. size")+
  facet_wrap(~pop,nrow=1)+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none",strip.background = element_rect(fill="white"),text=element_text(size=14))+
  geom_text(aes(label=round(value)),hjust=-0.3,vjust=0.5)+
  ggtitle(paste("Contraction time fixed at ", timeChoice, " generations",sep=""))
  
p4a
ggsave(paste(plot.dir,"Nanc.To.Ncur.PerPop.Step.pdf",sep=""),p4a,height=4,width=7)

p4b <- ggplot(toPlot,aes(x=label,y=value,color=pop))+
  geom_point(size=4)+
  geom_line(group=pop,size=1)+
  theme_bw()+
  xlab("")+
  ylab("Effective pop. size")+
  facet_wrap(~pop,nrow=1)+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none",strip.background = element_rect(fill="white"),text=element_text(size=14))+
  geom_text(aes(label=round(value)),hjust=-0.3,vjust=0.5)+
  ggtitle(paste("Contraction time fixed at ", timeChoice, " generations",sep=""))

p4b
ggsave(paste(plot.dir,"Nanc.To.Ncur.PerPop.Line.pdf",sep=""),p4b,height=4,width=7)
################################# p5: plot COM on own ###########################################
pop='COM'
input <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/COM/dadi.grid.search.COM.1D.3Epoch.LL.output.txt.gz",header=T,sep="\t")
head(input)
# INSTEAD want to do it relative to MLE 
input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
# scale results by Nanc that is inferred by dadi, rather than the given Nanc 
# get the Nanc from theta by dividing by 4*mu*L
# and then rescale based on that (nu * Nanc_from_theta) and (T*2*Nanc_from_theta)
input$Nanc_from_theta <- input$theta / (4*mu*Ls[[pop]])
input$TB_scaledByNanc_from_Theta_gen <- input$TB * 2* input$Nanc_from_theta
input$nuB_scaledByNanc_from_Theta_dip <- input$nuB * input$Nanc_from_theta
input$TF_scaledByNanc_from_Theta_gen <- input$TF * 2* input$Nanc_from_theta
input$nuF_scaledByNanc_from_Theta_dip <- input$nuF * input$Nanc_from_theta

#### restrict to within 1 pt of MLE ####
inputMLE <- input[input$deltaLL<=1,]
inputMLE_melt <- melt(inputMLE)
inputMLE_melt$pop <- "COM"

### pull out sfs ###
# pull out exp sfs
expSFS=unlist(strsplit(as.character(inputMLE$expectedSFS_fold_Theta1),", "))
obsSFS = unlist(strsplit(as.character(inputMLE$observedSFS_folded),", "))
# get rid of "None" values so that only folded SFS is left:
sfsDF=data.frame(countScaledByTheta1=as.numeric(expSFS[grep("None",expSFS,invert = T)]),obsCounts=as.numeric(obsSFS[grep("None",obsSFS,invert=T)]))
sfsDF$theta <- inputMLE$theta
# expSFS multiply by theta
sfsDF$countScaledByThetaFull <- sfsDF$countScaledByTheta1*sfsDF$theta
sfsDF$frequency <- row.names(sfsDF)
sfsDF$pop <- pop
allSFSes <- rbind(allSFSes,sfsDF)
# rename varialbes:
toPlot2 <- inputMLE_melt[inputMLE_melt$variable %in% c("nuB_scaledByNanc_from_Theta_dip","nuF_scaledByNanc_from_Theta_dip","Nanc_from_theta","TF_scaledByNanc_from_Theta_gen","TB_scaledByNanc_from_Theta_gen"),]
toPlot2$label <- NA
toPlot2[toPlot2$variable=="nuB_scaledByNanc_from_Theta_dip",]$label <- "Contracted"
toPlot2[toPlot2$variable=="Nanc_from_theta",]$label <- "Ancestral"
toPlot2[toPlot2$variable=="nuF_scaledByNanc_from_Theta_dip",]$label <- "Recovery"
toPlot2[toPlot2$variable=="TF_scaledByNanc_from_Theta_gen",]$label <- "RecTime"
toPlot2[toPlot2$variable=="TB_scaledByNanc_from_Theta_gen",]$label <- "BotTime"

toPlot2$pop <- "COM"

# but has long TF....
p5 <- ggplot(toPlot2[!toPlot2$variable %in% c("TF_scaledByNanc_from_Theta_gen","TB_scaledByNanc_from_Theta_gen"),],aes(x=label,y=value,color=pop))+
  geom_point(size=4)+
  geom_line(group=pop,size=1)+
  #geom_step(group=pop,direction = "vh",size=1)+
  theme_bw()+
  xlab("")+
  ylab("Effective pop. size")+
  facet_wrap(~pop,nrow=1)+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none",strip.background = element_rect(fill="white"),text=element_text(size=14))+
  geom_text(aes(label=round(value)),hjust=-0.3,vjust=0.5)+
  ggtitle(paste("Contraction time fixed: ", round(inputMLE$TB_scaledByNanc_from_Theta_gen), " gen\nRec. duration inferred: ",round(inputMLE$TF_scaledByNanc_from_Theta_gen)," gen",sep=""))

p5
ggsave(paste(plot.dir,"Nanc.To.Ncur.ToNrec.Line.COMOnly.pdf",sep=""),p5,height=4.5,width=3)

################### p6: try to combine COM with others ###############
allMLEResults_plusCOM <- rbind(allMLEResults_melt,inputMLE_melt)
write.table(allMLEResults_plusCOM,paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/",genotypeDate,"/grid.search/allMLEREsults_plusCOM.TFixed.",timeChoice,".txt",sep=""),quote=F,sep="\t",row.names=F)

allMLEResults_plusCOM$label <- "DontUse"
allMLEResults_plusCOM[allMLEResults_plusCOM$variable=="nu_scaledByNanc_from_Theta_dip",]$label <- "Ncur"
allMLEResults_plusCOM[allMLEResults_plusCOM$variable=="nuB_scaledByNanc_from_Theta_dip",]$label <- "Nbot"
allMLEResults_plusCOM[allMLEResults_plusCOM$variable=="Nanc_from_theta",]$label <- "Nanc"
allMLEResults_plusCOM[allMLEResults_plusCOM$variable=="nuF_scaledByNanc_from_Theta_dip",]$label <- "Nrec"
#allMLEResults_plusCOM[allMLEResults_plusCOM$variable=="TF_scaledByNanc_from_Theta_gen",]$label <- "RecTime"
#allMLEResults_plusCOM[allMLEResults_plusCOM$variable=="TB_scaledByNanc_from_Theta_gen",]$label <- "BotTime"

ToPlot3 <- allMLEResults_plusCOM[allMLEResults_plusCOM$label!="DontUse",]
ToPlot3$bestFitModel <- "Model"
ToPlot3[ToPlot3$pop !="COM",]$bestFitModel <- "2-Epoch" 
ToPlot3[ToPlot3$pop =="COM",]$bestFitModel <- "3-Epoch" 

p6 <- ggplot(ToPlot3,aes(x=label,y=value,color=pop))+
  geom_point(aes(size=value),alpha=0.9)+
  geom_line(group=pop,size=1,alpha=0.4)+
  theme_bw()+
  xlab("")+
  ylab("Effective pop. size")+
  facet_wrap(~pop~bestFitModel,scales="free_x",nrow=1)+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none",strip.background = element_rect(fill="white"),text=element_text(size=14))+
  geom_text(aes(label=round(value)),hjust=0.5,vjust=0.5,color="black",size=rel(5))+
  #ggtitle(paste("Contraction time fixed at ~", timeChoice, " generations",sep=""))+
  scale_size_continuous(range=c(3,18))

p6
ggsave(paste(plot.dir,"Nanc.To.Ncur.IncludesCOM.TFixedAt35.Line.pdf",sep=""),p6,height=4,width=7)

###################### p6b: plot bubbles, without text (Kirk request) ###########

p6b <- ggplot(ToPlot3,aes(x=label,y=value,color=pop))+
  geom_point(aes(size=value),alpha=1)+
  geom_line(group=pop,size=1,alpha=0.4)+
  theme_bw()+
  xlab("")+
  ylab("Effective pop. size")+
  facet_wrap(~pop~bestFitModel,scales="free_x",nrow=1)+
  scale_color_manual(values=unlist(colors))+
  theme(legend.position="none",strip.background = element_rect(fill="white"),text=element_text(size=14))+
  #geom_text(aes(label=round(value)),hjust=0.5,vjust=0.5,color="black",size=rel(5))+
  #ggtitle(paste("Contraction time fixed at ~", timeChoice, " generations",sep=""))+
  scale_size_continuous(range=c(3,18))+
  scale_y_continuous(breaks=c(seq(0,6000,by=750)))

p6b
ggsave(paste(plot.dir,"Nanc.To.Ncur.IncludesCOM.TFixedAt35.Line.NoText.pdf",sep=""),p6b,height=4,width=7) # 20200401 redid  Plot with Nbot relabeled as Ncur except for COM 

################### Pull out expected and obs SFSes from each (add 1D as well) ########
write.table(allSFSes,paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/BestFit.",timeChoice,".ExpSFS.ObsSFS.AllPops.txt",sep=""),row.names = F,col.names = T,quote=F,sep="\t")
### plot in another script: plot.EmpiricalObservedSFSesFromGridSearch.T35.AllPops.FORMANUSCRIPT.R