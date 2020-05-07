########## want to plot proportional SFSes #######
require(reshape2)

require(ggplot2)
########## empirical sfs ########
data <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/AK-15.plusMonomorphic.sfs",skip=1,stringsAsFactors = F)
numHap=dim(data)[2]-1
numDip=numHap/2
colnames(data) <- seq(0,dim(data)[2]-1)
data_melt <- melt(data[1,])
colnames(data_melt) <- c("freq","count")
data_melt$label <- "empirical"
data_melt_noMono <- data_melt[data_melt$freq!=0,]
data_melt_noMono$proportion <- data_melt_noMono$count/sum(data_melt_noMono$count)

######## models ############
models=c("nso_model_trim_20","nso_model_trim_20_plusContraction","nso_model_trim_20_simplify")
# nso_model_trim_20 = msmc model trimmed at pt 20
# nso_model_trim_20_plusContraction = msmc model trimmed at pt 20 with dadi contraction added on
# nso_model_trim_20_simplify = msmc model highly simplified with nu = 0.67 and T = 0.1 (mimicking an approximate model of 4500 long term (weighted avg of msmc curve), followed by 1000 gen at 3000 individuals (0.67 of origianl size))
########## get plots and LLs from python script output (/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/StitchPSMC_SFS_Together/southern.sea.otter.ComparePSMC.Dadi.trim.27.Nanc.5467.py)
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/StitchPSMC_SFS_Together/NSO-AK-Comparison/"
for(model in models){
  model_sfs <- read.table(paste(data.dir,model,".PROPORTIONAL.FOLDED.expSFS.txt",sep=""),skip=1,stringsAsFactors = F)
  # get prop of data and plot together and get multinom LL
  ### plot proportionals:
  colnames(model_sfs) <- seq(0,dim(model_sfs)[2]-1)
  model_sfs_melt <- melt(model_sfs[1,])
  colnames(model_sfs_melt) <- c("freq","proportion")
  model_sfs_melt$label <- "msmc_model"
  model_sfs_melt_noMono <- model_sfs_melt[model_sfs_melt$freq!=0,]
  combo <- rbind(data_melt_noMono[,c("freq","proportion","label")],model_sfs_melt_noMono)
  ###### get LL and add to plot 
  LLs <- read.table(paste(data.dir,model,".LLs.andOptimalTheta.txt",sep=""),header=T)
  p1 <- ggplot(combo[combo$proportion!=0,],aes(x=freq,y=proportion,group=label,color=label))+
    #geom_col(position="dodge")+
    geom_point(size=2)+
    geom_line()+
    theme_bw()+
    ggtitle(paste("Proportional SFS based on ",model,sep=""))+
    annotate("text",x=4,y=max(combo$proportion)-0.05,label=paste("dadiLL = ",round(LLs$dadiLL),"\nAB Multinomial LL = ",round(LLs$AnnabelLL),"\nNanc-based theta = ",round(LLs$NancTheta),"\nDadi optimal theta = ",round(LLs$dadiOptimalTheta),sep=""),size=4)+
    scale_color_manual(values=c("blue","red"))+
    scale_y_continuous(limits=c(0,0.4))
  p1 
  ggsave(paste(data.dir,model,".PROPORTIONAL.FOLDED.exp.sfs.plot.pdf",sep=""),height=5,width=7)
}

################## Get dadi models plotted too ##############
model="Dadi best fit model (T=35) from grid search"
twoEpochGrid <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/AllGridSearchResultsWithin1PtofMle.2Epoch.AK.CA.AL.KUR.txt",header=T,sep="\t")
# this has all pops
# want CA and T ~ 35
twoEpochGrid[twoEpochGrid$population=="AK" & twoEpochGrid$T_scaledByNanc_from_Theta_gen <= 36 & twoEpochGrid$T_scaledByNanc_from_Theta_gen >= 35,]
# gives you two to choose from. Want row 588 (nu = 328, T = 35.37)
AK_BestDadiModel <- twoEpochGrid[588,] # choose row 252 
AK_BestDadiModel_expSFS <- as.numeric(unlist(strsplit(as.character(AK_BestDadiModel$expectedSFS_fold_Theta1),", "))[2:(numDip+1)]) # this is folded and is relative to theta =1 , pull out the entries 1-6 (excluding first term which is monomorphic and masked)
AK_BestDadiModel_expSFS_freq <- data.frame(AK_BestDadiModel_expSFS/sum(AK_BestDadiModel_expSFS))
colnames(AK_BestDadiModel_expSFS_freq) <- "proportion"
AK_BestDadiModel_expSFS_freq$freq <- as.numeric(seq(1,numDip))
AK_BestDadiModel_expSFS_freq$label <- "best dadi model"

####### LLs for the best model (done at the end of that python script) ###
LLs <- read.table(paste(data.dir,"/bestFitDadiModel.T35.fromGridSearch.LLs.andOptimalTheta.txt",sep=""), header=T)
######### get LL:
combo2 <- rbind(AK_BestDadiModel_expSFS_freq,data_melt_noMono[,c("label","freq","proportion")])
#### plot with empirical:
p2 <- ggplot(combo2[combo2$proportion!=0,],aes(x=as.numeric(freq),y=proportion,color=label,group=label))+
  #geom_col(position="dodge")+
  geom_point(size=2)+
  geom_line()+
  theme_bw()+
  ggtitle(paste("Proportional SFS based on ",model,sep=""))+
  scale_color_manual(values=c("red","blue"))+
  annotate("text",x=4,y=max(combo$proportion)-0.05,label=paste("dadiLL = ",round(LLs$dadiLL),"\nAB Multinomial LL = ",round(LLs$AnnabelLL),"\nDadi optimal theta = ",round(LLs$dadiOptimalTheta),sep=""),size=4)+
  scale_y_continuous(limits=c(0,0.3))+
  scale_x_continuous(breaks=c(seq(1,numDip)))+
  xlab("freq")

p2
ggsave(paste(data.dir,"BestFitDadiModel.PROPORTIONAL.FOLDED.exp.sfs.plot.pdf",sep=""),p2,height=5,width=7)

