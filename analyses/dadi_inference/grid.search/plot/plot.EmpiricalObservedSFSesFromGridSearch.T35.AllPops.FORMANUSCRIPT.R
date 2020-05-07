######################### Plot best-fit SFSes from grid-search with T=35 gen ###########
# get sfses from grid.search (made from previous script plot.EmpiricalObservedSFSesFromGridSearch.T35.AllPops.FORMANUSCRIPT.R)
require(ggplot2)
require(reshape2)
allSFSes <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/BestFit.35.ExpSFS.ObsSFS.AllPops.txt",sep="\t",header=T)
head(allSFSes)
######### need to get 1D sfses #######
pops=c("AK","AL","CA","KUR","COM")
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/"
modelDate="inference_20190117/1D.1Epoch"
# need to fold the output because outputted it unfolded (exclude monomorphic)
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
all1EpochSFSes <- data.frame()
for(pop in pops){
    model="1D.1Epoch"
    # get sfs from run 1 (the same for all runs)
    sfs <- read.table(paste(data.dir,pop,"/",modelDate,"/",pop,".dadi.inference.",model,".runNum.1.20190117.expSFS",sep=""),header=T,sep="\t",fill = T)
    dadiOutput <- read.table(paste(data.dir,pop,"/",modelDate,"/",pop,".dadi.inference.",model,".runNum.1.20190117.output",sep=""),header=T,sep="\t",fill=T)
  # need to fold it 
    theta=dadiOutput$theta
    UnfoldedsfsDF <- data.frame(countScaledByTheta1=as.numeric(unlist(strsplit(as.character(sfs[1,])," "))),mask=unlist(strsplit(as.character(sfs[2,])," ")))
    # fold it:
    UnfoldedsfsDF$frequency <- as.numeric(seq(0,length(UnfoldedsfsDF$countScaledByTheta1)-1))
    FoldedsfsDF <- foldSFS(UnfoldedsfsDF)
    # checked it -- this worked
    FoldedsfsDF_noMono <- FoldedsfsDF[FoldedsfsDF$frequency!=0,]
    FoldedsfsDF_noMono$countScaledByThetaFull <- FoldedsfsDF_noMono$count*theta
    FoldedsfsDF_noMono$pop <- pop
    all1EpochSFSes <- rbind(all1EpochSFSes,FoldedsfsDF_noMono)
  }
all1EpochSFSes_melt <- melt(all1EpochSFSes[,c("frequency","pop","countScaledByThetaFull")],id.vars =  c("frequency","pop"))
all1EpochSFSes_melt$label <- "1 Epoch Model"
########################################

allSFSes_melt <- melt(allSFSes,id.vars = c("frequency","pop"))
allSFSes_melt$label <- NA
allSFSes_melt[allSFSes_melt$variable=="countScaledByThetaFull",]$label <- "Best Fit Model"
allSFSes_melt[allSFSes_melt$variable=="obsCounts",]$label <- "Empirical"

# combine with 1Epoch:
allSFSes_melt_plus1Epoch <- rbind(all1EpochSFSes_melt,allSFSes_melt[allSFSes_melt$variable %in% c("countScaledByThetaFull","obsCounts"),])

allSFSes_melt_plus1Epoch$label <- factor(allSFSes_melt_plus1Epoch$label, levels=c("1 Epoch Model","Best Fit Model","Empirical"))
allSFSes_melt_plus1Epoch$pop <- factor(allSFSes_melt_plus1Epoch$pop,levels=c("CA","AK","AL","KUR","COM"))
sfsplot1 <- ggplot(allSFSes_melt_plus1Epoch,aes(x=as.factor(frequency),y=value,fill=label))+
  geom_bar(stat="identity",position="dodge")+
  theme_bw()+
  facet_wrap(~pop,scales="free",ncol=1)+
  ylab("count")+
  xlab("frequency")+
  scale_fill_manual(values=c("blue","darkred","darkgray"))
sfsplot1
ggsave(paste(data.dir,"PlotComparingEmpirical.1Epoch.BestFitModel.AllPops.Comis3Epoch.T35ResultsfromGridSearch.ForSI.pdf",sep=""),sfsplot1,device="pdf",height=8,width=4)
