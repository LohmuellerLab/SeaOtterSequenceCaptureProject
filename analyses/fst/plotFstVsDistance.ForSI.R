require(ggplot2)
plotoutdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/fst/20181119/"
############# get approximate distances###########
#################### get approx distances #############
###### get combinations of all locations 
approxLatLongOfSamples <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/SamplingLatLong/approxLatLongOfSamples.txt",header=T,sep="\t",stringsAsFactors = F)
# add in Baja average -- averaging lat/long actually doesn't make sense
# need to average distances properly using geomean -- needs to take earth into account can't just average numbers. used geomean instead!
BajAvg <- c(Population="BAJ",centralLocation="Baja",latConverted=geomean(approxLatLongOfSamples[approxLatLongOfSamples$Population=="BAJ",][4:3])[2],longConverted=geomean(approxLatLongOfSamples[approxLatLongOfSamples$Population=="BAJ",][4:3])[1])
approxLatLongOfSamples2 <- rbind.data.frame(approxLatLongOfSamples,BajAvg,stringsAsFactors = F) # adding in baja average
# make numeric:
approxLatLongOfSamples2$latConverted <- as.numeric(approxLatLongOfSamples2$latConverted)
approxLatLongOfSamples2$longConverted <- as.numeric(approxLatLongOfSamples2$longConverted)

# get a distance matrix from geospheres: 
distances <- distm(approxLatLongOfSamples2[4:3],approxLatLongOfSamples2[4:3])/1000 # output is in meters, so am dividing by 1000 to get km.
colnames(distances) <- approxLatLongOfSamples2$centralLocation
rownames(distances) <- approxLatLongOfSamples2$centralLocation
distances_melt <- melt(distances)
colnames(distances_melt) <- c("pop1","pop2","approxDistance_km")
# okay so now I have a distance matrix! Need an fst matrix too.

 ##################### read in fsts ################
calldate=20181119 # date gt's were called in format YYYYMMDD
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/fst/",calldate,"/",sep="")

fstResults <- read.table(paste(fileoutdir,"pairwiseFst.allpopsWithDistances.km.20181219.txt",sep=""),sep="\t",header=T)

# make - fst 0:
fstResults$pairwiseWeightedFst_adjusted <- fstResults$pairwiseWeightedFst
fstResults[fstResults$pairwiseWeightedFst<0,]$pairwiseWeightedFst_adjusted <- 0
fstResults$bajaCA <- "no"
fstResults[fstResults$pop1 =="Baja" & fstResults$pop2 =="California" | fstResults$pop2 =="Baja" & fstResults$pop1 =="California", ]$bajaCA <- "Baja-California"

fstVsDistance1 <- ggplot(fstResults,aes(x=approxDistance_km,y=pairwiseWeightedFst_adjusted,color=bajaCA))+
  geom_point()+
  #ggtitle("Fst vs Distance (km)\nNote: CA distances calc'd from Monterey; Baja is averaged; neg. fst set to 0")+
  xlab("Approximate Distance Between Populations (km)")+
  ylab("Pairwise Fst (weighted)")+
  geom_smooth()+
  theme_bw()+
  theme(legend.position = "none")

fstVsDistance1

fstVsDistance2 <- ggplot(fstResults,aes(x=approxDistance_km,y=pairwiseWeightedFst_adjusted))+
  geom_point()+
  #ggtitle("Fst vs Distance (km)\nNote: CA distances calc'd from Monterey; Baja is averaged; neg. fst set to 0")+
  xlab("Approximate Distance Between Populations (km)")+
  ylab("Pairwise Fst (weighted)")+
  geom_smooth()+
  theme_bw()+
  theme(legend.position = "none")

fstVsDistance2

ggsave(paste(plotoutdir,"/fstVsDistance.NoBajalabels.ForSI.pdf",sep=""),fstVsDistance2,width=4,height=3,device="pdf")
