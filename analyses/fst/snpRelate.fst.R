#load R packages
library(gdsfmt)
library(SNPRelate)
require(ggplot2)
require(RColorBrewer)
requore(geosphere) # for distances : distm func
# tutorial: http://corearray.sourceforge.net/tutorials/SNPRelate/#ld-based-snp-pruning
# doesn't show ld pruning before hand, but I'm unsure if I should or not...
############### set up your colors -- keep this consistent across all plots ######
pops=c("CA","BAJ","AK","AL","COM","KUR")  # your populations
ALpops=c("Attu","Amchitka","Adak") # aleutian subpopulations
COMpops=c("Medny","Bering")
display.brewer.pal("Dark2",n=8)
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
# make Baja brown
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
##################### versions ###########
calldate=20181119 # date gt's were called in format YYYYMMDD
todaysdate=format(Sys.Date(),format="%Y%m%d")

########################## data.dir ##############
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/",sep="") 
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/fst/",calldate,"/",sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/fst/",calldate,"/",sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)


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

############## open gds file #############
#open the gds file
genofile <- snpgdsOpen(paste(data.dir,"/snp7/snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.gds",sep=""))

######################## population information ###########
###### *note* this popmap has bad samples in it (labeled as such)
# but don't worry, they get removed by interescting with the vcf file sample.id set!
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.PCA.txt",header=T) # this includes the RWAB samples
# exclude outliers, admixed,etc.
popmap <- popmap[popmap$notes=="good",]
sample.id=popmap$Sample
pop1_code=as.character(popmap$popCode)
pop2_code = as.character(popmap$SecondaryPop) # if you want to look at bering/medny and between aleutian islands

######### combinations of all population pairs #######
combos <- combn(pops,2)
combos

############## get sample ids from genofile ############
fstResults=data.frame(t(combos))
colnames(fstResults) <- c("pop1","pop2")
fstResults$pairwiseWeightedFst <- NA
fstResults$meanFst <- NA
for(i in seq(1,length(combos[1,]))){
  combo=combos[,i]
  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  # want to get fst for every pair
  # select samples
  selection <- popmap[(popmap$popCode %in% combo & popmap$Sample %in% sample.id),c("Sample","popCode")]
  samp.sel <- selection$Sample
  # select population ?
  pop.sel <- selection$popCode
  
  results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(as.character(pop.sel)),
                    method="W&C84",remove.monosnp = T,autosome.only = F,missing.rate = 0.2)
  fstResults[fstResults$pop1==combo[1] & fstResults$pop2==combo[2],]$pairwiseWeightedFst <- results$Fst
  fstResults[fstResults$pop1==combo[1] & fstResults$pop2==combo[2],]$meanFst <- results$MeanFst
}
############# write out results ################
write.table(fstResults,paste(fileoutdir,"pairwiseFst.MajorPops.",todaysdate,".txt",sep=""),row.names = F,quote=F,sep="\t")
# interpreting output: https://rdrr.io/bioc/SNPRelate/man/snpgdsFst.html
# output has 3 elements: test$fst -- which is weighted fst estimate (aka averaging Hs over all sites, averaging Ht over all sites, THEN dividing -- better way to estimate); meanFst which is averaged across all snps (less preferable); fstSNP a vector of fst for each snp

############# look within aleutians: #############
alcombos <- combn(ALpops,2)
alfstResults=data.frame(t(alcombos))
colnames(alfstResults) <- c("pop1","pop2")
alfstResults$pairwiseWeightedFst <- NA
alfstResults$meanFst <- NA
for(i in seq(1,length(alcombos[1,]))){
  combo=alcombos[,i]
  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  # want to get fst for every pair
  # select samples
  selection <- popmap[(popmap$SecondaryPop %in% combo & popmap$Sample %in% sample.id),c("Sample","SecondaryPop")]
  samp.sel <- selection$Sample
  # select population
  pop.sel <- selection$SecondaryPop
  results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(as.character(pop.sel)),
                       method="W&C84",remove.monosnp = T,autosome.only = F,missing.rate = 0.2)
  alfstResults[alfstResults$pop1==combo[1] & alfstResults$pop2==combo[2],]$pairwiseWeightedFst <- results$Fst
  alfstResults[alfstResults$pop1==combo[1] & alfstResults$pop2==combo[2],]$meanFst <- results$MeanFst
}
############# write out results ################
write.table(alfstResults,paste(fileoutdir,"pairwiseFst.Aleutian.",todaysdate,".txt",sep=""),row.names = F,quote=F,sep="\t")


############# look within commanders: #############
comcombos <- combn(COMpops,2)
comfstResults=data.frame(t(comcombos))
colnames(comfstResults) <- c("pop1","pop2")
comfstResults$pairwiseWeightedFst <- NA
comfstResults$meanFst <- NA
for(i in seq(1,length(comcombos[1,]))){
  combo=comcombos[,i]
  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  # want to get fst for every pair
  # select samples
  selection <- popmap[(popmap$SecondaryPop %in% combo & popmap$Sample %in% sample.id),c("Sample","SecondaryPop")]
  samp.sel <- selection$Sample
  # select population ?
  pop.sel <- selection$SecondaryPop
  results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(as.character(pop.sel)),
                       method="W&C84",remove.monosnp = T,autosome.only = F,missing.rate = 0.2)
  comfstResults[comfstResults$pop1==combo[1] & comfstResults$pop2==combo[2],]$pairwiseWeightedFst <- results$Fst
  comfstResults[comfstResults$pop1==combo[1] & comfstResults$pop2==combo[2],]$meanFst <- results$MeanFst
}
############# write out results ################
write.table(comfstResults,paste(fileoutdir,"pairwiseFst.Commander.",todaysdate,".txt",sep=""),row.names = F,quote=F,sep="\t")


##################################### try just CA and Baja with downsampled CA ###############
ca_baj = c("CA_down","BAJ_down")
downsampled_fstResults=data.frame(t(ca_baj))
colnames(downsampled_fstResults) <- c("pop1","pop2")
downsampled_fstResults$pairwiseWeightedFst <- NA
downsampled_fstResults$meanFst <- NA
# randomly picked 2 ca samples:
samp.sel <- c("115_Elut_CA_305","114_Elut_CA_214","168_Elut_BAJ_TS2","169_Elut_BAJ_R1")
pop.sel <- c("CA","CA","BAJ","BAJ")
results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(as.character(pop.sel)),
                     method="W&C84",remove.monosnp = T,autosome.only = F,missing.rate = 0.2)
downsampled_fstResults$pairwiseWeightedFst <- results$Fst
downsampled_fstResults$meanFst <- results$MeanFst
############# write out results ################
write.table(comfstResults,paste(fileoutdir,"pairwiseFst.CA_BAJ.randomlyDownsampledCa.toTwoIndividuals.115_Elut_CA_305.114_Elut_CA_214.",todaysdate,".txt",sep=""),row.names = F,quote=F,sep="\t")

############ try combos of all aleutians with all commanders -- see if something recolonized ###########
alcomcombos <- combn(c(ALpops,COMpops),2)
alcomfstResults=data.frame(t(alcomcombos))
colnames(alcomfstResults) <- c("pop1","pop2")
alcomfstResults$pairwiseWeightedFst <- NA
alcomfstResults$meanFst <- NA
for(i in seq(1,length(alcomcombos[1,]))){
  combo=alcomcombos[,i]
  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  # want to get fst for every pair
  # select samples
  selection <- popmap[(popmap$SecondaryPop %in% combo & popmap$Sample %in% sample.id),c("Sample","SecondaryPop")]
  samp.sel <- selection$Sample
  # select population ?
  pop.sel <- selection$SecondaryPop
  results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(as.character(pop.sel)),
                       method="W&C84",remove.monosnp = T,autosome.only = F,missing.rate = 0.2)
  alcomfstResults[alcomfstResults$pop1==combo[1] & alcomfstResults$pop2==combo[2],]$pairwiseWeightedFst <- results$Fst
  alcomfstResults[alcomfstResults$pop1==combo[1] & alcomfstResults$pop2==combo[2],]$meanFst <- results$MeanFst
}
write.table(alcomfstResults,paste(fileoutdir,"pairwiseFst.AleutianIslands.CommanderIslands.specific.",todaysdate,".txt",sep=""),row.names = F,quote=F,sep="\t")

####################### combos of all secondary locations ###################

secondaryPops <- unique(as.character(unlist(popmap$SecondaryPop)))
# want to treat Baja as one population: 
secondaryPops_combineBajasecondaryPops <- secondaryPops[!secondaryPops %in% c("ElRosario","TodosSantos")]
secondaryPops_combineBajasecondaryPops <- c(secondaryPops_combineBajasecondaryPops,"Baja")
secondaryCombos <- combn(secondaryPops_combineBajasecondaryPops,2)
fstResults=data.frame(t(secondaryCombos))
colnames(fstResults) <- c("pop1","pop2")
fstResults$pairwiseWeightedFst <- NA
fstResults$meanFst <- NA
fstResults$approxDistance_km <- NA

# set up popmap so that ElRosario and TodosSantos are Just Baja
# only use good samples:

popmap2 <- popmap
head(popmap2)
popmap2$locationPop <- as.character(popmap2$SecondaryPop) #population to use for distance calculation

popmap2[popmap2$SecondaryPop %in% c("TodosSantos","ElRosario"),]$locationPop <- ("Baja")

# and use SanNicolas for California locationPop


for(i in seq(1,length(secondaryCombos[1,]))){
  combo=secondaryCombos[,i]
  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  # want to get fst for every pair
  # select samples
  selection <- popmap2[(popmap2$locationPop %in% combo & popmap2$Sample %in% sample.id),c("Sample","locationPop")]
  samp.sel <- selection$Sample
  # select population 
  pop.sel <- selection$locationPop
  
  results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(as.character(pop.sel)),
                       method="W&C84",remove.monosnp = T,autosome.only = F,missing.rate = 0.2)
  fstResults[fstResults$pop1==combo[1] & fstResults$pop2==combo[2],]$pairwiseWeightedFst <- results$Fst
  fstResults[fstResults$pop1==combo[1] & fstResults$pop2==combo[2],]$meanFst <- results$MeanFst
  # add in approximate distance calculated from approx lat /long (note: using monterey as california representative)
  fstResults[fstResults$pop1==combo[1] & fstResults$pop2==combo[2],]$approxDistance_km <- as.numeric(distances_melt[distances_melt$pop1==combo[1] & distances_melt$pop2==combo[2],]$approxDistance_km)
}

head(fstResults)
write.table(fstResults,paste(fileoutdir,"pairwiseFst.allpopsWithDistances.km.",todaysdate,".txt",sep=""),row.names = F,quote=F,sep="\t")

###################### Plot FST vs Distance ##############################
# want to get FST for all the centralLocations
# make - fst 0:
fstResults$pairwiseWeightedFst_adjusted <- fstResults$pairwiseWeightedFst
fstResults[fstResults$pairwiseWeightedFst<0,]$pairwiseWeightedFst_adjusted <- 0
fstResults$bajaCA <- "no"
fstResults[fstResults$pop1 =="Baja" & fstResults$pop2 =="California" | fstResults$pop2 =="Baja" & fstResults$pop1 =="California", ]$bajaCA <- "Baja-California"
fstVsDistance <- ggplot(fstResults,aes(x=approxDistance_km,y=pairwiseWeightedFst_adjusted,color=bajaCA))+
  geom_point()+
  ggtitle("Fst vs Distance (km)\nNote: CA distances calc'd from Monterey; Baja is averaged; neg. fst set to 0")+
  xlab("Approximate Distance Between Populations (km)")+
  ylab("Pairwise Fst (weighted)")+
  geom_smooth()+
  geom_point(aes(y=0,x=distances_melt[distances_melt$pop1=="SanNicolas" & distances_melt$pop2=="Baja",]$approxDistance_km),color="red",shape=4)+
  geom_text(aes(y=-0.01,x=distances_melt[distances_melt$pop1=="SanNicolas" & distances_melt$pop2=="Baja",]$approxDistance_km),label="SanNic-Baja",color="red",size=5)+
  geom_text(aes(y=-0.01,x=distances_melt[distances_melt$pop1=="California" & distances_melt$pop2=="Baja",]$approxDistance_km+800),label="Monterey-Baja",color="black",size=5)+
  theme_bw()+
  theme(legend.position = "none")
  
fstVsDistance
ggsave(paste(plotoutdir,"/fstVsDistance.BajaAvgd.CaMonterey.andSanNic.labels.pdf",sep=""),fstVsDistance,width=10,height=6,device="pdf")
# need to figure out how to deal with the Bajas / San Nic distance. Currently it is Monterey_Baja distance. 

############## close gds file #############

snpgdsClose(genofile)
