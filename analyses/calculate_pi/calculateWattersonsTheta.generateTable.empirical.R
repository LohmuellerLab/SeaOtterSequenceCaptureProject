########## Want to calculate pi per site and theta for each population ########
# 20181026: this has been updated to NOT use Tanya's popgentools scripts which require data to be sep. by chromosome. INSTEAD, am using my "generate1DSFS.py" script for the SFS and for theta, and am using vcftools --site-pi calculations on my neutral population vcfs for pi.

require(ggplot2)
require(dplyr)
require(gridExtra)
pops=c("CA","AK","AL","COM","KUR")  # your populations
rundate=20180806 # genotype calling date
sfsdate=20181019 # date sfses were made
############### set up your colors -- keep this consistent across all plots ######
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=length(pops),name = "Dark2")
colors=list(CA=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
################### set up dir structure  ####################
todaysdate=format(Sys.Date(),format="%Y%m%d") # date you make plots

data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/",rundate,"/",sep="")
sfs.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",rundate,"/neutralSFS/",sep="")
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/PI_THETA/",rundate,"/",sep="")

dir.create(plot.dir,recursive = T,showWarnings = F)

prefix="neutral_all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled" # for pi files
sfsPrefix="all_9"
################ Get Sample Size from SFS ###############
########### THE SFS MUST BE FOLDED FOR THIS TO WORK! 
# get sample sizes from the SFS; you calculate pi and theta from the same pop-specific vcf files
allSFS_list=list()
allSFS_length_list=list()


# this is a nice for-loop -- want to EXCLUDE monomorphic bin!!!
for(i in (1:length(pops))){
  pop=pops[i]
  sfs <- read.table(paste(sfs.dir,pop,".",sfsPrefix,".unfolded.",suffix,sep=""),header=T,stringsAsFactors = F) 
  sfs$population <- pop
  popcolor=colors[pop]
  allSFS_list[[i]] <- sfs[sfs$frequency!=0,] # add the SFS to your list
}
allSFS <- bind_rows(allSFS_list)
#### make sure monomorphic bin has been excluded, otherwise SS will be +1 #### 
ss <- allSFS %>%
  group_by(population) %>%
  tally()
colnames(ss) <- c("population","numDiploidIndividuals")
#################### calculate pi per site #####################
#### based on vcftools per site output ### 
allPi <- data.frame()
for(i in (1:length(pops))){
  pop=pops[i]
  pi <- read.table(paste(data.dir,pop,"_",prefix,".sites.pi.gz",sep=""),header=T,stringsAsFactors = F) 
  piTotal=sum(pi$PI)
  callableSitesTotal <- length(pi$PI)
  #print(callableSitesTotal==length(pi$PI))
  piPerSite=(piTotal/callableSitesTotal)
  allPi <- rbind.data.frame(allPi,(cbind.data.frame(pop,piPerSite,piTotal,callableSitesTotal)))
  # add the pi per site to your dataframe
}
write.table(allPi,paste(data.dir,"allPops.neutral.piPerSitePerPopulation.",todaysdate,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
# temp plot (make a separate real plotting script)
piPlot <- ggplot(allPi,aes(x=pop,y=piPerSite,fill=pop))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=unlist(colors))+
  theme_bw()+
  theme(legend.title = element_blank())+
  xlab("Population")+
  ggtitle("Pi Per Neutral Site")
piPlot
ggsave(paste(plot.dir,"allPops.piPerSitePerPopulation.",todaysdate,".pdf",sep=""),piPlot,device = "pdf",height=5,width=7)



#################### calculate watterson's theta #####################
# note harmonic number is sum from 1 to n of 1/n. n in this case is the number of choromosomes, not individuals, so it's 2*N if N is number of diploid inds.
# define a function to calc harm number:

#### HARMONIC NUMBER

# An = SUM(i=1-->n-1)(1/i)
# Caclulate harmonic number for 2n - 1 : 
# this function will calculate the harmonic number for a 
harmonicNumber_2nMinus1 <- function(num_indv){
  n_dip=num_indv # the numver of diploid individuals
  n_hap=2*n_dip # number of haploids
  n_harm=n_hap - 1 # number of haps minus 1 for harmonic number
  harm_number=0
  for(i in 1:(n_hap-1)){
    harm_number <- harm_number + 1/i
  }
  return(harm_number)
}

# example for 26 individuals:
#harmonicNumber_2nMinus1(26) # is the harmonic # for 2*26-1 = 51 -->  4.518813 # checked it here for 51: https://www.math.utah.edu/~carlson/teaching/calculus/harmonic.html

#### WATTERSON'S THETA (uses harmonic number function)
# this function will calculate watterson's theta using the number of diploid inds, 
# number of SNPs, and the number of callable (HQ) sites used for the calculation of # of SNPs. 
wattersons_theta <- function(num_indv,numSNPs,callablesites){
  harm_num <- harmonicNumber_2nMinus1(num_indv)
  theta <- numSNPs / harm_num
  theta_divbyL <- theta/callablesites
  return(theta_divbyL)
}


allTheta <- data.frame()
for(i in (1:length(pops))){
  pop=pops[i]
  SampleSize=ss[ss$population==pop,]$numDiploidIndividuals # this is number of individuals NOT the number of chromosomes, because my harmonic number script multiplies it by 2 to get the number of chromosomes, so you don't need to do that yourself.
  totalSNPsSum=sum(allSFS[allSFS$population==pop & allSFS$frequency!=0 & allSFS$frequency!=SampleSize,]$count) # need to exclude monomorphic!
  # get from Pi thing for now, while Total SNPs is being wonky:
  #callableSitesTotal=sum(totalSNPs$siteCount) # this only works because I used called sites in neutral regions per pop as my regions
  callableSitesTotal=allPi[allPi$pop==pop,]$callableSitesTotal
  # sample size in chromosomes for this population (which is 2x number of individuals, calculated earlier in the script)
  # this uses my watterson's theta funciton; note that even though you input number of individuals, within the function it converts that to 2*numInd which is num of chromosomes.
  ThetaW <- wattersons_theta(num_indv = SampleSize,numSNPs = totalSNPsSum,callablesites = callableSitesTotal)
  allTheta <- rbind.data.frame(allTheta,(cbind.data.frame(pop,ThetaW,SampleSize,totalSNPsSum,callableSitesTotal))) # add the pi per site to your dataframe
}

thetaWPlot <- ggplot(allTheta,aes(x=pop,y=ThetaW,fill=pop))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=unlist(colors))+
  theme_bw()+
  theme(legend.title = element_blank())+
  xlab("Population")+
  ggtitle("Watterson's Theta Per Neutral Site")
thetaWPlot
write.table(allTheta,paste(data.dir,"allPops.WattersonsThetaPerSitePerPopulation.",todaysdate,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
ggsave(paste(plot.dir,"allPops.WattersonsTheta.PerSitePerPopulation.",todaysdate,".pdf",sep=""),thetaWPlot,device = "pdf",height=5,width=7)


############ Write out population information (work in progress) ########## 
for(pop in pops){
  sink(paste(data.dir,pop,".populationInformation.pi.theta.callableSites.txt",sep=""))
  cat(paste(todaysdate,pop,"\n## This is information about each population based on vcftools --site-pi; NOT from Tanya's scripts\n\n"))
  cat(paste("Population:",pop,"\n"))
  cat(paste("Sample Size (based on SFS, excluding monomorphic bin):",ss[ss$population==pop,]$numDiploidIndividuals,"\n"))
  cat(paste("Chromosome Count (2x sample size):",2*ss[ss$population==pop,]$numDiploidIndividuals,"\n"))
  cat(paste("Callable Neutral Sites:",allPi[allPi$pop==pop,]$callableSitesTotal,"\n"))
  cat(paste("pi (per site; neutral regions only):",round(allPi[allPi$pop==pop,]$piPerSite,5),"\n"))
  cat(paste("Watterson's Theta:",round(allTheta[allTheta$pop==pop,]$ThetaW,5),"\n"))
  sink()
}
