############ This calculates the SFS from fastsimcoal format
# need to figure out if input is folded or not
# look at other scripts to get file list by population 

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/neutralSimulations/"
out.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA//simulated/",sep="")
dir.create(out.dir,showWarnings = F)
models="CA.1D.2Epoch.35Gen.200Inds/20200224/"
states=c("preContraction","postContraction")
populations=""
numRep=10 # number of replicates
totalSites=6000000 # 6Mb were simulated 
############### REQUIRES SFS in SPECIFIC FORMAT ("R.format" from my scripts)
# which is UNFOLDED, and looks like:
#frequency count
#1         0     0
#2         1   817
#3         2   439
#4         3   272
#5         4   211
#6         5   166 ...
###################### fold sfs function ############
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
####################### calculate pi from folded SFS function #######################
# write a calculation function that takes in n, and a vector of frequencies (excluding monomorphic)
# then do the processing upstream
# unfoldedHaploidSampleSize is n
calculatePiFromSFS_simData <- function(unfoldedHaploidSampleSize,totalSimulatedSites,sfsFoldedExcludingMonomorphic){
  n=as.numeric(unfoldedHaploidSampleSize)
  totalSites=totalSimulatedSites
  sumTotal=0
  # need to skip over monomorphic so add 1 to each thing in the sequence
  # make sure monomorphic is excluded: (should already be, but this makes sure)
  sfsExclMono = sfsFoldedExcludingMonomorphic[sfsFoldedExcludingMonomorphic$frequency!=0,]
  for(i in seq(1,n/2)){
    # greek letter eta (pg 16 of Wakeley) is the folded count
    eta = sfsExclMono[i,2] # second column is count 
    # from the equation
    sumTotal=sumTotal+(i*(n-i)*eta)
    #print(sumTotal)
  }
  # then multiply sumTotal by 1/(n choose 2)
  pi_overall = (1/(choose(n,2))) * sumTotal
  # then divide by total sites to get pi per site (note that monomorphic sites were also projected)
  pi_perSite = pi_overall/totalSites
  return(pi_perSite)
}

####################### functions to calculate S and Wattersons Theta ############
# watterson's theta uses n-1 haploid individuals for harmonic number
# my harmonicNumber_2nMinus1 function takes in the number of *diploid* individuals
# and converts it to haploid 
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
########### Processes SFSes ###################
#df <-data.frame(model=character(),pi=numeric(),stringsAsFactors = F) # overall dataframe
for(model in models){
  indir=paste(data.dir,model,"/allSFSes/",sep="")
  out.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/simulated/",model,"/",sep="")
  dir.create(out.dir,showWarnings = F,recursive = T)
  modeldf <- data.frame(model=character(),pi=numeric(),S=numeric(),Wattersons_theta=numeric(),state=character(),stringsAsFactors = F) # specific to model
  modelname=as.character(lapply(strsplit(model,"/"),"[",1))
  for(state in states){
    files=list.files(indir,pattern=paste(state,".slim.output.unfolded.sfs.R.format",sep=""),full.names = T) # all replicates  ### make sure this is working
    statedf = data.frame(model=character(),pi=numeric(),S=numeric(),Wattersons_theta=numeric(),state=character(),stringsAsFactors = F)
    if(length(files)!=numRep){
      print("more or less than 10 replicates! stop what you're doing!")
      break
    }
    for(i in seq(1,length(files))){
      file=files[i]
      sfs <- read.table(file,header=T,stringsAsFactors = F) # input in "R format" by my simulations
      ###### SFS is to in format count and frequency
      # get rid of monomorphic sites and fold SFS 
      # fold SFS if it isn't (how to detect if it is/isn't?)
      sfsFolded <- foldSFS(sfs)
      nDip=max(sfsFolded$frequency)
      nHap=2*nDip # get unfolded sample size of total chromosomes from FOLDED SFS
      # need to exclude monomorphic!
      sfsFoldedExclMono <- sfsFolded[sfsFolded$frequency!=0,]
      pi=calculatePiFromSFS_simData(nHap,totalSites,sfsFoldedExclMono) # can include monomorphic or not as prefer
      # divide n by two to get diploid individuals (for input into function; gets converted to haploid within function) # no
      S=sum(sfsFoldedExclMono$count) #excludes monomorphic
      theta= wattersons_theta(nDip,S,totalSites) # note this uses Ndip in the function because the function itself converts nDip to nHap (this is a bit awkward but comes of merging functions across scripts)
      # add entry to data frame
      replicate = unlist(lapply(strsplit(file,"\\."),"[",6))
      statedf[i,"model"] <- as.character(lapply(strsplit(model,"/"),"[",1)) # initialize df
      statedf$pi <- pi
      statedf$S <- S
      statedf$Wattersons_theta <- theta
      statedf$state <- state
      statedf$rep <- replicate
      modeldf <- rbind(modeldf,statedf)
    }}
  # write out the df in the approrpriate indir
  write.table(modeldf,paste(out.dir,"pi.S.Theta.CalculatedFromSFSes.allreps.txt",sep=""),quote = F,row.names=F)
}
