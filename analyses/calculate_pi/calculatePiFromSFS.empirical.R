############ This calculates S, pi and Theta, and Tajima's D from SFS that is in fsc2 format from easy sfs
# it is a slightly odd format because it is in fact *FOLDED* but it has zeroes in half the bins 
# Want to get rid of those zeroes and convert to R format
######## be very careful not to use this script with an SFS in any other input format except from EASY SFS 
genotypeDate="20181119"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/fastsimcoal2-plusMonomorphic/",sep="")
out.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/empirical/",genotypeDate,"/",sep="")
dir.create(out.dir,showWarnings = T,recursive = T)
#populations=c("AK","CA","AL","BER","MED","KUR")
populations=c("AK","CA","AL","COM","KUR")
################# functions #################
############# convert sfs from fsc2 to R format ##############
convertToRFormat <- function(sfs_fscFormat){
  sfs=data.frame(t(sfs_fscFormat))
  colnames(sfs) <- c("count")
  # calculate pi: from sfs in fastsimcoal format (or dadi format?) or R format
  # input: frequency 
  # then make a function
  sfs$frequency <- as.numeric(lapply(strsplit(as.character(rownames(sfs)),"_"),"[",2))
  return(sfs)
}
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
# first function is total pi (not divided per site) for use in Tajima's D
# second is to get it per site by dividing by total sites
calculatePiFromSFS_empData_total <- function(unfoldedHaploidSampleSize,totalSimulatedSites,sfsFoldedExcludingMonomorphic){
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
  }
  # then multiply sumTotal by 1/(n choose 2)
  pi_overall = (1/(choose(n,2))) * sumTotal
  return(pi_overall)
}

calculatePiFromSFS_empData_persite <- function(unfoldedHaploidSampleSize,totalSimulatedSites,sfsFoldedExcludingMonomorphic){
  n=as.numeric(unfoldedHaploidSampleSize)
  totalSites=totalSimulatedSites
  sumTotal=0
  # need to skip over monomorphic so add 1 to each thing in the sequence
  # make sure monomorphic is excluded: (should already be, but this makes sure)
  sfsExclMono = sfsFoldedExcludingMonomorphic[sfsFoldedExcludingMonomorphic$frequency!=0,]
  #print(sfsExclMono)
  for(i in seq(1,n/2)){
    # greek letter eta (pg 16 of Wakeley) is the folded count
    eta = sfsExclMono[i,2] # second column is count 
    # from the equation 
    sumTotal=sumTotal+(i*(n-i)*eta)
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
#### also need another version of 'harm number' that has i^2 as denominator
# for use in Tajima's D denominator 
a2_forTajimasD_2nMinus1 <- function(num_indv){
  n_dip=num_indv # the numver of diploid individuals
  n_hap=2*n_dip # number of haploids
  n_harm=n_hap - 1 # number of haps minus 1 for harmonic number
  harm_number=0
  for(i in 1:(n_hap-1)){
    harm_number <- harm_number + 1/(i^2) # this is the only diff from harm number -- i^2
  }
  return(harm_number)
}
# example for 26 individuals:
#harmonicNumber_2nMinus1(26) # is the harmonic # for 2*26-1 = 51 -->  4.518813 # checked it here for 51: https://www.math.utah.edu/~carlson/teaching/calculus/harmonic.html

#### WATTERSON'S THETA (uses harmonic number function)
# this function will calculate watterson's theta using the number of diploid inds, 
# number of SNPs, and the number of callable (HQ) sites used for the calculation of # of SNPs. 
wattersons_theta_total <- function(num_indv,numSNPs,callablesites){
  harm_num <- harmonicNumber_2nMinus1(num_indv)
  theta <- numSNPs / harm_num
  return(theta)
}
wattersons_theta_perSite <- function(num_indv,numSNPs,callablesites){
  harm_num <- harmonicNumber_2nMinus1(num_indv)
  theta <- numSNPs / harm_num
  theta_divbyL <- theta/callablesites
  return(theta_divbyL)
}

################# get Tajima's D denominator ##############
TajimasDenominator <- function(num_indv,numSNPs){
  # from Wakeley EQ 4.35 pg 115
  n = num_indv*2 # haploid n -->  need haploid number of individuals in SFS (so if 19 diploids --> n=38 haploids) # note it's not n-1, it's n here! 
  S=numSNPs # number of segregating sites in SFS
  a1 <- harmonicNumber_2nMinus1(num_indv) # harmonic number -- note that you want to use num_indv becaues the harmonic number function converts from diploids to haploids for you(!!)
  a2 <- a2_forTajimasD_2nMinus1(num_indv) # similar to harm. number but with 1/(i^2) -- note that you want to use num_indv becaues the function converts from diploids to haploids for you(!!)
  e1=(1/a1)*(((n+1)/(3*(n-1)))-(1/a1))
  e2=(1/(a1^2+a2))*(((2*(n^2+n+3))/(9*n*(n-1)))-((n+2)/(n*a1))+(a2/a1^2))
  Var_pi_minus_S=(e1*S + e2*S*(S-1))
  TajimasDenominator=sqrt(Var_pi_minus_S) # is sqrt of Var_pi_minus_S
  return(TajimasDenominator)
}


####### read in SFS, convert to R format, 'fold' to get rid of the zero bins (even though it's already folded -- this is just eliminating the extra zero bins), and calculate stats ############ 
df <- data.frame(population=character(),diploid_individuals=numeric(),totalNeutralSites=numeric(),S=numeric(),pi_total=numeric(),pi_perSite=numeric(),Wattersons_theta_total=numeric(),Wattersons_theta_perSite=numeric(),Tajimas_Denominator_total=numeric(),Tajimas_D_total=numeric(),stringsAsFactors = F)
for(pop in populations){
  file=paste(data.dir,pop,"_MAFpop0.obs",sep="")
  popdf <- data.frame(population=character(),diploid_individuals=numeric(),totalNeutralSites=numeric(),S=numeric(),pi_total=numeric(),pi_perSite=numeric(),Wattersons_theta_total=numeric(),Wattersons_theta_perSite=numeric(),Tajimas_Denominator_total=numeric(),Tajimas_D_total=numeric(),stringsAsFactors = F) # specific to model 
  # specific to populatio
  input <- read.table(file,header=T,stringsAsFactors = F,skip=1) # input in "R format" by my simulations
  # convert to R format:
 sfs <-  convertToRFormat(input)
  # get rid of monomorphic sites and fold SFS 
  # fold SFS -- note it's actually already folded, but there are zeroes in the unfolded bins
 # so we are 'folding' it to eliminate those zeroes, but it's not actually changing counts
 # and it's not double-folding *** be careful here -- if those zero bins aren't present you might accidentally double-fold! 
  sfsFolded <- foldSFS(sfs)
  nDip=max(sfsFolded$frequency)
  nHap=2*nDip # get unfolded sample size of total chromosomes from FOLDED SFS
  totalSites=sum(sfsFolded$count) # must include monomorphic in this total
  # need to exclude monomorphic!
  sfsFoldedExclMono <- sfsFolded[sfsFolded$frequency!=0,]
  pi_total=calculatePiFromSFS_empData_total(nHap,totalSites,sfsFoldedExclMono) # 
  pi_perSite=calculatePiFromSFS_empData_persite(nHap,totalSites,sfsFoldedExclMono)
  # divide n by two to get diploid individuals (for input into function; gets converted to haploid within function) # no
  S=sum(sfsFoldedExclMono$count) #excludes monomorphic
  thetaW_total= wattersons_theta_total(nDip,S,totalSites)
  thetaW_perSite= wattersons_theta_perSite(nDip,S,totalSites) # note this uses Ndip in the function because the function itself converts nDip to nHap (this is a bit awkward but comes of merging functions across scripts)
  TsDenominator=TajimasDenominator(nDip,S) # note that number of diploids will get converted as part of function
  TajimasD=(pi_total-thetaW_total)/TsDenominator # denominator is sqrt of variance 
  # add entry to data frame
  popdf[1,"population"] <- pop # need to include model
  popdf[1,"diploid_individuals"] <- nDip
  popdf[1,"totalNeutralSites"] <- totalSites
  popdf[1,"pi_total"] <- pi_total # and value of pi
  popdf[1,"pi_perSite"] <- pi_perSite # and value of pi
  popdf[1,"S"] <- S
  popdf[1,"Wattersons_theta_total"] <- thetaW_total
  popdf[1,"Wattersons_theta_perSite"] <- thetaW_perSite
  popdf[1,"Tajimas_Denominator_total"] <- TsDenominator
  popdf[1,"Tajimas_D_total"] <- TajimasD
  df <- rbind(df,popdf)
# combine all popdfs:
}

# write out the df in the approrpriate indir
write.table(df,paste(out.dir,"pi.S.Theta.CalculatedFromSFSes.allpops.TajimasD.txt",sep=""),quote = F,row.names=F)

