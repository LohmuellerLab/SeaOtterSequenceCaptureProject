# the 2D SFS formatting for fsc may not be right (bug in EasySFS -- going to continue to troubleshoot) ; my logic internally here is okay, but coming from easy SFS and into fsc2 may be weird #
###################### Adding monomorphic sites to fastsimcoal SFSes ##########
# adding to fastsimcoal and dadi, and writing out total sites (L) for dadi
# this works on 1D and 2D sfses that are output from easySFS
# the wrapper for this script is in 
library("optparse")
option_list = list(
  make_option(c("-dataDir", "--dataDir"), type="character", default=NULL, 
              help="path to your data directory (projection output of easy SFS)", metavar="character"),
  make_option(c("-popFile", "--popFile"), type="character", default=NULL, 
              help="path to your data directory (projection output of easy SFS)", metavar="character"),
  make_option(c("-class", "--class"), type="character", default=NULL, 
              help="type of site (neutral, synonymous, missense)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

data.dir=opt$dataDir
pop=opt$popFile
popFile=read.table(pop,header=F)
class=opt$class
# for testing purposes only:

fsc.format.dir=paste(data.dir,"/fastsimcoal2/",sep="") # where easysfs output is
new.fsc.format.dir=paste(data.dir,"/fastsimcoal2-plusMonomorphic/",sep="") # where you'll put the new sfses that have monomorphic sites added in
dir.create(new.fsc.format.dir,showWarnings = F) 
dadi.format.dir = paste(data.dir,"/dadi/",sep="")
new.dadi.format.dir=paste(data.dir,"/dadi-plusMonomorphic/",sep="") # where you'll put the new sfses that have monomorphic sites added in
dir.create(new.dadi.format.dir,showWarnings = F) 
colnames(popFile) <- c("sample","population")
popOrder <- as.character(unique(popFile$population)) # this should be CA,AK,AL,COM,KUR  for my project

####################### read in file of monomorphic site counts that are monomorphic in all individuals and are called in at least [projection value] for each population ######
# note that sites that are variable in the whole sample, but monomorphic in one population, are already included in the zero bin of the sfs. # 

monomorph <- read.table(paste(data.dir,"/countsOfMonomorphicPassingProjectionThresholds.perPop.txt",sep=""),header = T)
monomorph2D <- read.table(paste(data.dir,"/countsOfMonomorphicPassingProjectionThresholds.perPair.txt",sep=""),header = T)


################################## 1D fastsimcoal format ################################
print("begining 1d fsc formatting")

#### go through SFSes and write them out with the monomorphic added in
for(pop in popOrder) {
  input <- list.files(fsc.format.dir,pattern=paste("^",pop,"_MAF",sep=""),full.names = T)
  sfs <- read.table(input,skip = 1,header = T) # skip first line: "1 observation"
  monoCount <- monomorph[monomorph$population==pop,]$HomREFcount
  sfs$d0_0 <- sfs$d0_0+monoCount
  sink(paste(new.fsc.format.dir,pop,"_MAFpop0.obs",sep=""))
  cat("1 observation\n")
  write.table(sfs,quote=F,row.names = F) # writes new sfs to the table
  sink()
}

################################## 1D dadi format ########################################
print("begining 1d dadi formatting")
for(pop in popOrder) {
  # pattern is pop-##.sfs
  input <- list.files(dadi.format.dir,pattern=paste("^",pop,"-[0-9]+.sfs",sep=""),full.names = T)
  #print(input)
  header <- readLines(input,n=1) # get first line of file
  sfs <- read.table(input,skip = 1,header = F) # skip first line: "13 folded "KUR""
  monoCount <- monomorph[monomorph$population==pop,]$HomREFcount
  # two rows; first is counts, second is mask; don't change the mask
  sfs[1,"V1"] <- sfs[1,"V1"]+monoCount
  totalSites <- sum(sfs[1,])
  sink(paste(new.dadi.format.dir,pop,"-",as.character(length(sfs)),".plusMonomorphic.sfs",sep=""))
  cat(header,"\n") # put the header back in
  write.table(sfs,quote=F,row.names = F,col.names = F) # writes new sfs to the table
  sink()
  sink(paste(new.dadi.format.dir,pop,"-",as.character(length(sfs)),".totalSiteCount.L.withMonomorphic.txt",sep=""))
  cat("pop\ttotalSites\n")
  cat(pop, totalSites,"\n")
  sink()
}

####################################### 2D SFS -- fastsimcoal ###################################
print("begining 2d fsc formatting")

# fsc files are numbered by pop0_1 where pop numbers are from my pop order (check this carefully)
# order should be CA,AK,AL,COM,KUR  for my project (0,1,2,3,4)
# 20190415: I think the labels may be backward for pop 0 and pop 1 , check it out! Yes confirmed: easySFS switched d0 and d1 labels (going to fix this in my script for the future) -- doesn't matter here because I am going to change the column and row names to d0 / d1 anyway. But is annoying.
for(i in seq(0,length(popOrder)-1)) { # pop i is the first pop listed will end up as pop d1 down the rows  (bc fsc wants pop1_0.obs with 0 along cols and 1 down rows)
  # the second listed population starts 
  # the numbers of j are >i until length of sequence (because once 0_1 is done, you don't do 1_0)
  for(j in seq(i+1,length(popOrder)-1)){ # pop j is second pop listed in filename -- will end up as pop d0 along the columns (bc fsc wants pop1_0.obs with 0 along cols and 1 down rows)
  # note counts are zero based in filename (eg CA is 0), but 1 based in R
  pop1=popOrder[i+1] # i is listed *first* in the filename which means it is down the rows and so is pop **1** not 0 (weird finicky fsc thing) 
  #i and j are zero based so add 1 to get right index
  pop0=popOrder[j+1] # j is listed *second* in filename which means it is across columns and is pop 0 (weird fsc formatting)
  input <- list.files(fsc.format.dir,pattern=paste(class,"_jointMAFpop",i,"_",j,sep=""),full.names = T)
  # if input is one file:
  if(length(input)==1){
    sfs <- read.table(input,skip = 1,header = T) # skip first line: "1 observation"
    # update dx and dy to be consistent, just 0 and 1 (0 along tops, 1 down sides)
    colnames(sfs) <- paste("d0_",seq(0,ncol(sfs)-1),sep="") # subtract 1 bc is zero based
    rownames(sfs) <- paste("d1_",seq(0,nrow(sfs)-1),sep="") # subtract 1 bc is zero based
    monoCount <- monomorph2D[(monomorph2D$population1==pop0 & monomorph2D$population2==pop1)|(monomorph2D$population1==pop1 & monomorph2D$population2==pop0),]$HomREFcountPassingBothProjThresholds
  # row and column names reflect pop numbers; first number in file name is colnames, second is rownames. so row is first and should be j, and col is second and should be i
    # 20190716: for fsc to run, you need the filename to be pop1_0.obs which pop d0  along the top and d1 down the side no matter what the pop IDs were from. Easy SFS messes up the d0 labels so I am fixing it so that it doesn't matter
    # in the easy SFS file name, it gives you the population numbers that are in the order you specified. This part is okay. The one that is listed first is down the side (row names) and the one that is listed second is along the rows. This is okay, but just need to be clear about which population is which. I want to keep the order of populations consistent with that so that the first listed population goes down the side (1) and second is along the top (0).
    #sfs[paste("d",j,"_0",sep=""),paste("d",i,"_0",sep = "")] <- sfs[paste("d",j,"_0",sep=""),paste("d",i,"_0",sep = "")]+monoCount # update monomorphic count --d oesn't matter
    # update monomorphic bin by adding in mono count: sfs["d1_0","d0_0"] is 0,0 bin 
    # new col/row names that are generic:
    sfs["d1_0","d0_0"] <- sfs["d1_0","d0_0"]+monoCount
    sink(paste(new.fsc.format.dir,"/neutral.",pop1,".",pop0,".NewNames_jointMAFpop1_0.obs",sep=""))
    cat("1 observation\n")
    write.table(sfs,quote=F,row.names = T) # writes new sfs to the table
    sink()
  }
  }
}


####################################### 2D SFS -- dadi ###################################
print("begining 2d dadi formatting")

# dadi files are labeled with populations, so don't need to indexing that i did for fsc above
# first entry in sfs is 0,0 bin
for(pop1 in popOrder){
  for(pop2 in popOrder){
    input <- list.files(dadi.format.dir,pattern=paste("^",pop1,"-",pop2,".sfs",sep=""),full.names = T)
    if(length(input)==1){
      header <- readLines(input,n=1) # get first line of file
      sfs <- read.table(input,skip = 1,header = F) # skip first line: "13 folded "KUR""
      monoCount <- monomorph2D[monomorph2D$population1==pop1 & monomorph2D$population2==pop2,]$HomREFcountPassingBothProjThresholds
      # two rows; first is counts, second is mask; don't change the mask
      # 0,0 bin is first entry
      sfs[1,"V1"] <- sfs[1,"V1"]+monoCount
      totalSites <- sum(sfs[1,])
      sink(paste(new.dadi.format.dir,pop1,"-",pop2,".plusMonomorphic.sfs",sep=""))
      cat(header,"\n") # put the header back in
      write.table(sfs,quote=F,row.names = F,col.names = F) # writes new sfs to the table
      sink()
      sink(paste(new.dadi.format.dir,pop1,"-",pop2,".totalSiteCount.L.withMonomorphic.txt",sep=""))
      cat("pop1\tpop2\ttotalSites\n")
      cat(pop1, pop2, totalSites,"\n")
      sink()
      
    }
  }
}

