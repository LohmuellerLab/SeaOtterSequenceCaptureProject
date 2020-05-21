require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/concattedSummaries/" # pooneh's simulations
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/loadCalcs/"
dir.create(plot.dir)
#pops=c("AK","AL","genericPop.LongerContract")
#models=c("1D.2Epoch.1.5Mb.cds")
#simdates=c(20190424,20190607)
# skipping AL "AL/1D.2Epoch.1.5Mb.cds/20190424/" and CA etc -- add those in next 
popModDates=c("AK/1D.3Epoch.LongerRecovery/20191202/","CA/1D.3Epoch.LongerRecovery/20191013/") # AK and AL have dadi parameters, genericPop has parameters based on AK MLE grid that is fur-trade relevant. ### need to come up with better classification system for this. 
#reps=c(seq(1,23))
reps=c(seq(1,25)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)
#states=c("PreContraction","PostContraction")
allAvgdInputs=data.frame()
allLoads=data.frame() 
# get avg homozygous derived per individual : (get hets too?)
for(popModDate in popModDates){
#for(rep in reps){
#print(rep)
for(h in hset){
  # check if rep exists (some have random hoffman failures)
  infile=paste(data.dir,popModDate,"/h_",h,"/replicate_",rep,".slim.output.allConcatted.summary.txt.gz",sep="")
  if(file.exists(infile)){
    input = read.table(infile,sep=",",header=T)
    # want to exclude sites that are at frequency 1 *prior* to the bottleneck
    pop= unlist(lapply(strsplit(popModDate,"/"),"[",1))
    # any sites that are at frequency 1 prior to the bottleneck (at gen 0) should be excluded from load calcs. If they are at frequency 1 only after the bottleneck they can stay to be part of load. This will be if geneartion <= bneckGen (model specific) and if the numhom==popsizeDIP. can't just go by mutid because those can be duplicated across chunks. Must go by chunk AND mutID within a replicate. Then need to remove them at all other time points.... 
    # give each mutation a unique ID that is their chunk (ie chromosome #) and mutid
    # you need this because mutIDs can be duplicated between chunks, but not within a chunk
    input$chunk.mutID <- paste(input$chunk,".",input$mutid,sep="")
    fixedToRemove <- input[(input$gen==0 & input$numhom==input$popsizeDIP),]$chunk.mutID # ~4000 sites per replicate. cool
    # should each be a unique value:
    length(unique(fixedToRemove))==length(fixedToRemove)
    #exclude those sites:
    inputWithFixedRemoved <- input[!(input$chunk.mutID %in% fixedToRemove),]
    # see how this changes:
    length(unique(input$chunk.mutID)) - length(unique(inputWithFixedRemoved$chunk.mutID))==length(fixedToRemove)# should equal total to remove TRUE
    # should be none left:
    #inputWithFixedRemoved[inputWithFixedRemoved$generation==0 & inputWithFixedRemoved$numhom==inputWithFixedRemoved$popsizeDIP,] # there should be no sites left.
    # want to get load:
    inputWithFixedRemoved$qFreq <- (inputWithFixedRemoved$numhet + (2*inputWithFixedRemoved$numhom)) / (2*inputWithFixedRemoved$popsizeDIP)
    inputWithFixedRemoved$pFreq <- 1 - inputWithFixedRemoved$qFreq
    inputWithFixedRemoved$loadComponent <- (2*h*abs(inputWithFixedRemoved$s)*inputWithFixedRemoved$qFreq*inputWithFixedRemoved$pFreq) + (abs(inputWithFixedRemoved$s)*((inputWithFixedRemoved$qFreq)^2))
    # want to categorize by s
    inputWithFixedRemoved$popModDate <- popModDate
    inputWithFixedRemoved$sCat <- NA
    # JAR categories from her 2018 paper
    inputWithFixedRemoved[inputWithFixedRemoved$s >= -1 & inputWithFixedRemoved$s < -0.01,]$sCat <- "strongly deleterious"
    inputWithFixedRemoved[inputWithFixedRemoved$s >= -0.01 & inputWithFixedRemoved$s < -0.001,]$sCat <- "moderately deleterious"
    inputWithFixedRemoved[inputWithFixedRemoved$s >= -0.001 & inputWithFixedRemoved$s < 0,]$sCat <- "weakly deleterious"
    inputWithFixedRemoved[inputWithFixedRemoved$s==0,]$sCat <- "neutral"
    model= unlist(lapply(strsplit(popModDate,"/"),"[",2))
    date= unlist(lapply(strsplit(popModDate,"/"),"[",3))
    inputWithFixedRemoved$population <- pop
    #input$state <- state
    inputWithFixedRemoved$h <- h
    inputWithFixedRemoved$model <- model
    inputWithFixedRemoved$date <- date
    inputWithFixedRemoved$replicate <- rep
    # want to get total S per generation:
    # want to get totals and avgs across all chunks per generation
    avgHomPerIndPersCat <- inputWithFixedRemoved %>%
      group_by(generation,h,population,sCat,model,replicate,popModDate,popsizeDIP) %>%
      summarise(totalNumHom=sum(numhom),totalHet=sum(numhet)) %>%
      mutate(avgHomPerInd=totalNumHom/popsizeDIP) %>%
      mutate(avgHetPerInd=totalHet/popsizeDIP) %>%
      mutate(avgDerivedAllelesPerInd=((2*totalNumHom)+totalHet)/(2*popsizeDIP))
    # don't group by sCat for load calcs:
    LoadPerGeneration <- inputWithFixedRemoved %>%
      group_by(generation,h,population,model,replicate,popModDate,popsizeDIP) %>%
      summarise(totalS=sum(loadComponent))
    LoadPerGeneration$W <- exp(-LoadPerGeneration$totalS)
    LoadPerGeneration$L  = 1 - LoadPerGeneration$W # mutation load 
    
    allAvgdInputs=rbind(allAvgdInputs,data.frame(avgHomPerIndPersCat))
    allLoads = rbind(allLoads,data.frame(LoadPerGeneration))
  }}}
  
# label h:
allAvgdInputs$hLabel <- ""
allAvgdInputs[allAvgdInputs$h==0,]$hLabel <- "h = 0 (rec.)"
allAvgdInputs[allAvgdInputs$h==0.5,]$hLabel <- "h = 0.5 (add.)"
allLoads$hLabel <- ""
allLoads[allLoads$h==0,]$hLabel <- "h = 0 (rec.)"
allLoads[allLoads$h==0.5,]$hLabel <- "h = 0.5 (add.)"
## want to write this out as a table so don't have to do it multiple times ###
write.table(allAvgdInputs,paste(outdir,"AvgHomozygousDerivedGTs.PerInd.ThroughTime.AllReps.RemovedBurninFixedVar.txt",sep=""),row.names = F,col.names = T,quote=F,sep="\t")
write.table(allLoads,paste(outdir,"LoadPerGeneration.ThroughTime.AllReps.RemovedBurninFixedVar.txt",sep=""),row.names = F,col.names = T,quote=F,sep="\t")

