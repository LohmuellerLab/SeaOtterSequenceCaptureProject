############# Likelihood ratio test for nested models from dadi (do fsc separately) ###########
pops=c("AK","AL","CA","KUR","COM")
modelDates=c("inference_20190110/1D.1Bottleneck.TB20gen","inference_20190110/1D.2Epoch","inference_20190117/1D.1Epoch")
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/"
########## get MLE for each pop-model #######
MLEs <- data.frame()
for(pop in pops){
  for(modelDate in modelDates){
    model=unlist(lapply(strsplit(modelDate,"/"),"[",2))
    input <- read.table(paste(data.dir,pop,"/",modelDate,"/",pop,".dadi.inference.",model,".all.output.concatted.txt",sep=""),header=T,sep="\t",fill = T)
    MLEdf <- data.frame(maxLL=max(input$LL)) # if there are multiple runs with same max LL, just picks top one here
    MLEdf$pop <- pop
    MLEdf$model <- model
    MLEs = rbind(MLEs,MLEdf)
}}


######## Now do Lhood ratio test for more complex to more simple models ##########
LRTResults <- data.frame()
############### 1Epoch vs 2 Epoch ############
model0="1D.1Epoch"
k0=1 # parameters for model (Nanc)
model1="1D.2Epoch"
k1=3 # parameters of model (Nanc,nu,T)
for(pop in pops){
  ln0=MLEs[MLEs$pop==pop & MLEs$model==model0,]$maxLL # ll of simpler model
  ln1=MLEs[MLEs$pop==pop & MLEs$model==model1,]$maxLL # ll of more complex model 
  df=k1-k0 # degrees of freedom is difference in number of parameters between models
  LRT=-2*(ln0-ln1)
  pval=pchisq(LRT,df,lower.tail = F) # use lower.tail = F to get the p-value otherwise it's the cumulative distribution
  LRTdf=data.frame(pop=pop,model0=model0,model1=model1,ln0=ln0,ln1=ln1,k0=k0,k1=k1,df=df,LRT=LRT,pval=pval)
  LRTResults = rbind(LRTdf,LRTResults)
}
################# 2 Epoch vs 3 Epoch ###############
model0="1D.2Epoch"
k0=3 # parameters for model (Nanc)
model1="1D.1Bottleneck.TB20gen"
k1=4 # parameters of model (Nanc,nu,T)
for(pop in pops){
  ln0=MLEs[MLEs$pop==pop & MLEs$model==model0,]$maxLL # ll of simpler model
  ln1=MLEs[MLEs$pop==pop & MLEs$model==model1,]$maxLL # ll of more complex model 
  df=k1-k0 # degrees of freedom is difference in number of parameters between models
  LRT=-2*(ln0-ln1)
  pval=pchisq(LRT,df,lower.tail = F) # use lower.tail = F to get the p-value otherwise it's the cumulative distribution
  LRTdf=data.frame(pop=pop,model0=model0,model1=model1,ln0=ln0,ln1=ln1,k0=k0,k1=k1,df=df,LRT=LRT,pval=pval)
  LRTResults = rbind(LRTdf,LRTResults)
}

# get number of tests:
numTests=dim(LRTResults)[1]
alpha_adjusted=0.05/numTests
alpha_adjusted
LRTResults$alpha_adjustedByNumTests <- alpha_adjusted
LRTResults$significant <- NA
LRTResults[LRTResults$pval < LRTResults$alpha_adjustedByNumTests,]$significant <- "yes"
LRTResults[LRTResults$pval >= LRTResults$alpha_adjustedByNumTests,]$significant <- "no"
#View(LRTResults)

write.table(LRTResults,paste(data.dir,"/LikelihoodRatioTest.LRT.ComparingAllModels.AllPops.dadi.txt",sep=""),row.names=F,col.names = T,sep="\t",quote = F)
