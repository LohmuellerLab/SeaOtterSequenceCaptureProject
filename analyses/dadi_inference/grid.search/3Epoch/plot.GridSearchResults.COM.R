### Plot commander 3D grid search
require(ggplot2)
require(scales)
populations=c("COM")
genotypeDate="20181119"
model="1D.3Epoch"
mu=8.64e-09
# later try rayshader
#################### 1. TB fixed at ~30 gen, TF nuF nuB inferred #######
input <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/COM/dadi.grid.search.COM.1D.3Epoch.LL.output.txt.gz",header=T,sep="\t")
head(input)
names(input)
### need to get L and mu:
L=6424414  # from file COM-35.totalSiteCount.L.withMonomorphic.txt 
input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
input$Nanc_from_theta <- input$theta / (4*mu*L)
input$TF_scaledByNanc_from_Theta_gen <- input$TF * 2* input$Nanc_from_theta
input$nuB_scaledByNanc_from_Theta_dip <- input$nuB * input$Nanc_from_theta
input$nuF_scaledByNanc_from_Theta_dip <- input$nuF * input$Nanc_from_theta

# get grid.search MLE coordinates: 
MLE_nuB=input[input$LL_model==max(input$LL_model),]$nuB
MLE_nuF=input[input$LL_model==max(input$LL_model),]$nuF
MLE_TF=input[input$LL_model==max(input$LL_model),]$TF
MLE_Nanc=input[input$LL_model==max(input$LL_model),]$Nanc_from_theta
MLE_nuRatio=input[input$LL_model==max(input$LL_model),]$nuB/input[input$LL_model==max(input$LL_model),]$nuF
## rescale them to real units:
MLE_nuB_rescale = MLE_nuB*MLE_Nanc
MLE_nuF_rescale = MLE_nuF*MLE_Nanc
MLE_TF_rescale = MLE_TF*2*MLE_Nanc

######## get ratio of nuB over nuF :
input$nuB_over_nuF_ratio = input$nuB / input$nuF

p1a <- ggplot(input,aes(nuB_over_nuF_ratio,TF_scaledByNanc_from_Theta_gen,color=deltaLL))+
  geom_point(size=2,alpha=0.5)+
  scale_x_log10()+
  scale_y_log10(label=comma)+
  scale_color_gradientn(breaks=c(30,10,1),limits=c(0,30),colors =c("green","yellow","darkorange","gray","darkgray"))+
  geom_point(aes(x=MLE_nuRatio,y=MLE_TF_rescale),shape=0,size=6,color="red")+
  labs(fill = element_text("deltaLL\nrelative to MLE")) #+
p1a

# get range of valid parameters:
allTops <- input[input$deltaLL <= 5,]
minMax <- allTops %>%
  summarise(minNanc=min(Nanc_from_theta),maxNanc=max(Nanc_from_theta),minTFgen=min(TF_scaledByNanc_from_Theta_gen),maxTFgen=max(TF_scaledByNanc_from_Theta_gen),minNBdip=min(nuB_scaledByNanc_from_Theta_dip),maxNBdip=max(nuB_scaledByNanc_from_Theta_dip),minNFdip=min(nuF_scaledByNanc_from_Theta_dip),maxNFdip=max(nuF_scaledByNanc_from_Theta_dip),minNuB=min(nuB),maxNu=max(nuB),minNuF=min(nuF),maxNuF=max(nuF),minTF=min(TF),maxT=max(TF))

write.table(minMax,paste(data.dir,"AllGridSearchResultsWithin1PtofMle.3Epoch.MinMaxRanges.COM.txt",sep=""),quote=F,row.names=F,sep="\t")
######### forcing TB and TF to be recent blows up the ancestral size like crazy!!! not good. #### could this be due to kamchatka admixture?
################# 2. try with fixing TB (~20 gen) and TF (~10 gen) #######
input <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/COM/dadi.grid.search.COM.1D.3Epoch.FixedTimes.LL.output.txt.gz",header=T,sep="\t")
input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
input$Nanc_from_theta <- input$theta / (4*mu*L)
input$TF_scaledByNanc_from_Theta_gen <- input$TF * 2* input$Nanc_from_theta
input$TB_scaledByNanc_from_Theta_gen <- input$TB * 2* input$Nanc_from_theta
# these are the estimated parameters:
input$nuB_scaledByNanc_from_Theta_dip <- input$nuB * input$Nanc_from_theta
input$nuF_scaledByNanc_from_Theta_dip <- input$nuF * input$Nanc_from_theta

# get grid.search MLE coordinates: 
MLE_nuB=input[input$LL_model==max(input$LL_model),]$nuB
MLE_nuF=input[input$LL_model==max(input$LL_model),]$nuF
MLE_TF=input[input$LL_model==max(input$LL_model),]$TF
MLE_Nanc=input[input$LL_model==max(input$LL_model),]$Nanc_from_theta
MLE_nuRatio=input[input$LL_model==max(input$LL_model),]$nuB/input[input$LL_model==max(input$LL_model),]$nuF
## rescale them to real units:
MLE_nuB_rescale = MLE_nuB*MLE_Nanc
MLE_nuF_rescale = MLE_nuF*MLE_Nanc
MLE_TF_rescale = MLE_TF*2*MLE_Nanc

################# 3. try with letting TB vary and fixing TF (~10 gen) #######
input <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/COM/dadi.grid.search.COM.1D.3Epoch.FixTFOnly.LL.output.txt.gz",header=T,sep="\t")

### need to get L and mu:
L=6424414  # from file COM-35.totalSiteCount.L.withMonomorphic.txt 
input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
input$Nanc_from_theta <- input$theta / (4*mu*L)
input$TB_scaledByNanc_from_Theta_gen <- input$TB * 2* input$Nanc_from_theta
input$nuB_scaledByNanc_from_Theta_dip <- input$nuB * input$Nanc_from_theta
input$nuF_scaledByNanc_from_Theta_dip <- input$nuF * input$Nanc_from_theta

# get grid.search MLE coordinates: 
MLE_nuB=input[input$LL_model==max(input$LL_model),]$nuB
MLE_nuF=input[input$LL_model==max(input$LL_model),]$nuF
MLE_TB=input[input$LL_model==max(input$LL_model),]$TB
MLE_Nanc=input[input$LL_model==max(input$LL_model),]$Nanc_from_theta
MLE_nuRatio=input[input$LL_model==max(input$LL_model),]$nuB/input[input$LL_model==max(input$LL_model),]$nuF
## rescale them to real units:
MLE_nuB_rescale = MLE_nuB*MLE_Nanc
MLE_nuF_rescale = MLE_nuF*MLE_Nanc
MLE_TB_rescale = MLE_TB*2*MLE_Nanc

######## get ratio of nuB over nuF :
input$nuB_over_nuF_ratio = input$nuB / input$nuF

p3a <- ggplot(input,aes(nuB_over_nuF_ratio,TB_scaledByNanc_from_Theta_gen,color=deltaLL))+
  geom_point(size=2,alpha=0.5)+
  scale_x_log10()+
  scale_y_log10(label=comma)+
  scale_color_gradientn(breaks=c(30,10,1),limits=c(0,30),colors =c("green","yellow","darkorange","gray","darkgray"))+
  geom_point(aes(x=MLE_nuRatio,y=MLE_TB_rescale),shape=0,size=6,color="red")+
  labs(fill = element_text("deltaLL\nrelative to MLE")) #+
p3a

###################### 4. 2 epoch with singletons ###############
input <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/COM/dadi.grid.search.COM.1D.2Epoch.LL.output.txt.gz",header=T,sep="\t")

### need to get L and mu:
L=6424414  # from file COM-35.totalSiteCount.L.withMonomorphic.txt (make sure you didn't use COM-31)
input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
input$Nanc_from_theta <- input$theta / (4*mu*L)
input$T_scaledByNanc_from_Theta_gen <- input$T * 2* input$Nanc_from_theta
input$nu_scaledByNanc_from_Theta_dip <- input$nu * input$Nanc_from_theta

# get grid.search MLE coordinates: 
MLE_nu=input[input$LL_model==max(input$LL_model),]$nu
MLE_T=input[input$LL_model==max(input$LL_model),]$T
MLE_Nanc=input[input$LL_model==max(input$LL_model),]$Nanc_from_theta

## rescale them to real units:
MLE_nu_rescale = MLE_nu*MLE_Nanc
MLE_T_rescale = MLE_T*2*MLE_Nanc

######## get ratio of nuB over nuF :

p4a <- ggplot(input,aes(x=nu_scaledByNanc_from_Theta_dip,y=T_scaledByNanc_from_Theta_gen,color=deltaLL))+
  geom_point(size=2,alpha=0.5)+
  scale_x_log10()+
  scale_y_log10(label=comma)+
  scale_color_gradientn(breaks=c(30,10,1),limits=c(0,30),colors =c("green","yellow","darkorange","gray","darkgray"))+
  labs(fill = element_text("deltaLL\nrelative to MLE")) #+
p4a
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/COM/2Epoch.RealUnits.pdf",p4a,height=5,width=7)


###################### 4. 2 epoch with singletons ###############
input <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/COM-noSingletons/dadi.grid.search.COM.1D.2Epoch.LL.output.txt.gz",header=T,sep="\t")

### need to get L and mu:
L=6424414  # from file COM-35.totalSiteCount.L.withMonomorphic.txt (make sure you didn't use COM-31)
input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
input$Nanc_from_theta <- input$theta / (4*mu*L)
input$T_scaledByNanc_from_Theta_gen <- input$T * 2* input$Nanc_from_theta
input$nu_scaledByNanc_from_Theta_dip <- input$nu * input$Nanc_from_theta

# get grid.search MLE coordinates: 
MLE_nu=input[input$LL_model==max(input$LL_model),]$nu
MLE_T=input[input$LL_model==max(input$LL_model),]$T
MLE_Nanc=input[input$LL_model==max(input$LL_model),]$Nanc_from_theta

## rescale them to real units:
MLE_nu_rescale = MLE_nu*MLE_Nanc
MLE_T_rescale = MLE_T*2*MLE_Nanc

######## get ratio of nuB over nuF :

p5a <- ggplot(input,aes(x=nu_scaledByNanc_from_Theta_dip,y=T_scaledByNanc_from_Theta_gen,color=deltaLL))+
  geom_point(size=2,alpha=0.5)+
  scale_x_log10()+
  scale_y_log10(label=comma)+
  scale_color_gradientn(breaks=c(30,10,1),limits=c(0,30),colors =c("green","yellow","darkorange","gray","darkgray"))+
  labs(fill = element_text("deltaLL\nrelative to MLE")) +
  ggtitle("Singletons Masked")
p5a
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/grid.search/COM-noSingletons/2Epoch.RealUnits.noSingletons.pdf",p5a,height=5,width=7)
