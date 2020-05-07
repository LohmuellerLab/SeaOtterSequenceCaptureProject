# Grid search
# also want to input fsc parameters
require(ggplot2)
populations=c("CA","AK","KUR","AL")
#populations="CA"
genotypeDate="20181119"
model="1D.2Epoch"
mu=8.64e-09

# modification from Kirk: plot relative to MLE LL rather than to data:data LL
############################### set up vectors of inferred paramaters as reference points ###########
### from dadi and fsc results, set reference points (figure out good way to make this standardized)
# draw from: /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/20181119/excelResultsComparisons/20190110/excelResultsSummary.dadi.FSCSummaryAdded.AllPops.20190110.xlsx
# set up a dictionary-like list of with the MLE results of dadi and fsc inference from above file
# contraction sizes in diploids
######## eventually write out this out as a table in a different script to automate #############
# total sites assessed: (these values are also in excel file)
Ls = vector(mode="list",length=length(populations))
Ls["CA"] <- 5989967
Ls["AK"] <- 6335196
Ls["AL"] <- 6379260
Ls["KUR"] <- 6416481

dadiNus = vector(mode="list",length=length(populations))
names(dadiNus) <- populations
dadiNus["CA"] <- 0.8
dadiNus["AK"] <- 29.5
dadiNus["AL"] <- 114.2
dadiNus["KUR"] <- 1.7

# sizes in DIPLOIDS
fscNus = vector(mode="list",length=length(populations))
names(fscNus) <- populations
fscNus["CA"] <- 309.5
fscNus["AK"] <- 389.5
fscNus["AL"] <- 1811
fscNus["KUR"] <- 1054.5

# times in generations 
dadiTs = vector(mode="list",length=length(populations))
names(dadiTs) <- populations
dadiTs["CA"] <- 0.2
dadiTs["AK"] <- 3.7
dadiTs["AL"] <- 10.2
dadiTs["KUR"] <- 0.1

fscTs = vector(mode="list",length=length(populations))
names(fscTs) <- populations
fscTs["CA"] <- 73
fscTs["AK"] <- 56
fscTs["AL"] <- 306
fscTs["KUR"] <- 80

dadiNancs = vector(mode="list",length=length(populations))
names(dadiNancs) <- populations
dadiNancs["CA"] <- 3585.3
dadiNancs["AK"] <- 4381.4
dadiNancs["AL"] <- 4834.3
dadiNancs["KUR"] <- 4348.6


fscNancs = vector(mode="list",length=length(populations))
names(fscNancs) <- populations
fscNancs["CA"] <- 3576.5
fscNancs["AK"] <- 4411.5
fscNancs["AL"] <- 4869.5
fscNancs["KUR"] <- 4289.5


############################### read in results and plot #####################

for(pop in populations){
  indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/",genotypeDate,"/grid.search/",pop,"/",sep="")
  input <- read.table(paste(indir,"dadi.grid.search.",pop,".",model,".LL.output.txt.gz",sep=""),header=T,sep="\t")
  #input$deltaLL = input$LL_model - input$LL_data # this is relative to data:data
  # INSTEAD want to do it relative to MLE 
  input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
    # scale results by Nanc that is inferred by dadi, rather than the given Nanc 
  # get the Nanc from theta by dividing by 4*mu*L
  # and then rescale based on that (nu * Nanc_from_theta) and (T*2*Nanc_from_theta)
  input$Nanc_from_theta <- input$theta / (4*mu*Ls[[pop]])
  input$T_scaledByNanc_from_Theta_gen <- input$T * 2* input$Nanc_from_theta
  input$nu_scaledByNanc_from_Theta_dip <- input$nu * input$Nanc_from_theta
  # get grid.search MLE coordinates: 
  MLE_nu=input[input$LL_model==max(input$LL_model),]$nu
  MLE_T=input[input$LL_model==max(input$LL_model),]$T
  MLE_Nanc=input[input$LL_model==max(input$LL_model),]$Nanc_from_theta
  # PLOT LL surface #
  p1 <- ggplot(input,aes(nu,T,fill=deltaLL))+
    geom_tile()+
    scale_x_log10()+
    scale_y_log10()+
    #scale_fill_gradientn(breaks=c(-50,-25,-10,-5),limits=c(-50,0),colors =c("darkgray","gray","darkorange","yellow","green"))+
    scale_fill_gradientn(breaks=c(30,10,1),limits=c(0,30),colors =c("green","yellow","darkorange","gray","darkgray"))+
    geom_point(aes(x=dadiNus[[pop]]/dadiNancs[[pop]],y=dadiTs[[pop]]/(2*dadiNancs[[pop]])),shape=0,size=6,color="red")+
    geom_point(aes(x=fscNus[[pop]]/fscNancs[[pop]],y=fscTs[[pop]]/(2*fscNancs[[pop]])),shape=2,size=6,color="red")+
    geom_point(aes(x=MLE_nu,y=MLE_T),shape=8,size=6,color="red")+
    #annotate(geom="text",x=(MLE_nu),y=(MLE_T+0.0004),label=paste("Grid-Search MLE"),color="black")+
    #ggtitle(paste(pop, " grid Search of parameters in dadi\nRed star is grid search MLE\nRed square indicates original dadi inferred MLE params\nRed triangle is fsc2 inferred MLE parameters",sep="")) +
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20),legend.title=element_text(size=16),legend.text = element_text(size=16))+
    labs(fill = element_text("deltaLL\nrelative to MLE")) #+
    #geom_vline(xintercept = dadiNancs[[pop]],linetype="dashed")
  p1
  ggsave(paste(indir,pop,".nu.Vs.T.dadiUnits.gridSearch.pdf",sep=""),p1,height=4,width=7)
  ggsave(paste(indir,pop,".nu.Vs.T.dadiUnits.gridSearch.png",sep=""),p1,device="png",height=4,width=7)
  
  # PLOT nu and T scaled by the inferred Nancs from each dadi model#
  p3 <- ggplot(input,aes(nu_scaledByNanc_from_Theta_dip,T_scaledByNanc_from_Theta_gen,color=deltaLL))+
    geom_point(size=3)+
    scale_x_log10()+
    scale_y_log10(breaks=c(1,10,100,1000,10000))+
    scale_color_gradientn(breaks=c(30,10,1),limits=c(0,30),colors =c("green","yellow","darkorange","gray","darkgray"))+
    geom_point(aes(x=dadiNus[[pop]],y=dadiTs[[pop]]),shape=0,size=6,color="red")+
    geom_point(aes(x=fscNus[[pop]],y=fscTs[[pop]]),shape=2,size=6,color="red")+
    geom_point(aes(x=MLE_nu*MLE_Nanc,y=MLE_T*2*MLE_Nanc),shape=8,size=6,color="red")+
    annotate(geom="text",x=(MLE_nu*MLE_Nanc),y=(MLE_T*2*MLE_Nanc+20),label=paste("Grid-Search MLE\nnu: ",round(MLE_nu*MLE_Nanc,2),"\nT: ",round(MLE_T*2*MLE_Nanc,2)),color="black")+
    #annotate(geom="text",x=dadiNus[[pop]],y=dadiTs[[pop]],label=paste("Dadi Inference MLE\nnu: ",round(dadiNus[[pop]]),"\nT: ",round(dadiTs[[pop]])),color="black",size=4)+
    ggtitle(paste(pop, " grid Search of parameters in dadi\nRed star is grid search MLE\nRed square indicates original dadi inferred MLE params\nRed triangle is fsc2 inferred MLE parameters",sep="")) +
    #geom_vline(xintercept = dadiNancs[[pop]],linetype="dashed")+
    xlab("nu: Diploid invidiuals (scaled by inferred Nancs)")+
    ylab("T: Generations")+
    labs(fill = element_text("deltaLL\nrelative to MLE")) #+
  
  p3
  ggsave(paste(indir,pop,".nu.Vs.T.ScaledByInferredNanc.gridSearch.pdf",sep=""),p3,height=5,width=7)

  
}



########################## also plot BER and MED (skip fsc and dadi marks) #################
populations=c("BER","MED")
Ls["MED"] <- 6325257 # from /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/splitCOM/projection-20190301-hetFilter-0.75/dadi-plusMonomorphic/MED-21.totalSiteCount.L.withMonomorphic.txt
Ls["BER"] <- 6408246
genotypeDate="20181119"
model="1D.2Epoch"
mu=8.64e-09

for(pop in populations){
  indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/dadi_inference/",genotypeDate,"/grid.search/",pop,"/",sep="")
  input <- read.table(paste(indir,"dadi.grid.search.",pop,".",model,".LL.output.txt.gz",sep=""),header=T,sep="\t")
  #input$deltaLL = input$LL_model - input$LL_data # this is relative to data:data
  # INSTEAD want to do it relative to MLE 
  input$deltaLL = abs(input$LL_model - max(input$LL_model)) # this is relative to MLE
  # scale results by Nanc that is inferred by dadi, rather than the given Nanc 
  # get the Nanc from theta by dividing by 4*mu*L
  # and then rescale based on that (nu * Nanc_from_theta) and (T*2*Nanc_from_theta)
  input$Nanc_from_theta <- input$theta / (4*mu*Ls[[pop]])
  input$T_scaledByNanc_from_Theta_gen <- input$T * 2* input$Nanc_from_theta
  input$nu_scaledByNanc_from_Theta_dip <- input$nu * input$Nanc_from_theta
  # get grid.search MLE coordinates: 
  MLE_nu=input[input$LL_model==max(input$LL_model),]$nu
  MLE_T=input[input$LL_model==max(input$LL_model),]$T
  MLE_Nanc=input[input$LL_model==max(input$LL_model),]$Nanc_from_theta
  # PLOT LL surface #
  p1 <- ggplot(input,aes(nu,T,fill=deltaLL))+
    geom_tile()+
    scale_x_log10()+
    scale_y_log10()+
    #scale_fill_gradientn(breaks=c(-50,-25,-10,-5),limits=c(-50,0),colors =c("darkgray","gray","darkorange","yellow","green"))+
    scale_fill_gradientn(breaks=c(30,10,1),limits=c(0,30),colors =c("green","yellow","darkorange","gray","darkgray"))+
    geom_point(aes(x=MLE_nu,y=MLE_T),shape=8,size=6,color="red")+
    #annotate(geom="text",x=(MLE_nu),y=(MLE_T+0.0004),label=paste("Grid-Search MLE"),color="black")+
    labs(fill = element_text("deltaLL\nrelative to MLE")) #+
  #geom_vline(xintercept = dadiNancs[[pop]],linetype="dashed")
  p1
  ggsave(paste(indir,pop,".nu.Vs.T.dadiUnits.gridSearch.pdf",sep=""),p1,height=5,width=7)
  # PLOT nu and T scaled by the inferred Nancs from each dadi model#
  p3 <- ggplot(input,aes(nu_scaledByNanc_from_Theta_dip,T_scaledByNanc_from_Theta_gen,color=deltaLL))+
    geom_point(size=3)+
    scale_x_log10()+
    scale_y_log10(breaks=c(1,10,100,1000,10000))+
    scale_color_gradientn(breaks=c(30,10,1),limits=c(0,30),colors =c("green","yellow","darkorange","gray","darkgray"))+
    #annotate(geom="text",x=dadiNus[[pop]],y=dadiTs[[pop]],label=paste("Dadi Inference MLE\nnu: ",round(dadiNus[[pop]]),"\nT: ",round(dadiTs[[pop]])),color="black",size=4)+
    ggtitle(paste(pop, " grid Search of parameters in dadi\nRed star is grid search MLE\nRed square indicates original dadi inferred MLE params\nRed triangle is fsc2 inferred MLE parameters",sep="")) +
    #geom_vline(xintercept = dadiNancs[[pop]],linetype="dashed")+
    xlab("nu: Diploid invidiuals (scaled by inferred Nancs)")+
    ylab("T: Generations")+
    labs(fill = element_text("deltaLL\nrelative to MLE")) #+
  
  p3
  ggsave(paste(indir,pop,".nu.Vs.T.ScaledByInferredNanc.gridSearch.pdf",sep=""),p3,height=5,width=7)
}

