require(ggplot2)
require(scales)
require(RColorBrewer)
require(scales)

mu=8.64e-09
model1="dadiModel1.OldShallowContraction"
model2="dadiModel2.RecentExtremeContraction"
model3="dadiModel3.BothContractions"
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/StitchPSMC_SFS_Together/CA/runMSMCOnSimluations/"

########## function to scale results #######
# input is .final.txt file, read into a df
###### info from msmc guide ons caling:
"######################## scale results #################
# to convert: from the msmc guide
# MSMC outputs times and rates scaled by the mutation rate per basepair per generation. First, scaled times are given in units of the per-generation mutation rate. This means that in order to convert scaled times to generations, divide them by the mutation rate. In humans, we used mu=1.25e-8 per basepair per generation.To convert generations into years, multiply by the generation time, for which we used 30 years.
# 
# To get population sizes out of coalescence rates, first take the inverse of the coalescence rate, scaledPopSize = 1 / lambda00. Then divide this scaled population size by 2*mu (yes, this factor 2 is different from the time scaling, sorry)."
# common plotting scale:
commonScaleMin=100
commonScaleMax=20000


scaleMSMC <- function(input,mu,gen,spp,label,category="main"){
  input$Ne <- (1/input$lambda_00)/(2*mu) # note the factor of 2! (not in time scaling) confirmed correct: https://github.com/stschiff/msmc-tools/blob/master/plot_utils.py
  input$LeftYears <- gen*(input$left_time_boundary/mu)
  input$RightYears <- gen*(input$right_time_boundary/mu)
  input$label <- label
  input$spp <- spp
  input$category <- "main" # main result (not bootstrap)
  input$Left_generations <- (input$left_time_boundary/mu)
  input$Right_generations <- (input$right_time_boundary/mu)
  return(input)
}


############################### model 1 #############################
# this is the actual demographic model
model=model1
modelRect1=data.frame(xmin=0,xmax=1000,ymin=0,ymax=2000)
modelRect2=data.frame(xmin=1000,xmax=30000,ymin=0,ymax=4500)
model_allReps <- data.frame()
for(replicate in c(seq(1,10))){
  # check if file is empty; skip if it is
  if(file.exists(paste(data.dir,model,"/rep_",replicate,"/msmc.RunOn.",model,".sims.out.final.txt",sep=""))){
    results <- read.table(paste(data.dir,model,"/rep_",replicate,"/msmc.RunOn.",model,".sims.out.final.txt",sep=""),header=T)
    results_scaled <- scaleMSMC(results,mu,gen=7,model,paste("rep_",replicate,sep=""),category="main")
# note now this result is labeled with its rep number (in the function), so can tell apart
    model_allReps = rbind(model_allReps,results_scaled)
  }
  else{ print(paste("rep",replicate," does not exist!"))}
}
p1 <- ggplot(model_allReps)+
  geom_rect(data=modelRect1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="purple",alpha=0.4)+
  geom_rect(data=modelRect2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="purple",alpha=0.4)+
  geom_step(stat="identity",size=0.5,aes(x=Left_generations,y=Ne,group=label),color="black")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(model)+
  xlab("Generations")+
  ylab("IICR") +
  scale_y_log10(labels=comma,limits=c(commonScaleMin,commonScaleMax))+
  scale_x_log10() +
  theme(text=element_text(size=14))

p1
ggsave(paste(data.dir,model,"/MSMC.RunOnSimulations.CommonScale.",model,".pdf",sep=""),p1,height=5,width=7)
############################### model 2 #############################
model=model2
modelRect1=data.frame(xmin=0,xmax=35,ymin=0,ymax=195)
modelRect2=data.frame(xmin=35,xmax=30000,ymin=0,ymax=3500)
model_allReps <- data.frame()
for(replicate in c(seq(1,10))){
  # check if file is empty; skip if it is
  if(file.exists(paste(data.dir,model,"/rep_",replicate,"/msmc.RunOn.",model,".sims.out.final.txt",sep=""))){
    results <- read.table(paste(data.dir,model,"/rep_",replicate,"/msmc.RunOn.",model,".sims.out.final.txt",sep=""),header=T)
    results_scaled <- scaleMSMC(results,mu,gen=7,model,paste("rep_",replicate,sep=""),category="main")
    # note now this result is labeled with its rep number (in the function), so can tell apart
    model_allReps = rbind(model_allReps,results_scaled)
  }
  else{ print(paste("rep",replicate," does not exist!"))}
}
p2 <- ggplot(model_allReps)+
  geom_rect(data=modelRect1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="purple",alpha=0.4)+
  geom_rect(data=modelRect2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="purple",alpha=0.4)+
  geom_step(stat="identity",size=0.5,aes(x=Left_generations,y=Ne,group=label),color="black")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(model)+
  xlab("Generations")+
  ylab("IICR") +
  scale_y_log10(labels=comma,limits=c(commonScaleMin,commonScaleMax))+
  scale_x_log10() +
  theme(text=element_text(size=14))

p2
ggsave(paste(data.dir,model,"/MSMC.RunOnSimulations.CommonScale.",model,".pdf",sep=""),p2,height=5,width=7)



############################### model 3 #############################
model=model3
modelRect1=data.frame(xmin=0,xmax=22,ymin=0,ymax=270)
modelRect2=data.frame(xmin=22,xmax=1000,ymin=0,ymax=2000)
modelRect3=data.frame(xmin=1000,xmax=30000,ymin=0,ymax=4500)

model_allReps <- data.frame()
for(replicate in c(seq(1,10))){
  # check if file is empty; skip if it is
  if(file.exists(paste(data.dir,model,"/rep_",replicate,"/msmc.RunOn.",model,".sims.out.final.txt",sep=""))){
    results <- read.table(paste(data.dir,model,"/rep_",replicate,"/msmc.RunOn.",model,".sims.out.final.txt",sep=""),header=T)
    results_scaled <- scaleMSMC(results,mu,gen=7,model,paste("rep_",replicate,sep=""),category="main")
    # note now this result is labeled with its rep number (in the function), so can tell apart
    model_allReps = rbind(model_allReps,results_scaled)
  }
  else{ print(paste("rep",replicate," does not exist!"))}
}
p3 <- ggplot(model_allReps)+
  geom_rect(data=modelRect1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="purple",alpha=0.4)+
  geom_rect(data=modelRect2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="purple",alpha=0.4)+
  geom_rect(data=modelRect3,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="purple",alpha=0.4)+
  geom_step(stat="identity",size=0.5,aes(x=Left_generations,y=Ne,group=label),color="black")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(model)+
  xlab("Generations")+
  ylab("IICR") +
  scale_y_log10(labels=comma,limits=c(commonScaleMin,commonScaleMax))+
  scale_x_log10() +
  theme(text=element_text(size=14))

p3
ggsave(paste(data.dir,model,"/MSMC.RunOnSimulations.CommonScale.",model,".pdf",sep=""),p3,height=5,width=7)




