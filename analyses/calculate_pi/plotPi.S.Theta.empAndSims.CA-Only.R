require(ggplot2)
require(reshape2)
require(dplyr)
require(scales)
######### Plot pi, theta and S for empirical and simulated data
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PI_THETA/"

out.dir=paste(data.dir,"/plots/",sep="")
dir.create(out.dir,recursive = T)

models="CA.1D.2Epoch.35Gen.200Inds/20200224/"

genotypeDate="20181119" # for empirical
# input of all populations
emp <- read.table(paste(data.dir,"/empirical/",genotypeDate,"/pi.S.Theta.CalculatedFromSFSes.allpops.txt",sep=""),header=T)

# for each model read in simulated data:
simdata <- data.frame()
for(model in models){
  input <- read.table(paste(data.dir,"/simulated/",model,"/pi.S.Theta.CalculatedFromSFSes.allreps.txt",sep=""),header=T)
  # combine each model together:
  simdata <- rbind(simdata,input)
}

# make colnames match:
simdata$population <- unlist(lapply(strsplit(as.character(simdata$model),"\\."),"[",1))
simdata$label <- simdata$state
simdata$color <- "simulated"
emp$model <- "empirical"
emp$label <- "empirical"
emp$color <- "empirical"
emp$rep <- "empirical"
emp$state <- "postContraction"
# restric temp to just AK and CA:
emp <- emp[emp$population %in% c("CA"),]
emp_melt <- melt(emp)

allData <- rbind(simdata,emp)
allData_melt <- melt(allData)
#################### Plot empirical and simulated data together ##########
require(ggplot2)
# exclude S, doesnt make much sense since number of called sites differ
# get an average across replicates 
simdata_avgAcrossReps <- simdata %>%
  group_by(model,state,population,label,color) %>%
  summarise(pi=mean(pi),Wattersons_theta=mean(Wattersons_theta),sdPi=sd(pi),sdTheta=sd(Wattersons_theta))

simdata_avgAcrossReps_melt <- melt(simdata_avgAcrossReps)
simdata_avgAcrossReps_melt$label <- factor(simdata_avgAcrossReps_melt$label,levels=c("preContraction","postContraction"))
p1 <- ggplot(simdata_avgAcrossReps_melt[simdata_avgAcrossReps_melt$variable %in% c("pi","Wattersons_theta"),],aes(x=label,y=value,fill=color))+
  geom_col()+
  geom_col(data=emp_melt[emp_melt$variable %in% c("pi","Wattersons_theta"),],aes(x=label,y=value))+
  facet_wrap(population~variable,scales="free")+
  theme_bw()
p1
ggsave(paste(out.dir,"NeutralPiFromSimulated.vsEmpirical.AK.CA.pdf",sep=""),p1,height=5,width=9)

simdata_melt <- melt(simdata)
simdata_melt$label <- factor(simdata_melt$label,levels=c("preContraction","postContraction"))
simdata_melt$variable2 <- as.character(simdata_melt$variable)
simdata_melt[simdata_melt$variable=="Wattersons_theta",]$variable2 <- "Watterson's \u03b8"
simdata_melt[simdata_melt$variable=="pi",]$variable2 <- "\u03c0"
simdata_melt$variable2 <- factor(simdata_melt$variable2,levels=c("\u03c0","Watterson's \u03b8"))

emp_melt$variable2 <- as.character(emp_melt$variable)
emp_melt[emp_melt$variable=="Wattersons_theta",]$variable2 <- "Watterson's \u03b8"
emp_melt[emp_melt$variable=="pi",]$variable2 <- "\u03c0"
emp_melt$variable2 <- factor(emp_melt$variable2,levels=c("\u03c0","Watterson's \u03b8"))

p2 <- ggplot(simdata_melt[simdata_melt$variable %in% c("pi","Wattersons_theta"),],aes(x=label,y=value,fill=label))+
  geom_boxplot(alpha=0.75)+
  geom_point()+
  geom_hline(data=emp_melt[emp_melt$variable %in% c("pi","Wattersons_theta"),],aes(yintercept=value),linetype="solid",size=1.5,color="darkgray")+
  facet_wrap(population~variable2)+
  theme_bw()+
  scale_fill_manual(values=c("blue","darkred"))+
  theme(legend.position = "none",text=element_text(size=18),strip.background = element_rect(fill="white"))+
  xlab("")+
  scale_y_continuous(labels = scales::scientific)
  
p2
ggsave(paste(out.dir,"NeutralPiFromSimulated.vsEmpirical.CA.Only.NiceBoxPlot.FORSI.GrayEmpiricalLine.usethis.pdf",sep=""),p2,height=4,width=8,device=cairo_pdf)
