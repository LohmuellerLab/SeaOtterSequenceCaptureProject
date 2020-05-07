### try with root of CA
require(RCircos) # have to install from this source file to be the exact right version: /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/BITE/BITE/BITE/extrapackages/RCircos_1.1.3.tar.gz
require(reshape2)
require(RColorBrewer)
require(GenABEL)
require(car)
require(rmarkdown)
require(GenABEL)
require(BITE)
require(ape)
source("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/TREEMIX/plotting_funcs.R")  # for plotting residuals
source("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/BITE/BITE/BITE/R/treemix_boostrap.R") # note, I modified this source code to increase axis label size using axis.cex
pops="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/poporder.noBAJA.noFerret.combineAL.txt" # for plotting residuals
# pop names in order for residual plotting
########## Plotting treemix consensus tree from my modified BITE
# recall: when I modified the BITE script I added the -bootstrap flag to treemix so it would resample my data 
#### to see how to plot, look at:
# ?treemix.bootstrap from BITE
# treemix.bootstrap(in.file, out.file = "tmp", phylip.file, 
##pop.color.file = NULL, 
#nboot, 
#cex = 1, disp = 0.003, plus = 0.01, flip = vector(), 
#arrow = 0.05, ybar = 0.1, xmin = 0, lwd = 1, font = 1, 
#scale = T, mbar = T, plotboot = T, plotmig = T, plotnames = T, 
#...)
migrations=c(seq(1,5)) # converges fast after 3 so don't need >5.
k=100
nboot=100 ## make sure this is right!!! v important!
fileFlag="snp7" # or snp9
root="CA" # other option : NoOutgroup
global="global" # or ""
alFlag="combine-AL" # or sep-AL if you want aleutian islands to be separate
##### first get mle likelihoods : 
LLs <- data.frame(mig=numeric(),ll=numeric())
for(mig in migrations){
  data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/BITE-treemix/",fileFlag,"/",alFlag,"/root.",root,"/",global,"/mig",mig,"/",sep="")
  llik <- read.table(paste(data.dir,"treemix.global.k.",k,".root.",root,".resampednBoots.",nboot,".tfConsensus.llik",sep=""))
  df=llik[,c("V4","V7")]
  colnames(df) <- c("mig","ll")
  LLs <- rbind(LLs,df)
}
LLs$AIC <- 2*LLs$mig-2*LLs$ll
# write out all:
write.table(LLs,paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/BITE-treemix/",fileFlag,"/",alFlag,"/root.",root,"/global/MigEdgesandLLs.all.txt",sep=""),row.names = F,quote=F)
# pick mig that maximizes LL:
maxLL=max(LLs$ll)
# and that minimizes AIC
minAIC=min(LLs$AIC)
# get LLs that are within 2 pts of MLE: 
outLL1 = LLs[LLs$ll >= maxLL-2 & LLs$ll <= maxLL+2,]
write.table(outLL1,paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/BITE-treemix/",fileFlag,"/",alFlag,"/root.",root,"/global/MigEdgesThatMaxLL.2ptsOfMle.txt",sep=""),row.names = F,quote=F)
# get LLs that minimize AIC (accounts for extra params)
outLL2 = LLs[LLs$AIC==minAIC,]
write.table(outLL2,paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/BITE-treemix/",fileFlag,"/",alFlag,"/root.",root,"/global/MigEdgesThatMinAIC.txt",sep=""),row.names = F,quote=F)



######## make plots ##############
for(mig in migrations){
  data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/BITE-treemix/",fileFlag,"/",alFlag,"/root.",root,"/",global,"/mig",mig,"/",sep="")
  input <- paste(data.dir,"treemix.global.k.",k,".root.",root,".resampednBoots.",nboot,".tfConsensus",sep="") 
  # add in lhood:
  newick <- paste(data.dir,"treemix.global.k.",k,".root.",root,".resampednBoots.",nboot,"_outtree.newick",sep="") # from consensus tree; branch lengths are how many boots the node showed up in (support of the node)
  pdf(paste(data.dir,"bootstrappedTree.pdf"),height=7,width=10)
  treemix.bootstrap(in.file = input,phylip.file = newick,nboot=nboot,boot.cex = 2,boot.legend.location = "bottomleft",boot.legend.cex=0.9,ybar=0.3,xbar=0.0,cex=2.5,disp = 0.001)
  dev.off()
  # also want to plot consensus tree topology:
  newickTree <- read.tree(newick)
  pdf(paste(data.dir,"consensus.Tree.Topology.WithBootValues.pdf"))
  plot(newickTree,use.edge.length=F)
  edgelabels(newickTree$edge.length,bg="transparent",frame="none",col="purple")
  title("Consensus Tree topology\nbootstrap values shown; no branch lengths")
  dev.off()
  # plot residuals:
  pdf(paste(data.dir,"redisuals.pdf"))
  plot_resid(input,pops)
  dev.off()
  

}

#### Refine mig=3  plot for manuscript: #####
mig=3
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/TREEMIX/20181119/BITE-treemix/",fileFlag,"/",alFlag,"/root.",root,"/",global,"/mig",mig,"/",sep="")
input <- paste(data.dir,"treemix.global.k.",k,".root.",root,".resampednBoots.",nboot,".tfConsensus",sep="") 
# add in lhood:
newick <- paste(data.dir,"treemix.global.k.",k,".root.",root,".resampednBoots.",nboot,"_outtree.newick",sep="") # from consensus tree; branch lengths are how many boots the node showed up in (support of the node)
pdf(paste(data.dir,"bootstrappedTree.ForManuscript.pdf"),height=7,width=10)
# taking out mig bar; 
treemix.bootstrap(in.file = input,phylip.file = newick,nboot=nboot,boot.cex = 5,boot.legend.location = "bottomleft",boot.legend.cex=2,scale=T,mbar = F,ybar=0.3,xbar=0.025,cex=2.5,disp = 0.001)
dev.off()

