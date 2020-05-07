#load R packages
library(gdsfmt)
library(SNPRelate)
require(ggplot2)
# guide: https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#principal-component-analysis-pca
# For relatedness analysis, identity-by-descent (IBD) estimation in SNPRelate can be done by either the method of moments (MoM) (Purcell et al., 2007) or maximum likelihood estimation (MLE) (Milligan, 2003; Choi et al., 2009). For both of these methods it is preffered to use a LD pruned SNP set.

calldate=20181119 # date gt's were called in format YYYYMMDD
todaysdate=format(Sys.Date(),format="%Y%m%d")
# file locations:
pops=c("CA","BAJ","AK","AL","COM","KUR")
indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/snp7/",sep="") # this is where your snp vcf file is and where you will save your gds file (downloaded from hoffman)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/Relatedness/",calldate,sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/Relatedness/",calldate,sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)

#open the gds file
genofile <- snpgdsOpen(paste(indir,"/snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.gds",sep=""))
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# LD snp pruning: (1min); r2 threshold : 0.2; recommended by SNPRelate tutorial
# https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#ld-based-snp-pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
head(snpset)
# Get all selected snp id
snpset.id <- unlist(snpset)
head(snpset.id)

#population information
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.PCA.txt",header=T)
head(popmap)
sample.id = as.character(popmap$Sample)
pop1_code = as.character(popmap$PrimaryPop)
pop2_code = as.character(popmap$SecondaryPop)

# separate by population and make sure the samples are in dataset (removed some low cov ind)
AK.id <- popmap[popmap$PrimaryPop == "Alaska" & popmap$Sample %in% samp.id,]$Sample
CA.id <- popmap[popmap$PrimaryPop == "California" & popmap$Sample %in% samp.id,]$Sample
KUR.id <- popmap[popmap$PrimaryPop == "Kuril" & popmap$Sample %in% samp.id,]$Sample
COM.id <- popmap[popmap$PrimaryPop == "Commander" & popmap$Sample %in% samp.id,]$Sample
AL.id <- popmap[popmap$PrimaryPop == "Aleutian" & popmap$Sample %in% samp.id,]$Sample
BAJ.id <- popmap[popmap$PrimaryPop == "Baja" & popmap$Sample %in% samp.id,]$Sample
########## RElationship inference using Plink IBD #########
# each population separately
# Estimate IBD coefficients
# removing monosnps (20181211)
AK.ibd <- snpgdsIBDMoM(genofile, sample.id=AK.id, snp.id=snpset.id, remove.monosnp = T,
                    maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
CA.ibd <- snpgdsIBDMoM(genofile, sample.id=CA.id, snp.id=snpset.id,remove.monosnp = T,
                       maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
KUR.ibd <- snpgdsIBDMoM(genofile, sample.id=KUR.id, snp.id=snpset.id,remove.monosnp = T,
                       maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
COM.ibd <- snpgdsIBDMoM(genofile, sample.id=COM.id, snp.id=snpset.id,remove.monosnp = T,
                       maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
AL.ibd <- snpgdsIBDMoM(genofile, sample.id=AL.id, snp.id=snpset.id,remove.monosnp = T,
                       maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
BAJ.ibd <- snpgdsIBDMoM(genofile, sample.id=BAJ.id, snp.id=snpset.id,remove.monosnp = T,
                        maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
# baja + CA:
CA.BAJ.ibd <- snpgdsIBDMoM(genofile, sample.id=c(as.character(BAJ.id),as.character(CA.id)), snp.id=snpset.id,
                           maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
# baja + CA + AK:
CA.AK.BAJ.ibd <- snpgdsIBDMoM(genofile, sample.id=c(as.character(BAJ.id),as.character(CA.id),as.character(AK.id)), snp.id=snpset.id,
                              maf=0.05, missing.rate=0.05, num.thread=2,autosome.only = F)
# Make a data.frame of them all
AK.ibd.coeff <- snpgdsIBDSelection(AK.ibd)
AK.ibd.coeff$population <- "AK"
CA.ibd.coeff <- snpgdsIBDSelection(CA.ibd)
CA.ibd.coeff$population <- "CA"
KUR.ibd.coeff <- snpgdsIBDSelection(KUR.ibd)
KUR.ibd.coeff$population <- "KUR"
COM.ibd.coeff <- snpgdsIBDSelection(COM.ibd)
COM.ibd.coeff$population <- "COM"
AL.ibd.coeff <- snpgdsIBDSelection(AL.ibd)
AL.ibd.coeff$population <- "AL"
BAJ.ibd.coeff <- snpgdsIBDSelection(BAJ.ibd)
BAJ.ibd.coeff$population <- "BAJ"
# put into one dataframe
all.ibd.coeff <- rbind(CA.ibd.coeff,AK.ibd.coeff,AL.ibd.coeff,COM.ibd.coeff,KUR.ibd.coeff,BAJ.ibd.coeff)
# color column
all.ibd.coeff$KinshipGreaterThan0.2 <- "no"
all.ibd.coeff[all.ibd.coeff$kinship > 0.2,]$KinshipGreaterThan0.2 <- "yes"
all.ibd.coeff$pop1 <- sapply(strsplit(all.ibd.coeff$ID1,"_"),"[",3)
all.ibd.coeff$pop2 <- sapply(strsplit(all.ibd.coeff$ID2,"_"),"[",3)
# fix BER/MED to COM
all.ibd.coeff[all.ibd.coeff$pop1=="BER" | all.ibd.coeff$pop1=="MED",]$pop1 <- "COM"
all.ibd.coeff[all.ibd.coeff$pop2=="BER" | all.ibd.coeff$pop2=="MED",]$pop2 <- "COM"
all.ibd.coeff$pop1==all.ibd.coeff$pop2 # should all be true now (no cross population comparisons)
p0 <- ggplot(all.ibd.coeff,aes(x=k0,y=k1,color=KinshipGreaterThan0.2))+
  geom_point()+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits=c(0,1))+
  geom_abline(intercept = 1,slope = -1,color="red",linetype="dashed")+
  facet_wrap(~population)+
  theme_bw()+
  scale_color_manual(values=c("dodgerblue","orange"))
p0
ggsave(paste(plotoutdir,"/relatedness.PlinkMoM.allPops.",todaysdate,".pdf",sep=""),p0,device="pdf",width = 8,height=5)

write.table(all.ibd.coeff[all.ibd.coeff$KinshipGreaterThan0.2=="yes",],paste(fileoutdir,"/relatedness.plinkMoM.",todaysdate,".individualsKinshipGT0.2.txt",sep=""),quote=F,row.names = F)

################### plot all kinship coefficients in a heatmap ##########
pops=c("CA" , "BAJ", "AK" , "AL" ,  "KUR")

for(pop in pops){
  # note: pop1 and pop2 should be the same, but this is defensive coding in case they aren't
  info = all.ibd.coeff[all.ibd.coeff$pop1==pop & all.ibd.coeff$pop2==pop,]
  p0b <- ggplot(info,aes(x=ID1,y=ID2,fill=kinship))+
    geom_tile()+
    scale_fill_gradientn(colours=c("white","gray","pink"))+
    geom_text(aes(label = round(kinship,2)))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste(pop," kinship coefficients from Plink\ncontains relatives/admixed/outliers"))
  p0b
  ggsave(paste(plotoutdir,"/",pop,".kinshipHeatMap.PlinkMoM.",todaysdate,".pdf",sep=""),p0b,device="pdf",width = 10,height=8)
}
# do Com separately because it's so big: 
pops="COM"
for(pop in pops){
  # note: pop1 and pop2 should be the same, but this is defensive coding in case they aren't
  info = all.ibd.coeff[all.ibd.coeff$pop1==pop & all.ibd.coeff$pop2==pop,]
  p0b <- ggplot(info,aes(x=ID1,y=ID2,fill=kinship))+
    geom_tile()+
    scale_fill_gradientn(colours=c("white","gray","pink"))+
    geom_text(aes(label = round(kinship,2)))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste(pop," kinship coefficients from Plink\ncontains relatives/admixed/outliers"))
  p0b
  ggsave(paste(plotoutdir,"/",pop,".kinshipHeatMap.PlinkMoM.",todaysdate,".pdf",sep=""),p0b,device="pdf",width = 20,height=16) # much bigger
}

############# California and Baja together ###########
CA.BAJ.ibd.coeff <- snpgdsIBDSelection(CA.BAJ.ibd)
CA.BAJ.ibd.coeff$pop1 <- sapply(strsplit(CA.BAJ.ibd.coeff$ID1,"_"),"[",3)
CA.BAJ.ibd.coeff$pop2 <- sapply(strsplit(CA.BAJ.ibd.coeff$ID2,"_"),"[",3)
CA.BAJ.ibd.coeff$population <- "CA-BAJ"
CA.BAJ.ibd.coeff$KinshipGreaterThan0.2 <- "no"
CA.BAJ.ibd.coeff[CA.BAJ.ibd.coeff$kinship > 0.2,]$KinshipGreaterThan0.2 <- "yes"

p1a <- ggplot(CA.BAJ.ibd.coeff,aes(x=k0,y=k1,color=interaction(pop1,pop2)))+
  geom_point()+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits=c(0,1))+
  geom_abline(intercept = 1,slope = -1,color="red",linetype="dashed")+
  theme_bw()
  #scale_color_manual(values=c("dodgerblue","orange"))
p1a
ggsave(paste(plotoutdir,"/relatedness.PlinkMoM.CA.Baja.",todaysdate,".pdf",sep=""),p1a,device="pdf",width = 8,height=5)


############## CA + Baja plot kinship coefficients: heat map ################
p1b <- ggplot(CA.BAJ.ibd.coeff,aes(x=ID1,y=ID2,fill=kinship))+
  geom_tile()+
  scale_fill_gradientn(colours=c("white","gray","pink"))+
  geom_text(aes(label = round(kinship,2)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("CA-BAJ Kinship coefficients from Plink")
p1b
ggsave(paste(plotoutdir,"/CA.BAJ.kinshipHeatmap.PlinkMoM.CA.Baja.",todaysdate,".pdf",sep=""),p1b,device="pdf",width = 10,height=8)

################ WRITE OUT A TABLE WITH EACH POPULATION and WITH CA+BAJ #############
allPopsExceptCABAJ <- all.ibd.coeff[!(all.ibd.coeff$population %in% c("CA","BAJ")),]
# combine with the CA-BAJ combined analysis:
allPops_CABAJCombo <- rbind(allPopsExceptCABAJ,CA.BAJ.ibd.coeff)
write.table(allPops_CABAJCombo,paste(fileoutdir,"/AllKinships.EachPopAnalyzedSeparately.CA-BAJCombined.FullTable.txt",sep=""),row.names = F,col.names = T,sep="\t",quote=F)
write.table(allPops_CABAJCombo[,c("population","ID1","ID2","kinship")],paste(fileoutdir,"/AllKinships.EachPopAnalyzedSeparately.CA-BAJCombined.KinshipOnly.TableForManuscript.txt",sep=""),row.names = F,col.names = T,sep="\t",quote=F)

