#load R packages
library(gdsfmt)
library(SNPRelate)
require(ggplot2)
require(GGally)
require(ggrepel)
############### set up your colors -- keep this consistent across all plots ######
pops=c("CA","BAJ","AK","AL","COM","KUR")  # your populations

require(RColorBrewer)
display.brewer.pal("Dark2",n=8)
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
# make Baja brown
colors=list(California=colorPal[1],Baja=colorPal[7],Alaska=colorPal[2],Aleutian=colorPal[3],Commander=colorPal[4],Kuril=colorPal[5]) # your population colors
##################################################################

calldate=20181119 # date gt's were called in format YYYYMMDD
todaysdate=format(Sys.Date(),format="%Y%m%d")
# file locations:

indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds/",calldate,"/snp9a/",sep="") # this is where your snp vcf file is and where you will save your gds file (downloaded from hoffman)
plotoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/PCA/",calldate,"/snp9a/",sep="")
fileoutdir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/PCA/snp9a/",calldate,"/",sep="")
dir.create(plotoutdir,recursive = T)
dir.create(fileoutdir,recursive = T)

####### Now I'm redoing with just my "clean" snp9 gds that has all those individuals removed
# snp9a: this has removed admixed/relatives, but does contain the full COM set (am going to downsample during this script)
genofile <- snpgdsOpen(paste(indir,"snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.gds",sep="")) # this 
# 20180802: adding LD snp pruning: (1min); r2 threshold : 0.2; recommended by SNPRelate tutorial
# https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#ld-based-snp-pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
# I double checked, and it doesn't exclude any chromosomes due to any 22 chr cutoffs. cool
head(snpset)
# Get all selected snp id
snpset.id <- unlist(snpset)
head(snpset.id)

########### pop map ########
#population information
popmap = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/allSamples.Populations.PCA.txt",header=T) # this includes the RWAB samples
sample.id = as.character(popmap$Sample)
pop1_code = as.character(popmap$PrimaryPop)
pop2_code = as.character(popmap$SecondaryPop)
seq_code = as.character(popmap$sequencer)
note_code= as.character(popmap$notes) 
############################## First PCA: maf 0.06; all individuals (except low coverage which were removed before) ###########
pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only = F, maf=0.06, missing.rate=0.2)
#sink()
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc
###### must do this pc.percent after every PCA!! otherwise the xlab and ylab will be wrong.


#################### Look at FASTSTRUCTURE and relatedness and label individuals in the popMap file  with admixed, relative, outlier ... #########


#make a data.frame
tab1a <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2=factor(pop2_code)[match(pca$sample.id,sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
head(tab1a)

#plot first 2 pc coloring by primary population
p1a <- ggplot(tab1a,aes(x=EV1,y=EV2,color=pop1,shape=sequencer))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16))+
  scale_color_manual(values=unlist(colors))
p1a

#################### Look at FASTSTRUCTURE and relatedness and label individuals in the popMap file  with admixed, relative, outlier ... #########
#plot first 2 pc coloring by primary population with shape colored by note (admixed, outlier, etc.)
# add that note section to the tab dataframe: 
tab1b <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2=factor(pop2_code)[match(pca$sample.id,sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)

p1b <- ggplot(tab1b,aes(x=EV1,y=EV2,color=pop1,shape=note))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,14,2))+
  scale_color_manual(values=unlist(colors))
p1b
ggsave(paste(plotoutdir,"/PCA.inclCA.LDPruned.BAJA.note.",todaysdate,".pdf",sep=""),p1b,device="pdf",width = 8,height=5)

##################### Plot additional PCs ################
#plot pc pairs for the first four pc 
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digit=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop1,labels=lbls)

######################  Second PCA: without MAF filter (to look at impact of low-freq variants) ################################
pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only = F, maf=0.00, missing.rate=0.2)
#population information
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

tab2 <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2=factor(pop2_code)[match(pca$sample.id,sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)


p2 <- ggplot(tab2,aes(x=EV1,y=EV2,color=pop1,shape=note,label=sample.id))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,8,2))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  scale_color_manual(values=unlist(colors))
p2
ggsave(paste(plotoutdir,"/PCA.inclCA.LDPruned.BAJA.noMAF.",todaysdate,".pdf",sep=""),p2,device="pdf",width = 7,height=5)

######################  Third PCA: exclude southern sea otter populaitons (California + baja) ################################
sub <- tab1b$sample.id[tab1b$pop1!="California" & tab1b$pop1!="Baja"]
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc
# check output! 
tab3 <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)

#plot first 2 pc coloring by primary population, excluding California and Baja
p3a <- ggplot(tab3,aes(x=EV1,y=EV2,color=pop1,shape=note,label=sample.id))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,8))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  scale_color_manual(values=unlist(colors))
p3a
ggsave(paste(plotoutdir,"/PCA.excludeCA.BAJ.LDPruned.note.",todaysdate,".pdf",sep=""),p3a,device="pdf",width = 7,height=5)


# plot with names
p3b <- ggplot(tab3,aes(x=EV1,y=EV2,color=pop1,shape=note,label=sample.id))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,8))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  geom_text_repel(aes(label=sample.id))+
  scale_color_manual(values=unlist(colors))
p3b
ggsave(paste(plotoutdir,"/PCA.excludeCA.LDPruned.withNames.",todaysdate,".pdf",sep=""),p3b,device="pdf",width = 7,height=5)

############### Fourth PCA : Just Southern sea otter (CA and Baja) #########
sub <- tab1b$sample.id[(tab1b$pop1=="California" |tab1b$pop1=="Baja") & tab1b$sequencer!="HiSeq4000" & tab1b$note=="good"]
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

tab4 <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#close gds file

p4 <- ggplot(tab4,aes(x=EV1,y=EV2,color=pop1,shape=note,label=sample.id))+
  geom_point()+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(16,2))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  scale_color_manual(values=unlist(colors))+
 geom_text_repel(aes(label=sample.id))
  
p4
ggsave(paste(plotoutdir,"/PCA.CA.BAJ.only.LDPruned.",todaysdate,".pdf",sep=""),p4,device="pdf",width = 7,height=5)

################# fifth PCA : downsample Commanders; plot with and without CA/BAJ ########################
# this is a mostly random sample of commanders individuals to exclude to downsample
excludeCOM = c("132_Elut_BER_35","133_Elut_BER_37","160_Elut_MED_18","62_Elut_MED_16","46_Elut_BER_29","48_Elut_BER_33","162_Elut_MED_21","47_Elut_BER_32","49_Elut_BER_36","51_Elut_BER_39","53_Elut_BER_93","59_Elut_MED_13","157_Elut_MED_30","63_Elut_MED_17","71_Elut_BER_44","73_Elut_BER_50","82_Elut_MED_20","84_Elut_MED_23","85_Elut_MED_24","87_Elut_MED_26","88_Elut_MED_29")

sub <- tab1b$sample.id[!(tab1b$sample.id %in% excludeCOM) & tab1b$sequencer!="HiSeq4000" & tab1b$note=="good"]
#sink(paste(fileoutdir,"/PCA.record.only.CA.Baja.",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
#sink()
# check output! 
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

tab5a <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#close gds file

p5a <- ggplot(tab5a,aes(x=EV1,y=EV2,color=pop1,shape=note,label=sample.id))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,8,2))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs",sep=""))+
  scale_color_manual(values=unlist(colors))
p5a
ggsave(paste(plotoutdir,"/PCA.downsampleCOM.only.LDPruned.",todaysdate,".pdf",sep=""),p5a,device="pdf",width = 7,height=5)

################# fifth PCA (c) : downsample Commanders; plot with and without CA/BAJ/AK so just focusing on KUR,COM,AL ########################
# this is a mostly random sample of commanders individuals to exclude to downsample
# including the ones that sergio flagged as clustering oddly in his first NJ tree (though later more randomly group differently)
excludeCOM = c("132_Elut_BER_35","133_Elut_BER_37","160_Elut_MED_18","62_Elut_MED_16","46_Elut_BER_29","48_Elut_BER_33","162_Elut_MED_21","47_Elut_BER_32","49_Elut_BER_36","51_Elut_BER_39","53_Elut_BER_93","59_Elut_MED_13","157_Elut_MED_30","63_Elut_MED_17","71_Elut_BER_44","73_Elut_BER_50","82_Elut_MED_20","84_Elut_MED_23","85_Elut_MED_24","87_Elut_MED_26","88_Elut_MED_29")

sub <- tab1b$sample.id[!(tab1b$sample.id %in% excludeCOM) & tab1b$pop1!="California" & tab1b$pop1!="Baja" & tab1b$sequencer!="HiSeq4000" & tab1b$note=="good"]
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
# check output!
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

tab5b <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#close gds file

p5b <- ggplot(tab5b,aes(x=EV1,y=EV2,color=pop1,shape=note,label=sample.id))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,8,2))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs with MAF cutoff of >= 0.6",sep="")) +
  scale_color_manual(values=unlist(colors))
p5b
ggsave(paste(plotoutdir,"/PCA.downsampleCOM.noCA.BAJ.only.LDPruned.",todaysdate,".pdf",sep=""),p5b,device="pdf",width = 7,height=5)



################# fifth PCA (c) : downsample Commanders; plot with and without CA/BAJ/AK so just focusing on KUR,COM,AL ########################
# this is a mostly random sample of commanders individuals to exclude to downsample
excludeCOM = c("132_Elut_BER_35","133_Elut_BER_37","160_Elut_MED_18","62_Elut_MED_16","46_Elut_BER_29","48_Elut_BER_33","162_Elut_MED_21","47_Elut_BER_32","49_Elut_BER_36","51_Elut_BER_39","53_Elut_BER_93","59_Elut_MED_13","157_Elut_MED_30","63_Elut_MED_17","71_Elut_BER_44","73_Elut_BER_50","82_Elut_MED_20","84_Elut_MED_23","85_Elut_MED_24","87_Elut_MED_26","88_Elut_MED_29")

sub <- tab1b$sample.id[!(tab1b$sample.id %in% excludeCOM) & tab1b$pop1!="California" & tab1b$pop1!="Alaska" & tab1b$pop1!="Baja" & tab1b$sequencer!="HiSeq4000" & tab1b$note=="good"]
#sink(paste(fileoutdir,"/PCA.record.only.CA.Baja.",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
#sink()
# check output! 
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

tab5c <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#close gds file

p5c <- ggplot(tab5c,aes(x=EV1,y=EV2,color=pop2,shape=note,label=sample.id))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,8,2))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs with MAF cutoff of >= 0.6",sep="")) # +
  #scale_color_manual(values=unlist(colors))
p5c
ggsave(paste(plotoutdir,"/PCA.downsampleCOM.noCA.BAJ.AK.GOODONE.only.LDPruned.",todaysdate,".pdf",sep=""),p5c,device="pdf",width = 7,height=5)

################### Sixth PCA : just Aleutian Islands ####################
sub <- tab1b$sample.id[tab1b$pop1=="Aleutian" & tab1b$sequencer!="HiSeq4000" & tab1b$note=="good"]
#sink(paste(fileoutdir,"/PCA.record.only.CA.Baja.",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
#sink()
# check output! 
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

tab6 <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#close gds file

p6 <- ggplot(tab6,aes(x=EV1,y=EV2,color=pop2,shape=note,label=sample.id))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,8,2))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs with MAF cutoff of >= 0.6",sep="")) # +
#scale_color_manual(values=unlist(colors))
p6
ggsave(paste(plotoutdir,"/PCA.AleutianOnly.LDPruned.",todaysdate,".pdf",sep=""),p6,device="pdf",width = 7,height=5)


################### Sixth PCA : just Aleutian Islands ####################
sub <- tab1b$sample.id[tab1b$pop1=="Commander" & tab1b$sequencer!="HiSeq4000" & tab1b$note=="good"]
#sink(paste(fileoutdir,"/PCA.record.only.CA.Baja.",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
#sink()
# check output! 
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

tab7 <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#close gds file

p7<- ggplot(tab7,aes(x=EV1,y=EV2,color=pop2,shape=note,label=sample.id))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,16,8,2))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs with MAF cutoff of >= 0.6",sep="")) # +
#scale_color_manual(values=unlist(colors))
p7
ggsave(paste(plotoutdir,"/PCA.CommandersOnly.LDPruned.",todaysdate,".pdf",sep=""),p7,device="pdf",width = 7,height=5)

##################### FINAL PCA: only the "good" samples (not outliers, relatives, admixed, etc) ################
sub <- tab1b$sample.id[tab1b$note=="good"] # no other notes allowed 
#sink(paste(fileoutdir,"/PCA.record.only.CA.Baja.",todaysdate,".txt",sep=""))
pca <- snpgdsPCA(genofile, snp.id=snpset.id,autosome.only = F, sample.id=sub, maf=0.06)
#sink()
pc.percent <- pca$varprop*100
pc = head(round(pc.percent, 2))
pc

# check output! 
tab_final <- data.frame(sample.id = pca$sample.id, pop1 = factor(pop1_code)[match(pca$sample.id, sample.id)], pop2 = factor(pop2_code)[match(pca$sample.id, sample.id)],sequencer=factor(seq_code)[match(pca$sample.id, sample.id)],note=factor(note_code)[match(pca$sample.id,sample.id)],EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#close gds file

p_final <- ggplot(tab_final,aes(x=EV1,y=EV2,color=pop1,shape=note,label=sample.id))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC2 (", pc[2],"%)")) +
  xlab(paste("PC1 (", pc[1],"%)"))+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size=14))+
  scale_shape_manual(values=c(1,2,16))+
  ggtitle(paste("PCA based on ",as.character(length(pca$snp.id))," LD Pruned SNPs with MAF cutoff of >= 0.6 ",sep=""))+
  scale_color_manual(values=unlist(colors))
  
p_final
ggsave(paste(plotoutdir,"/FINAL.CLEAN.PCA.sfsQualityIndsOnly.plusBAJ.",todaysdate,".pdf",sep=""),p_final,device="pdf",width = 7,height=5)





############################ CLOSE THE GDS FILE ########################
# if you don't do this, it'll mess things up
snpgdsClose(genofile)

