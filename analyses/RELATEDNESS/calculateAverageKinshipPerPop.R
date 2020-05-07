######### get average kinship 
kinship <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/Relatedness/20181119/AllKinships.EachPopAnalyzedSeparately.CA-BAJCombined.KinshipOnly.TableForManuscript.txt",sep="\t",header=T)
head(kinship)
require(dplyr)
avgKinship <- kinship %>%
  group_by(population) %>%
  summarize(avgKinship = mean(kinship))
write.table( avgKinship, "/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/Relatedness/20181119/averageKinshipPerPopulation.txt",sep="\t",row.names = F,quote=F)
mean(avgKinship$avgKinship)
