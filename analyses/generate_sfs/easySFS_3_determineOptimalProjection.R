require(dplyr)
require(ggplot2)
require(purrr) # for map
require(reshape2)
require(readr)
require(tidyr)
## these files are the result of easySFS_2_parseProjection
genotypeDate="20181119"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/projection_preview",sep="")
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/",genotypeDate,"/easySFS_projection/projection_preview/",sep="")
dir.create(plot.dir,recursive = T)
# get list of files:
fileList=list.files(pattern=paste("R.format",sep=""),path = data.dir,full.names = F )
fileList
data1 <- data_frame(filename = fileList) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_delim(delim=",",file.path(data.dir, .)))# a new data column 
         
  )  


data <- unnest(data1)
data$population <-  sapply(strsplit(as.character(data$filename),"\\."),'[',1)

maxima <- data %>%
  group_by(population) %>%
  filter(snps==max(snps))



p0 <- ggplot(data,aes(x=projection,y=snps))+
  geom_point()+
  facet_wrap(~population,scales="free_x")+
  theme_bw()+
  scale_x_continuous(breaks=seq(2,80,2))
p0
# want to find the maximum # of snps that still has sample size >10
ggsave(paste(plot.dir,"/easySFS.projections.allPops.pdf",sep=""),p0,height=5,width=10)

# just commanders:
p1 <- ggplot(data[data$population=="COM",],aes(x=projection,y=snps))+
  geom_point()+
  facet_wrap(~population,scales="free_x")+
  theme_bw()+
  scale_x_continuous(breaks=seq(2,80,2))
p1
ggsave(paste(plot.dir,"/easySFS.projections.COM.pdf",sep=""),p0,height=5,width=10)
