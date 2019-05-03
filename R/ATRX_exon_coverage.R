###################################
# Author: Komal S Rathi
# Function: ATRX Deletion Analysis
# Date: 08/27/2018
###################################

setwd('~/Projects/Maris-lab/ATRX/')

library(ggplot2)
library(tidyr)
library(dplyr)

lf <- list.files(pattern = '*_coverage.txt', path = 'data', full.names = T)
for(i in 1:length(lf)){
  print(i)
  fname <- gsub('data/|_coverage.txt','',lf[i])
  if(i == 1){
    dat <- read.delim(lf[i], stringsAsFactors = F, header = F)
    dat$Sample <- fname
  } else {
    x <- read.delim(lf[i], stringsAsFactors = F, header = F)
    x$Sample <- fname
    dat <- rbind(dat, x)
  }
}
dat$exon_length <- dat$V3-dat$V2
dat <- dat[,c(4,7,8,9)]
colnames(dat) <- c("id","read_count","Sample","exon_length")
dat$transcript <- gsub('_.*','',dat$id)
dat$exon_n <- gsub('.*_','',gsub('_ENSE.*','',dat$id))
dat$exon_name <- gsub('.*_','',dat$id)

# add lib sizes
lib.size <- read.delim('data/flagstat_output.txt', header = F)
dat <- merge(dat, lib.size, by.x = 'Sample', by.y = 'V1')
colnames(dat)[8] <- "library_size"

# lib.size <- 79125352
dat$per.million.scaling.factor <- dat$library_size/10^6
dat$FPKM <- ((dat$read_count+1)/dat$per.million.scaling.factor)/dat$exon_length
res <- dat[,c(1,2,5,6,7,8,3,4,9,10)]
# upper <- (dat$read_count+1)*(10^6)
# lower <- (dat$exon_length)*(dat$library_size)
# dat$FPKM <- upper/lower

# reorder boxes
dat <- dat %>% group_by(Sample) %>% mutate(median = median(FPKM))
dat$Sample <- reorder(dat$Sample, dat$median)

# plot exon coverage boxplot
ggplot(dat, aes(x = Sample, y = FPKM)) + geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) + 
  ggtitle('Distribution of Exon Coverage (ATRX)')
ggsave(filename = 'ATRX_exon_coverage.pdf', width = 30, height = 10)  
