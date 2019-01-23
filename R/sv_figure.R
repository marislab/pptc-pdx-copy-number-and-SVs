setwd('~/Projects/Maris-lab/PPTC/')

library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)

source('R/pubTheme.R')

# clinical file 
clin <- read.delim('data/2018-12-28-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
clin <- clin[,c('Model','Histology.Detailed')]
# old <- read.delim('data/2018-08-23-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
# setdiff(old$Model, clin$Model)
# setdiff(clin$Model, old$Model)
# ICb-0614EPN was changed to ICb-10614EPN
# IC-2664PNET was removed

# colors
cols <- read.delim('data/figures/2018-08-23-all-hist-colors', header = F, stringsAsFactors = F)
cols <- cols[which(cols$V1 %in% clin$Histology.Detailed),]

############# Fig 1
segdata <- read.delim('data/figures/2018-09-20-255pdx-final-pptc-SNPRANK-noXY.seg', stringsAsFactors = F)
segdata <- segdata[which(segdata$Sample != "Rh-28"),]
segdata <- segdata[which(!segdata[,"Chromosome"] %in% c("X","Y","M")),]

cutoff <- 0.152		# which represents a 10% copy number change for a log2 ratio
cutoff <- 0.2		# which represents a 10% copy number change for diploid regions

cn_breaks <- list()
for(i in unique(segdata[,1])){
  sample_id <- i
  message(i)
  segdata_i <- segdata[which(segdata[,1] == i),]
  for(chr in c(1:22)){
    segdata_i_shr <- segdata_i[which(segdata_i[,2] == chr),]
    message(paste(chr,nrow(segdata_i_shr)) )
    if(nrow(segdata_i_shr) > 1){
      brkpos <- segdata_i_shr[which( abs(segdata_i_shr[1:nrow(segdata_i_shr)-1,6] - segdata_i_shr[2:nrow(segdata_i_shr),6]) > cutoff) +1,3]
      if(length(brkpos) > 1) {
        cn_breaks[[paste(i,chr)]]<- cbind(paste("chr",rep(chr,length(brkpos)),sep=""),brkpos,rep(sample_id,length(brkpos)))
      }else if(length(brkpos) == 1){
        cn_breaks[[paste(i,chr)]]<- cbind(paste("chr",chr,sep=""),brkpos,sample_id)
      }
    }
  }
}
breakpoints <- as.data.frame(do.call(rbind, cn_breaks))
colnames(breakpoints) <- c("chr","breakpoints","sample")
breakpoints <- merge(clin, breakpoints, by.y = 'sample', by.x = 'Model')
breakpoints$breakpoints <- paste0(breakpoints$chr,':',breakpoints$breakpoints)
breakpoints <- breakpoints %>% group_by(Model, Histology.Detailed) %>% summarise(n.breakpoints = n()) %>% as.data.frame()
breakpoints <- breakpoints %>% group_by(Histology.Detailed) %>% mutate(median = median(n.breakpoints)) %>% as_data_frame()
breakpoints$Histology.Detailed <- reorder(breakpoints$Histology.Detailed, breakpoints$median)
write.table(breakpoints, file = 'results/Breakpoints_rawdata.txt', quote = F, sep = "\t", row.names = F)

cols <- cols[match(levels(breakpoints$Histology.Detailed), cols$V1),] # reorder colors to match histology
group.colors <- setNames(cols$V2, nm = cols$V1)

# boxplot
p <- ggplot(breakpoints, aes(x = Histology.Detailed, y = n.breakpoints, color = Histology.Detailed, alpha = 0.5)) + 
  geom_boxplot(outlier.shape = 21) + 
  geom_jitter(position=position_jitter(width=.1, height=0), shape = 21) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  theme_bw() +
  guides(alpha = FALSE, fill = FALSE) + 
  theme_Publication() + xlab("Histology") + ylab('Number of Breakpoints') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  guides(color = F) +
  scale_color_manual(values = group.colors)
ggsave(filename = 'results/BreakpointsPlot.pdf', plot = p, device = 'pdf', height = 6, width = 10)

breakpoint.medians <- unique(breakpoints[,c("Histology.Detailed","median")])
write.table(breakpoint.medians, file = 'results/Breakpoints_medians.txt', quote = F, sep = "\t", row.names = F)

################### Fig 2
# mutations (ins/dels only)
mut <- read.delim('data/figures/2018-12-28-mutations-per-model.txt', stringsAsFactors = F)
mut$Histology.Detailed <- NULL
mut <- unique(mut[,1:5])
mut$value <- rowSums(mut[,2:5])
mut <- merge(mut, clin, by = "Model", all.y = T)
mut <- unique(mut[,c("Model","Histology.Detailed","value")])
mut[is.na(mut)] <- 0
mut$value <- log2(mut$value + 1)
mut <- mut %>% group_by(Histology.Detailed) %>% mutate(median = median(value)) %>% as.data.frame()
mut$Histology.Detailed <- reorder(mut$Histology.Detailed, mut$median)
write.table(mut, file = 'results/Indels_rawdata.txt', quote = F, sep = "\t", row.names = F)

# boxplot
q <- ggplot(mut, aes(x = Histology.Detailed, y = value, color = Histology.Detailed, alpha = 0.5)) + 
  geom_boxplot(outlier.shape = 21, fill = 'white') + 
  geom_jitter(position=position_jitter(width=.1, height=0), shape = 21) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  theme_bw() +
  guides(alpha = FALSE, fill = FALSE) + 
  theme_Publication() + xlab("Histology") + ylab('Number of Indels (log2)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  guides(color = F) +
  scale_color_manual(values = group.colors)

# merge
ggsave(filename = 'results/SVplot_log2.pdf', plot = q, device = 'pdf', height = 6, width = 13)

sv.medians <- unique(mut[,c("Histology.Detailed","median")])
write.table(sv.medians, file = 'results/SVplot_medians.txt', quote = F, sep = "\t", row.names = F)
