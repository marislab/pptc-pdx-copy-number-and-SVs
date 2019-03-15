##########################################################################
# Author: Komal S Rathi & Jo Lynne Rokita
# Function: Boxplots for Indels, Breakpoints, and High Breakpoint Density
# Date: 02/14/2019
##########################################################################


if (!require("tidyr")){
  install.packages("tidyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("dplyr")){
  install.packages("dplyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("reshape2")){
  install.packages("reshape2", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("ggplot2")){
  install.packages("ggplot2", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (!require("plyr")){
  install.packages("plyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)

# Setting working directory
mainDir <- "~/copy-number-and-SVs/"
setwd(mainDir)
repo <- "~/pptc-pdx-copy-number-and-SVs/R/"
dataDir <- "~/copy-number-and-SVs/data/"

dir.create(file.path(mainDir,"sv_figure_result"))
result <- "~/copy-number-and-SVs/sv_figure_result/"

source(paste0(repo,"pubTheme.R"))

# clinical file 
clin.long <- read.delim(paste0(dataDir,'pptc-pdx-clinical-web.txt'), stringsAsFactors = F)
clin <- clin.long[,c('Model','Histology.Detailed')]

# colors
cols <- read.delim(paste0(dataDir,'2019-02-09-all-hist-colors.txt'), header = F, stringsAsFactors = F)
cols <- cols[which(cols$V1 %in% clin$Histology.Detailed),]

############# Fig 1
segdata <- read.delim(paste0(dataDir,'2019-02-10-252pdx-final-pptc-SNPRANK-noXY.seg'), stringsAsFactors = F)
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
for.chrom <- breakpoints %>% group_by(Model, chr, Histology.Detailed) %>% dplyr::summarise(n.breakpoints = n()) %>% as.data.frame()
setdiff(segdata$Sample, unique(breakpoints$Model)) # missing two models - likely no breakpoints?? ALL-59 and "NCH-WT-5

####### breakpoints overall
breakpoints <- breakpoints %>% group_by(Model, Histology.Detailed) %>% dplyr::summarise(n.breakpoints = n()) %>% as.data.frame()
breakpoints <- breakpoints %>% group_by(Histology.Detailed) %>% dplyr::mutate(median = median(n.breakpoints)) %>% as.data.frame()
to.include <- setdiff(segdata$Sample, breakpoints$Model) # add missing models
df <- data.frame(Model = to.include, 
                 Histology.Detailed = clin[which(clin$Model %in% to.include),'Histology.Detailed'], 
                 n.breakpoints = c(rep(0,length(to.include))), 
                 median = c(rep(0,length(to.include))))
breakpoints <- rbind(breakpoints, df)
breakpoints$Histology.Detailed <- reorder(breakpoints$Histology.Detailed, breakpoints$median)
write.table(breakpoints, file = paste0(result,'Breakpoints_rawdata.txt'), quote = F, sep = "\t", row.names = F)

cols <- cols[match(levels(breakpoints$Histology.Detailed), cols$V1),] # reorder colors to match histology
group.colors <- setNames(cols$V2, nm = cols$V1)

# boxplot
p <- ggplot(breakpoints, aes(x = Histology.Detailed, y = n.breakpoints, color = Histology.Detailed, alpha = 0.5)) + 
  geom_boxplot(outlier.shape = 21) + 
  geom_jitter(position=position_jitter(width=.1, height=0), shape = 21) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  guides(alpha = FALSE, fill = FALSE) + 
  theme_Publication() + xlab("Histology") + ylab('Number of Breakpoints') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  guides(color = F) +
  scale_color_manual(values = group.colors)
p
ggsave(filename = paste0(result,'BreakpointsPlot.pdf'), plot = p, device = 'pdf', height = 6, width = 10)

breakpoint.medians <- unique(breakpoints[,c("Histology.Detailed","median")])
write.table(breakpoint.medians, file = paste0(result,'Breakpoints_medians.txt'), quote = F, sep = "\t", row.names = F)

################### Fig 2
# mutations (ins/dels only)
mut <- read.delim(paste0(dataDir,'2019-02-13-mutations-per-model.txt'), stringsAsFactors = F)
mut$Histology.Detailed <- NULL
mut <- unique(mut[,1:5])
mut$value <- rowSums(mut[,2:5])
mut <- merge(mut, clin, by = "Model", all.y = T)
mut <- unique(mut[,c("Model","Histology.Detailed","value")])
mut[is.na(mut)] <- 0
mut$value <- log2(mut$value + 1)
mut <- mut %>% group_by(Histology.Detailed) %>% dplyr::mutate(median = median(value)) %>% as.data.frame()
mut$Histology.Detailed <- reorder(mut$Histology.Detailed, mut$median)
write.table(mut, file = paste0(result,'Indels_rawdata.txt'), quote = F, sep = "\t", row.names = F)

# boxplot
q <- ggplot(mut, aes(x = Histology.Detailed, y = value, color = Histology.Detailed, alpha = 0.5)) + 
  geom_boxplot(outlier.shape = 21, fill = 'white') + 
  geom_jitter(position=position_jitter(width=.1, height=0), shape = 21) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  guides(alpha = FALSE, fill = FALSE) + 
  theme_Publication() + xlab("Histology") + ylab('Number of Indels (log2)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  guides(color = F) +
  scale_color_manual(values = group.colors)
q
ggsave(filename = paste0(result,'SVplot_log2.pdf'), plot = q, device = 'pdf', height = 6, width = 13)

sv.medians <- unique(mut[,c("Histology.Detailed","median")])
write.table(sv.medians, file = paste0(result,'SVplot_medians.txt'), quote = F, sep = "\t", row.names = F)


####CHROMOTHRIPSIS######

# breakpoints by chromosome for breakpoint density analysis
for.chrom$breakpoints <- paste0(for.chrom$chr,':',for.chrom$breakpoints)
for.chrom$chr <- gsub(pattern = "chr", replacement = "", x = for.chrom$chr) 
# relevel chromosomes
for.chrom$chr <- factor(for.chrom$chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                                                  "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"))

# subset models with N >= 10 per one chromosome
chromo <- subset(for.chrom, n.breakpoints >=10)
chromo$chr.bp <- paste0("Chr", chromo$chr, ":", chromo$n.breakpoints)
length(unique(chromo$Model))/length(unique(for.chrom$Model)) # percent of models with HDB
head(chromo)

# create table for supplement
gathered <- chromo[,c("Model", "chr.bp", "Histology.Detailed")]

unite.chrom <- gathered %>% 
  group_by(Histology.Detailed, Model) %>% 
  dplyr::summarise_all(funs(trimws(paste(., collapse = ', '))))
names(unite.chrom) <- c("Histology", "Model", "Chromosomes:Breakpoints")
# export for main table
write.table(unite.chrom, paste0(result,'high-density-breakpoints.txt'), sep = "\t", quote = F, row.names = F, col.names = T)

# collect n per each platform
snp <- as.data.frame(table(clin.long$Histology.Detailed, clin.long$Have.snp.file))
names(snp) <- c("Histology.Detailed", "snp", "snp.N")
snp <- subset(snp, snp != "no")
rna <- as.data.frame(table(clin.long$Histology.Detailed, clin.long$RNA.Part.of.PPTC))
names(rna) <- c("Histology.Detailed", "rna", "rna.N")
rna <- subset(rna, rna != "no")
wes <- as.data.frame(table(clin.long$Histology.Detailed, clin.long$Have.maf))
names(wes) <- c("Histology.Detailed", "wes", "wes.N")
wes <- subset(wes, wes != "no")
numbers <- merge(merge(snp, rna), wes)
numbers$snp <- NULL
numbers$rna <- NULL
numbers$wes <- NULL
# output numbers to file
total <- as.data.frame(table(clin.long$Histology.Detailed))
write.table(numbers, paste0(result,"numbers.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(total, paste0(result,"total.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

# percent of total models with high breakpoint density 
length(unique(chromo$Model))/sum(numbers$snp.N)

# unique per histology
chromo.hist <- unique(chromo[c("Model", "Histology.Detailed")])
bpdens <- as.data.frame(table(chromo.hist$Histology.Detailed))
names(bpdens) <- c("Histology.Detailed", "n.chrom")
n.chrom <- merge(numbers, bpdens, all.x = T)
n.chrom$perc.chrom <- round(n.chrom$n.chrom/n.chrom$snp.N*100, 1)
# make NA == 0
n.chrom$perc.chrom <- ifelse(is.na(n.chrom$perc.chrom), 0, n.chrom$perc.chrom)
n.chrom$perc.no <- 100-n.chrom$perc.chrom
melt.chrom <- melt(n.chrom, id.vars = "Histology.Detailed", measure.vars = c("perc.chrom", "perc.no"))

# reorder levels by perc.chrom
levels(melt.chrom$Histology.Detailed)
newlevels <- subset(melt.chrom, variable == "perc.chrom")
melt.chrom$Histology.Detailed <- factor(melt.chrom$Histology.Detailed, levels=newlevels$Histology.Detailed[order(-newlevels$value)])

################### Fig 3
# chromothripsis figure
q <- ggplot(melt.chrom, aes(y = value, x = Histology.Detailed, fill = variable, 
                            label=paste0(round(value,0),"%"))) + 
  geom_bar(stat="identity", position = position_fill(reverse = TRUE))+
  scale_fill_manual(values = c("firebrick3", "gray89"), 
                    name="Breakpoint\nDensity",
                    labels=c("high", "low"))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme_Publication()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(x = "Histology", y = "Percent of Models")
q
ggsave(filename = paste0(result,'HighBreakpointDensity.pdf'), plot = q, device = 'pdf', height = 6, width = 10)
