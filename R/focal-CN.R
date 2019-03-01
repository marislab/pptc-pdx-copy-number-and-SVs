#source("https://bioconductor.org/biocLite.R")
#biocLite("Homo.sapiens")

#install.packages("rmatio")
require(rmatio)
require(BSgenome.Hsapiens.UCSC.hg19)
require(dplyr)
require(tidyr)
library(IRanges)
library(reshape2)
library(data.table)

#set directories for saving files, specify histology of interest
###create directories for saving files
cnDir <- "~/Box Sync/PPTC-genomics-collaboration/Manuscript/scripts/focal-cn/"
pptc.folder <- "~/Box Sync/PPTC-genomics-collaboration/"
script.folder <- "~/Box Sync/PPTC-genomics-collaboration/Manuscript/scripts/oncoprint-r-scripts/"

###load RNA expression matrix
load("~/Box Sync/PPTC-genomics-collaboration/Pedcbio-upload/2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda", verbose = T) 
###load gistic output
gistic.out <- read.delim(paste0(pptc.folder, "Data/GISTIC-results/all-pdx/2018-08-09-gistic-results-256pdx-noXY-snpfast2-nomirna/all_thresholded.by_genes.txt"),as.is=TRUE,check.names=FALSE)
colnames(gistic.out)[colnames(gistic.out) == "IC-2264PNET"] <- "IC-2664PNET"

clin <- read.delim(paste0(pptc.folder, "Data/clinical/2019-02-09-pdx-clinical-final-for-paper.txt"), as.is = T, header = T)

###read in hugo file with ENS ids
hugo <- read.delim(paste0(pptc.folder,"/Data/Hugo_Symbols/2019-02-14-Hugo-Symbols-approved.txt"),
                   sep = "\t", header = T, as.is = T)
###read in gtf file, skip 5 rows of header
gtf <- read.table(paste0(pptc.folder,"Data/Hugo_Symbols/Homo_sapiens.GRCh37.87.gtf"),
                  sep = "\t", header = F, as.is = T)[-5,]


##rename columns
colnames(gtf) <- c("chr", "source", "feature", "start", "end", "NA", "strand", "frame", "attribute")
###reshape GTF and merge with hugo based on ENS IDs
gtf.genes <- subset(gtf, feature == "gene")
#remove colnames
gtf.genes$attribute <- gsub('gene_id |gene_version |gene_name |gene_source |gene_biotype ', "", gtf.genes$attribute)  
##split attribute column
gtf.genes.split<- gtf.genes %>% separate(attribute, sep = ";", 
                                         into = c("Ensembl.gene.ID", "gene_ver", "gene_symbol", "gene_source", "biotype"), remove = T)

#head(gtf.genes.split)
gtf.hugo <- merge(gtf.genes.split, hugo, by = "Ensembl.gene.ID")
gtf.hugo = gtf.hugo %>%
  mutate(chr=paste0('chr', chr))

#head(gtf.hugo)
#names(gtf.hugo)
gtf.short <- gtf.hugo[,c("chr", "start", "end", "Approved.symbol")]
gtf.short <- subset(gtf.short, chr != "chrX" & chr != "chrY")
#write.table(gtf.hugo[,c(2,5,6,1,3,4,8:ncol(gtf.hugo))], paste0("~/Box Sync/PPTC-genomics-collaboration/Data/Hugo_Symbols/", Sys.Date(), "-gtf-hugo.txt"),
 #                                                             sep = "\t", quote = F, col.names = T, row.names = F)

###load seg file and use LRR for focal CN
seg <- read.delim(paste0(pptc.folder, "Pedcbio-upload/2019-02-10-252pdx-final-pptc-SNPRANK-noXY.seg"), check.names = F) 
seg.noxy <- subset(seg, Chromosome != "Y" & Chromosome != "X")
seg2 <- seg.noxy[,c(2:ncol(seg.noxy), 1)]
names(seg2) <- c("chr", "start", "end", "markers", "CN", "Model")
seg2 = seg2 %>%
  mutate(chr=paste0('chr', chr))

#make genomic ranges object
seg.gr = seg2 %>%
  mutate(chr=paste0('chr', chr)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = F)

hug.gr = gtf.short %>%
  mutate(chr=paste0('chr', chr)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = F)

overlaps <- mergeByOverlaps(seg.gr, hug.gr)
#remove redundant chr
overlaps$seg.gr <- gsub("^chr", "", overlaps$seg.gr)
overlaps$hug.gr <- gsub("^chr", "", overlaps$hug.gr)
overlaps$markers <- NULL
#head(overlaps)
colnames(overlaps) <- c("seg.int", "CN", "Model", "gene.int", "Hugo.Symbol")
#write.table(overlaps, paste0(pptc.folder, "Manuscript/scripts/focal-cn/", Sys.Date(), "-seg-gene-model-overlaps.txt"), col.names = T, quote = F, row.names = F, sep = "\t")

ov.df <- as.data.frame(overlaps)
cn.short <- ov.df[,c("Hugo.Symbol", "Model", "CN")]

###reformat DF to wide format to create CN matrix
data_wide <- acast(cn.short, Hugo.Symbol ~ Model, value.var = "CN", fun.aggregate = mean, drop = T)
#head(data_wide)
#write CN matrix
write.table(data_wide, paste0(pptc.folder, "Manuscript/scripts/focal-cn/", Sys.Date(), "-seg-CN-matrix.txt"), col.names = T, quote = F, row.names = T, sep = "\t")

###recode CN for focal amps, dels, and no change for pedcbio/oncoprints
genecn2 <- apply(data_wide, 2, function(x) ifelse(x>=0.538578182,"Amplification", 
                                                ifelse(x<=(-1.739),"Hom_Deletion",   
                                                       ifelse(x<0.538578182&x>(-1.739),"",x))) ) 

###need to add modifications and hemizygous deletions
###add C19MC amps to ICb-1343ENB due to manual inspection
if("ICb-1343ENB" %in% colnames(genecn2)) {
  genecn2["MIR519C","ICb-1343ENB"] <- "Amplification"
  genecn2["MIR515-1","ICb-1343ENB"] <- "Amplification"
  genecn2["MIR515-2","ICb-1343ENB"] <- "Amplification"
  genecn2["LAIR1","ICb-1343ENB"] <- "Amplification"
}

###add SMARCB1 deletion in 6753ATRT due to manual inspection
if("ICb-6753ATRT" %in% colnames(genecn2)) {
  genecn2["SMARCB1","ICb-6753ATRT"] <- "Hom_Deletion"
}
###add SMARCB1 hemi-deletion due to manual inspection
if("KT-12" %in% colnames(genecn2)) {
  genecn2["SMARCB1","KT-12"] <- "Hem_Deletion"
}
if("ICb-10593ATRT" %in% colnames(genecn2)) {
  genecn2["SMARCB1","ICb-10593ATRT"] <- "Hem_Deletion"
}
if("BT-29" %in% colnames(genecn2)) {
  genecn2["SMARCB1","BT-29"] <- "Hem_Deletion"
}

###fix CDKN2B dels 
if("ALL-50" %in% colnames(genecn2)) {
  genecn2["CDKN2B","ALL-50"] <- "Hom_Deletion"
}
if("ALL-30" %in% colnames(genecn2)) {
  genecn2["CDKN2B","ALL-30"] <- "Hem_Deletion"
}
###fix PALKTY - should be hemi del for both CDKN2B
if("PALKTY" %in% colnames(genecn2)) {
  genecn2["CDKN2B","PALKTY"] <- "Hem_Deletion"
}
if("ALL-19" %in% colnames(genecn2)) {
  genecn2["CDKN2B","ALL-19"] <- "Hem_Deletion"
}
if("ALL-84" %in% colnames(genecn2)) {
  genecn2["CDKN2B","ALL-84"] <- "Hem_Deletion"
  genecn2["CDKN2A","ALL-84"] <- "Hem_Deletion"
}
if("ALL-26" %in% colnames(genecn2)) {
  genecn2["CDKN2B","ALL-26"] <- "Hom_Deletion"
}
if("ALL-29" %in% colnames(genecn2)) {
  genecn2["CDKN2B","ALL-29"] <- "Hem_Deletion"
}
if("PAKHZT" %in% colnames(genecn2)) {
  genecn2["CDKN2B","PAKHZT"] <- "Hom_Deletion"
}
if("IC-6634GBM" %in% colnames(genecn2)) {
  genecn2["CDKN2A","IC-6634GBM"] <- "Hom_Deletion"
}
if("TC-71" %in% colnames(genecn2)) {
  genecn2["CDKN2B","TC-71"] <- "Hem_Deletion"
}
###NOTE: COG-N-557x and 471x are MYCN-amp by pathology, but the amp segment starts after MYCN and was not called
if("COG-N-557x" %in% colnames(genecn2)) {
  genecn2["MYCN","COG-N-557x"] <- "Amplification"
}
if("COG-N-471x" %in% colnames(genecn2)) {
  genecn2["MYCN","COG-N-471x"] <- "Amplification"
}
###add WT1 deletion due to manual inspection
if("KT-13" %in% colnames(genecn2)) {
  genecn2["WT1","KT-13"] <- "Hem_Deletion"
}
if("NCH-WT-4" %in% colnames(genecn2)) {
  genecn2["WT1","NCH-WT-4"] <- "Hem_Deletion"
}
if("NCH-WT-7" %in% colnames(genecn2)) {
  genecn2["WT1","NCH-WT-7"] <- "Hem_Deletion"
}
if("KT-6" %in% colnames(genecn2)) {
  genecn2["WT1","KT-6"] <- "Hem_Deletion"
}  
if("KT-18" %in% colnames(genecn2)) {
  genecn2["WT1","KT-18"] <- "Hem_Deletion"
}  
if("KT-11" %in% colnames(genecn2)) {
  genecn2["WT1","KT-11"] <- "Hem_Deletion"
}
if("NCH-WT-5" %in% colnames(genecn2)) {
  genecn2["WT1","NCH-WT-5"] <- "Hem_Deletion"
}
if("NCH-WT-6-S13-1506" %in% colnames(genecn2)) {
  genecn2["WT1","NCH-WT-6-S13-1506"] <- "Hem_Deletion"
}

###osteo TP53
if("NCH-OS-1" %in% colnames(genecn2)) {
  genecn2["TP53","NCH-OS-1"] <- "Hem_Deletion"
}
if("OS-45-TSX-pr1" %in% colnames(genecn2)) {
  genecn2["TP53","OS-45-TSX-pr1"] <- "Hom_Deletion"
}
if("OS-36-SJ" %in% colnames(genecn2)) {
  genecn2["TP53","OS-36-SJ"] <- "Hom_Deletion"
}


#####gather RNA evidence of homozygous deletions
df <- setDT(as.data.frame(genecn2), keep.rownames = T)
df[df ==""] <- NA

# Convert row names into first column
setDT(df, keep.rownames = T)

# Melting nbl dataframe
melt.df <- melt(df, id = "rn")
# renaming column names
colnames(melt.df) <- c("gene_short_name","variable","value")

tmp <- copy(rna.mat)

# melt tmp -> melted_dt
setDT(tmp)
# to turn off scientific notation
options(scipen = 999)
melted_dt <- melt(tmp, id = "gene_short_name")

## creating new column for merging
melted_dt <- setDT(melted_dt)[,concat_col:=paste0(gene_short_name,",",variable)]
melt.df <- setDT(melt.df)[,concat_col:=paste0(gene_short_name,",",variable)]

# merging data tables on concat_col
new <- merge(melt.df,melted_dt, by = "concat_col", all.x=T)

# Subsetting required columns from new
new.subset <- subset(new, select = c("gene_short_name.x","variable.x","value.x","value.y"))


# Getting all models which have no RNA expression
models_no_expr <- setdiff(unique(new.subset$variable.x), colnames(rna.mat))

# Writing dataframe to file
write.table(models_no_expr, paste0(pptc.folder, "Manuscript/scripts/focal-cn/", Sys.Date(), "-models-no-RNA-expr.txt"),
            sep = "\t", col.names = F, row.names = F, quote = F)


# Applying conditions to filter out Hom_Deletions which are greated than 5
filtered_new_less_than_5_HD <- setDT(new.subset)[value.y < 1 & value.x == "Hom_Deletion"| (is.na(value.y) & value.x == "Hom_Deletion"),]
filtered_new_remaining <- new.subset[value.x=="Amplification" | value.x=="Hem_Deletion" | is.na(value.x),]

filtered_new <- rbind(filtered_new_less_than_5_HD,filtered_new_remaining)


filtered_new <- setDT(filtered_new)[(!is.na(gene_short_name.x) & !is.na(variable.x)),]

# Renaming columns in filtered_new 
setnames(filtered_new,c("gene_short_name.x","variable.x", "value.x","value.y"),c("gene_short_name","variable","alteration","FPKM"))

# casting dataframe 
final <-as.data.frame(dcast(unique(filtered_new), gene_short_name ~ variable,value.var = "alteration"))
final[is.na(final)] <- NA

rownames(final) <- final$gene_short_name  
final$gene_short_name <- NULL
large.cn<- as.matrix(final)
######ADD HEMIZYGOUS DELETIONS FROM GISTIC OUTPUT

###remove rows with chr due to duplicates
gistic.cn <- gistic.out[grep("CDKN2A$|CDKN2B$|SMARCB1$|TP53$",gistic.out$`Gene Symbol`,invert=F),]
rownames(gistic.cn) <- gistic.cn$`Gene Symbol`
gistic.cn <- gistic.cn[,4:ncol(gistic.cn)]
##add only hemi dels
gistic.cn[gistic.cn == 0] <- ""
gistic.cn[gistic.cn == 2] <- ""
gistic.cn[gistic.cn == -2] <- ""
gistic.cn[gistic.cn == 1] <- ""


####add hemi del for CDKN2A/B for leukemia models    
leuk <- subset(clin, Histology.Oncoprints == "leukemia") 

for(i in which(colnames(gistic.cn) %in% as.list(leuk$Model))){ 
  gistic.cn["CDKN2A",i][gistic.cn["CDKN2A",i] == -1] = "Hem_Deletion"
  gistic.cn["CDKN2B",i][gistic.cn["CDKN2B",i] == -1] = "Hem_Deletion"
}

###add het for SMARCB1 for brain only
brain.df <- subset(clin, Histology.Oncoprints == "brain") 
for(i in which(colnames(gistic.cn) %in% as.list(brain.df$Model))){ 
  gistic.cn["SMARCB1",i][gistic.cn["SMARCB1",i] == -1] = "Hem_Deletion"
}

###add het for TP53 for osteo only
ost.df <- subset(clin, Histology.Oncoprints == "osteosarcoma") 
for(i in which(colnames(gistic.cn) %in% as.list(ost.df$Model))){ 
  gistic.cn["TP53",i][gistic.cn["TP53",i] == -1] = "Hem_Deletion"
}


##remove other -1s 
gistic.cn[gistic.cn == -1] <- ""

###merge hemi matrix with current CN matrix
source(paste0(script.folder, "merge-CN-gistic-matrices-new.R"))

###fix ALL-50 - should be hom del for both CDKN2A/B
if("ALL-50" %in% colnames(gistic.cn)) {
  New.gistic.cn["CDKN2B","ALL-50"] <- "Hom_Deletion"
}
###fix PALKTY - should be het del for both CDKN2B
if("PALKTY" %in% colnames(gistic.cn)) {
  New.gistic.cn["CDKN2B","PALKTY"] <- "Hem_Deletion"
}
###fix OS-52-SJ - should be hom del for TP53
if("OS-52-SJ" %in% colnames(gistic.cn)) {
  New.gistic.cn["TP53","OS-52-SJ"] <- "Hom_Deletion"
}


###subset for only models in PPTC (theoretically should be done in master files, but just in case)
New.gistic.cn <- New.gistic.cn[,which(colnames(New.gistic.cn) %in% clin$Model)]

###remove empty rows
short.cn.mat <- New.gistic.cn[!apply(New.gistic.cn == "", 1, all),] 
##Write new CN matrix for oncoprints
write.table(as.data.frame(short.cn.mat), paste0(cnDir, Sys.Date(), "-short_cn_matrix_fpkm1.txt"), quote = F,col.names = T,row.names = T, sep = "\t")

####code back to numeric for pedcbio
pedc.cn <- New.gistic.cn
pedc.cn[pedc.cn == ""] <- 0
pedc.cn[pedc.cn == "Hom_Deletion"] <- -2
pedc.cn[pedc.cn == "Hem_Deletion"] <- -1
pedc.cn[pedc.cn == "Amplification"] <- 2
##Write new CN matrix for oncoprints
write.table(as.data.frame(pedc.cn), paste0(pptc.folder, "Pedcbio-upload/", Sys.Date(), "-focal-cn-fpkm1.txt"), quote = F,col.names = T,row.names = T, sep = "\t")
