###### create cohesive oncoprints to regenerate Figures 2 and 3########
#
#     Authors: Jo Lynne Rokita, Gonzalo Lopez     
#     Updated 2019-12-27
################################################################

# working directory (created with git clone)
mainDir <- "~/pptc-pdx-copy-number-and-SVs/"
dataDir <- paste0(mainDir,"data/")
# set path to your git cloned repo
script.folder <- paste0(mainDir, "R/") 
# create new directories in mainDir
subDir <- paste0(mainDir,"results/")
ifelse(!dir.exists(file.path(subDir)), dir.create(file.path(subDir)), "Directory exists!")


##set wd
setwd(mainDir)

####Dependencies
source(paste0(script.folder, "install-packages.R"))

# obtain chromosome lengths and the cumulative position of every chr in the genome; no X, Y
chrlengths <- seqlengths(Hsapiens)[paste("chr",c(1:22) ,sep="")]

chr_genome_pos<-list()
for(i in 1:length(chrlengths)) chr_genome_pos[[i]] <- sum(as.numeric(c(0,chrlengths)[1:(i) ]))
chr_genome_pos<-unlist(chr_genome_pos)
names(chr_genome_pos) <-names(chrlengths)

# read gene annotation matrix from matlab (the same used to run GISTIC)
genemat <- read.mat(paste0(dataDir, "hg19.mat"))
genepos <- unlist(genemat$rg$start)
genechr <- unlist(genemat$rg$chr)
names(genepos) <-names(genechr) <-unlist(genemat$rg$symb)

clin <- read.delim(paste0(mainDir, "data/pptc-pdx-clinical-web.txt"), as.is = T, header = T)
clin <- clin[,c("Model", "Histology.Detailed")]
clin$Histology <- ifelse(clin$Histology.Detailed == "Ph+-ALL" | clin$Histology.Detailed == "Ph-likeALL", "Ph+/-like-ALL", clin$Histology.Detailed)
hists <- as.list(unique(clin$Histology))

# load gene level copy number matrix from GISTIC - located in data folder
df <- read.delim(paste0(mainDir, "data/all_data_by_genes.txt"),as.is=TRUE,check.names=FALSE)
names(df)

#convert to matrix
genecn.df <- df[,c(4:ncol(df))]
rownames(genecn.df)<-df[,1]
genecn_annot <- genecn.df[,1:3]

###use only histologies with N >=10
hist.n <- as.data.frame(table(clin$Histology))
hist.10 <- subset(hist.n, Freq >=10)

##leukemia
hist.all <- subset(hist.10, Var1 == "T-ALL" | Var1 == "MLL-ALL" | Var1 == "Ph+/-like-ALL" | Var1 == "BCP-ALL")

##solid tumors
hist.10.solid <- subset(hist.10, Var1 != "T-ALL" & Var1 != "MLL-ALL" & Var1 != "Ph+/-like-ALL" &
                            Var1 != "BCP-ALL")

###df used
#df <- hist.all
df <- hist.10.solid

###plot
#pdf(paste(mainDir, "/results/", Sys.Date(), "-leukemkia-CN-plots.pdf", sep = ""), width = 6, height = 8)
pdf(paste(mainDir, "/results/", Sys.Date(), "-solid-CN-plots.pdf", sep = ""), width = 6, height = 8)
par(mfrow=c(5,1))
for (each in df$Var1){
  sub.clin <- subset(clin, Histology == each) 
  genecn.df2 <- genecn.df[,which(colnames(genecn.df) %in% sub.clin$Model)]
  genecn <- as.matrix(genecn.df2)
  dim(genecn)
  genecn <- genecn[grep("chr",rownames(genecn),invert=T),]

  # create vectors for every gene's chromosome and position relative to the genome
  chrvector <- genechr[rownames(genecn)]
  gencoord <- chr_genome_pos[chrvector] + genepos[rownames(genecn)]
  names(gencoord)<-rownames(genecn)

  # create a vector of colors alternating by chromosome for both gains and losses
  gaincol<-"red"
  losscol<-"blue"
  listcolgain <- listcolloss <- list()
  for(chr in unique(chrvector)){
  listcolgain[[chr]] <- rep(gaincol,length(which(chrvector == chr)))
  listcolloss[[chr]] <- rep(losscol,length(which(chrvector == chr)))
    if(gaincol == "red"){
  	gaincol<-"salmon"
  	losscol<-"lightblue"
  	}else{
	gaincol<-"red"
	losscol<-"blue"
	}
  }
  
  colgain <- unlist(listcolgain)
  colloss <- unlist(listcolloss)
  names(colgain)<-names(colloss) <-names(chrvector)

  # create vector with gain and loss frequencies for each gene
  gains <- apply(genecn,1,function(x) length(which(x > 0.1))/ncol(genecn) )
  losses <- apply(genecn,1,function(x) length(which(x < -0.1))/ncol(genecn) )

  # plot data
  plot(x=NULL,y=NULL,xlim=c(1,max(gencoord)),ylim=c(min(-losses),max(gains)),xaxt='n',bty='n', ann = F, yaxt='n')
  points(gencoord,gains[names(gencoord)],type='h',col=colgain[names(gencoord)])
  points(gencoord,-losses[names(gencoord)],type='h',col=colloss[names(gencoord)])
  grid(nx=NA,ny=NULL)
  title(main = paste(each, " (n = ", ncol(genecn), ")", sep = ""), xlab = "", ylab = "Log R Ratio")
}
dev.off()
