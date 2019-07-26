################################################
# Author: Komal S Rathi
# Function: Pedcbio formatting scripts
# Date: 01/23/2019
################################################

setwd('~/Projects/Maris-lab/PPTC')
library(biomaRt)
library(dplyr)
detach('package:plyr')

####### clinical - patient and samples data (n = 261)
clin <- read.delim('data/pedcbio/2019-07-25-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
expr <- clin[which(clin$RNA.Part.of.PPTC == "yes"),'Model']
clin <- clin[,c('PersonID','Model','Model',
        'Histology.pedcbio','Histology.Detailed',
        'PI','Sex','Phase','Age',
        'Site.of.Initial.Tumor','Site.of.Specimen','New.Histopathology',
        'Reported_Ethnicity','Inferred_Ethnicity','Molecular.Subtype')]
colnames(clin) <- c("PATIENT_ID","SAMPLE_ID","TUMOR_SAMPLE_BARCODE",
                    "CANCER_TYPE","CANCER_TYPE_DETAILED",
                    "PI","GENDER","PHASE","AGE",
                    "SITE_OF_INITIAL_TUMOR","SITE_OF_SPECIMEN","HISTOPATHOLOGY",
                    "REPORTED_ETHNICITY","INFERRED_ETHNICITY","MOLECULAR_SUBTYPE")
clin$AGE[clin$AGE %in% c("Unknown","Uknown","Unavailable")] <- ''

patient <- unique(clin[,c('PATIENT_ID','GENDER','AGE','REPORTED_ETHNICITY','INFERRED_ETHNICITY')])
patient <- patient %>% group_by(PATIENT_ID, GENDER, REPORTED_ETHNICITY, INFERRED_ETHNICITY) %>%
  summarise(AGE = min(AGE)) %>%
  as.data.frame()
samples <- clin[,c('PATIENT_ID','SAMPLE_ID','TUMOR_SAMPLE_BARCODE',
                   'CANCER_TYPE','CANCER_TYPE_DETAILED','PI','PHASE',
                   'SITE_OF_INITIAL_TUMOR','SITE_OF_SPECIMEN','HISTOPATHOLOGY','MOLECULAR_SUBTYPE')]

# add expression
samples$EXPRESSION <- 'NO'
samples$EXPRESSION[samples$PATIENT_ID %in% expr] <- "YES"

# add drug response
drug <- read.delim('data/pedcbio/2019-01-21-47-drugresponse-pptc.txt', stringsAsFactors = F, check.names = F)
colnames(drug)[1] <- 'TUMOR_SAMPLE_BARCODE'
colnames(drug) <- toupper(colnames(drug))
colnames(drug)[2:ncol(drug)] <- paste0('TX_',colnames(drug)[2:ncol(drug)])
samples <- merge(samples, drug, by = 'TUMOR_SAMPLE_BARCODE', all.x = TRUE)
samples[,13:ncol(samples)][is.na(samples[,13:ncol(samples)])] <- ''
samples <- samples[,c(3, 4, 5, 7, 8, 9, 10, 11, 1, 6, 2, 12, 13:ncol(samples))]

patient <- rbind(colnames(patient), colnames(patient), 
                 c("STRING","STRING","STRING","STRING","NUMBER"), 
                 c(2,2,2,1,2), 
                 colnames(patient), patient)
samples <- rbind(colnames(samples), colnames(samples), 
                 rep("STRING", ncol(samples)), 
                 c(6, 5, 5, 4, 4, 4, 3, 3, 2, 2, 2, 1, rep(0, length(13:ncol(samples)))), 
                 colnames(samples), samples)
samples[is.na(samples)] <- ''
patient[is.na(patient)] <- ''
patient[1:4,1] <- paste0('#',patient[1:4,1])
samples[1:4,1] <- paste0('#',samples[1:4,1])

# write out the data
write.table(samples, file = 'data/pedcbio/pptc/data_clinical_samples.txt', quote = F, sep = "\t", row.names = F, col.names = F)
write.table(patient, file = 'data/pedcbio/pptc/data_clinical_patients.txt', quote = F, sep = "\t", row.names = F, col.names = F)

####### mrna-sequencing (n = 244) 
mat <- get(load('data/pedcbio/2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda'))
mat$gene_short_name <- as.character(mat$gene_short_name)
mat$gene_short_name[mat$gene_short_name == "MIR219-1"] <- 'MIR219A1'
mat$gene_short_name[mat$gene_short_name == "MIR219-2"] <- 'MIR219A2'
genes <- read.delim('~/Projects/cavatica/cBio/genes.tsv')
genes <- genes[,c('HUGO_GENE_SYMBOL','ENTREZ_GENE_ID')]
colnames(genes) <- c('Hugo_Symbol','Entrez_Gene_Id')
mat <- merge(genes, mat, by.x = 'Hugo_Symbol', by.y = 'gene_short_name', all.y = T)
mat$Entrez_Gene_Id[is.na(mat$Entrez_Gene_Id)] <- ""
dim(mat[which(mat$Entrez_Gene_Id == ""),])
setdiff(colnames(mat), clin$PATIENT_ID)

write.table(mat, file = 'data/pedcbio/pptc/data_rna_seq_mrna.txt', quote = F, sep = "\t", row.names = F)
system('mkdir -p data/pedcbio/pptc/case_lists')
write.table(colnames(mat)[3:ncol(mat)], file = 'data/pedcbio/pptc/case_lists/cases_rna_seq_mrna.txt', quote = F, row.names = F, sep = "\t", col.names = F)

# z-score the matrix
zscore <- function(x){
  z <- (x - mean(x)) / sd(x)
  return(z)
}
log.trans <- log2(mat[,3:ncol(mat)]+1)
res.zscore <- t(apply(log.trans, MARGIN = 1, zscore))
res.zscore <- cbind(mat[,1:2], res.zscore)
write.table(res.zscore, file = 'data/pedcbio/pptc/data_rna_seq_mrna_median_Zscores.txt', quote = F, sep = "\t", row.names = F)

####### Fusions (n = 166)
dat <- read.delim('data/pedcbio/pptc/data_fusions.txt', stringsAsFactors = F)
setdiff(dat$Tumor_Sample_Barcode, clin$SAMPLE_ID)
# write.table(dat, file = 'data/pedcbio/pptc/data_fusions.txt', quote = F, sep = "\t", row.names = F)

####### segmentation (n = 252)
seg <- read.delim('data/pedcbio/2019-02-10-252pdx-final-pptc-SNPRANK-noXY.seg', stringsAsFactors = F)
colnames(seg) <- c('ID','chrom','loc.start','loc.end','num.mark','seg.mean')
setdiff(seg$ID, clin$PATIENT_ID)
length(unique(seg$ID))
seg <- seg[-which(seg$chrom %in% c("X","Y")),]
write.table(seg, file = 'data/pedcbio/pptc/data_cna_hg19.seg', quote = F, sep = "\t", row.names = F)

####### MAF  (n = 240)
# remove all non-mandatory columns
load('data/pedcbio/2019-02-14-allpdx-clean-maf-240.rda')
all.cols <- c("Hugo_Symbol","Entrez_Gene_Id","Center","NCBI_Build",
              "Chromosome","Start_Position","End_Position","Strand",
              "Variant_Classification","Variant_Type","Reference_Allele",
              "Tumor_Seq_Allele1","Tumor_Seq_Allele2","dbSNP_RS",
              "dbSNP_Val_Status","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode",
              "Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2",
              "Tumor_Validation_Allele1","Tumor_Validation_Allele2",
              "Match_Norm_Validation_Allele1","Match_Norm_Validation_Allele2",
              "Verification_Status","Validation_Status","Mutation_Status",
              "Sequencing_Phase","Sequence_Source","Validation_Method",
              "Score","BAM_File","Sequencer","HGVSp_Short",
              "t_alt_count","t_ref_count","n_alt_count","n_ref_count","SWISSPROT","Protein_position")
req.cols <- c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode","HGVSp_Short")
pptc.merge <- pptc.merge[,grep('dbNSFP|HGNC|ESP|CGC|X1000|COSMIC|UniProt|CCLE|ClinVar|RefSeq|EXAC|AF|_rna|GO|TCGAscape_|Tumorscape|transcript|Familial|Ensembl', colnames(pptc.merge), invert = T)]
setdiff(req.cols, colnames(pptc.merge))
colnames(pptc.merge)[which(colnames(pptc.merge) == "Protein_Change")] <- "HGVSp_Short"
setdiff(all.cols, colnames(pptc.merge))
colnames(pptc.merge)[which(colnames(pptc.merge) == "Start_position")] <- "Start_Position"
colnames(pptc.merge)[which(colnames(pptc.merge) == "End_position")] <- "End_Position"
colnames(pptc.merge)[which(colnames(pptc.merge) == "Center.x")] <- "Center"
colnames(pptc.merge)[which(colnames(pptc.merge) == "SwissProt_acc_Id")] <- "SWISSPROT"
colnames(pptc.merge)[which(colnames(pptc.merge) %in% c("TVarCov", "NVarCov"))] <- c("t_alt_count","n_alt_count")
pptc.merge$t_ref_count <- pptc.merge$TTotCov-pptc.merge$t_alt_count
pptc.merge$n_ref_count <- pptc.merge$NTotCov-pptc.merge$n_alt_count
pptc.merge <- pptc.merge[,which(colnames(pptc.merge) %in% c(all.cols, req.cols))]
setdiff(pptc.merge$Tumor_Sample_Barcode,clin$PATIENT_ID)
pptc.merge[is.na(pptc.merge)] <- ''

# remove quotation marks
x <- apply(pptc.merge, MARGIN = 2, FUN = function(x) length(grep('"', x)))
colnames <- names(x[x > 0])
if(length(colnames) == 0){
  print("Do nothing")
} else {
  pptc.merge$DrugBank <- gsub('"','',pptc.merge$DrugBank)
  pptc.merge$HGNC_Gene.family.description <- gsub('"','',pptc.merge$HGNC_Gene.family.description)
  pptc.merge$HGNC_Name.Synonyms <- gsub('"','',pptc.merge$HGNC_Name.Synonyms)
  pptc.merge$HGNC_Previous.Names <- gsub('"','',pptc.merge$HGNC_Previous.Names)
}

# convert variant classification
vc <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Silent","Splice_Site","Translation_Start_Site","Nonstop_Mutation","3'UTR","3'Flank","5'UTR","5'Flank","IGR","Intron","RNA","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame","Splice_Region and Unknown")
pptc.merge$Variant_Classification <- as.character(pptc.merge$Variant_Classification)
pptc.merge$Variant_Classification[pptc.merge$Variant_Classification == "Start_Codon_SNP"] <- "Translation_Start_Site"
pptc.merge$Variant_Classification[pptc.merge$Variant_Classification == "lincRNA"] <- "RNA"
pptc.merge$Variant_Classification[pptc.merge$Variant_Classification == "Stop_Codon_Del"] <- "Nonstop_Mutation"
pptc.merge$Variant_Classification[pptc.merge$Variant_Classification == "Start_Codon_Ins"] <- "De_novo_Start_InFrame"
pptc.merge$Variant_Classification[pptc.merge$Variant_Classification == "Stop_Codon_Ins"] <- "Nonsense_Mutation"
pptc.merge$Variant_Classification[pptc.merge$Variant_Classification == "Start_Codon_Del"] <- "In_Frame_Del"
setdiff(pptc.merge$Variant_Classification, vc)
pptc.merge[is.na(pptc.merge)] <- ''
write.table(pptc.merge, file = 'data/pedcbio/pptc/data_mutations.txt', quote = F, sep = "\t", row.names = F)
write.table(unique(pptc.merge$Tumor_Sample_Barcode), file = 'data/pedcbio/pptc/case_lists/cases_sequenced.txt', quote = F, row.names = F, sep = "\t", col.names = F)

####### copy number (n = 252)
cn <- read.delim('data/pedcbio/2019-07-10-focal-cn-fpkm1.txt', stringsAsFactors = F, check.names = F)
cn$Hugo_Symbol <- rownames(cn)
genes <- read.delim('~/Projects/cavatica/cBio/genes.tsv')
genes <- unique(genes[,c(2,1)])
cn <- merge(cn, genes, by.x = 'Hugo_Symbol','HUGO_GENE_SYMBOL', all.x = T)
colnames(cn)[c(1,254)] <- c('Hugo_Symbol', 'Entrez_Gene_Id')
setdiff(colnames(cn), clin$SAMPLE_ID)
cn <- cn[grep('chrX|chrY|^hsa-mir|^hsa-let', cn$Hugo_Symbol, invert = T),]
cn$Hugo_Symbol2 <- gsub("[|]|chr.*","",cn$Hugo_Symbol)
# cn <- unique(cn)
# cn$Entrez_Gene_Id[cn$Entrez_Gene_Id < 0] <- ''
# cn$Hugo_Symbol <- ifelse(grepl("orf", cn$Hugo_Symbol), toupper(cn$Hugo_Symbol), cn$Hugo_Symbol)
# genes <- read.delim('~/Projects/cavatica/cBio/genes.tsv')
# cn.annot <- cn[,1:2]
# setdiff(cn.annot$Hugo_Symbol, genes$HUGO_GENE_SYMBOL)

# biomart of hg38
tmp <- plyr::count(cn$Hugo_Symbol2)
tmp <- tmp[which(tmp$freq > 1),]
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
results = getBM(attributes=c("hgnc_symbol","chromosome_name"), 
                filters="hgnc_symbol",
                values=tmp$x, 
                mart=ensembl);
results <- results[grep('CHR',results$chromosome_name,invert = T),]
results$chromosome_name <- paste0('chr',results$chromosome_name)
results$hgnc_symbol <- paste0(results$hgnc_symbol,'|',results$chromosome_name)
remove.genes <- cn[which(cn$Hugo_Symbol2 %in% tmp$x),'Hugo_Symbol']
remove.genes <- remove.genes[!remove.genes %in% results$hgnc_symbol]
cn <- cn[-which(cn$Hugo_Symbol %in% remove.genes),]
cn$Hugo_Symbol <- cn$Hugo_Symbol2
cn$Hugo_Symbol2 <- NULL
cn <- unique(cn)
cn$Entrez_Gene_Id[cn$Entrez_Gene_Id < 0] <- ''
cn$Hugo_Symbol <- ifelse(grepl("orf", cn$Hugo_Symbol), toupper(cn$Hugo_Symbol), cn$Hugo_Symbol)
tmp <- plyr::count(cn$Hugo_Symbol)
cn <- cn[,c(1,254,2:253)]
write.table(cn, file = 'data/pedcbio/pptc/data_CNA.txt', quote = F, sep = "\t", row.names = F)
write.table(colnames(cn)[3:ncol(cn)], file = 'data/pedcbio/pptc/case_lists/cases_cna.txt', quote = F, sep = "\t", row.names = F, col.names = F)

