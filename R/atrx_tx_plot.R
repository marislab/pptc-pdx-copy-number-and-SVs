set.seed(12345)
library(ggplot2)
library(ggbeeswarm)
library(plotly)
library(ggforce)

data <- read.delim("2018-11-13_FPKM_3_transcripts_per_exon.txt")
options(scipen = 999)
# Number of Models 
length(unique(data$Model))


# total number of Models = 228
n_pages_needed = 76
pdf("2019-02-08-ALL-PDX-FPKM-exon-number-ATRX.pdf", width=20, height=9)

for (i in seq_len(n_pages_needed)) {
  print(ggplot(data,aes(data$exon_n, data$FPKM, group = transcript, color = transcript)) +
          geom_bar(aes(fill = transcript), stat = "identity", width = 0.5) +
          geom_hline(aes(yintercept = 1),linetype = "dashed") +
          scale_x_continuous(breaks = pretty(data$exon_n, n = 36)) +
          theme_bw() +
          labs(title = "FPKM for exons in all PDX models", x = "exon number", y = "FPKM") + 
          theme(legend.text=element_text(size=15), 
                axis.title.x = element_text(size=15), 
                axis.title.y = element_text(size=15),
                plot.title = element_text(hjust = 0.5,size = 20), 
                legend.title = element_text(size=15)) +
          facet_grid_paginate(data$transcript ~ data$Model, scales = "fixed", ncol = 3 , nrow = 3, page = i) + 
          theme(strip.text.x = element_text(size = 12)) +
          theme(strip.text.y = element_blank()) )
  
}
dev.off()