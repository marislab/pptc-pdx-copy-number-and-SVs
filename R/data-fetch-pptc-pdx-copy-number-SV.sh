#!/bin/bash
# Data files for pptc-pdx-copy-number and SVs
# INPUT files
# Expression matrix = 2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda
# gistic output = 2018-08-09-gistic-results-256pdx-noXY-snpfast2-nomirna/all_thresholded.by_genes.txt
# clinical file = 2019-02-09-pdx-clinical-final-for-paper.txt
# Approved hugo symbols file = 2019-02-14-Hugo-Symbols-approved.txt
# gtf = Homo_sapiens.GRCh37.87.gtf
# seg file = 2019-02-10-252pdx-final-pptc-SNPRANK-noXY.seg
# sourcing a script = merge-CN-gistic-matrices-new.R
# hist colors = 2019-02-09-all-hist-colors.txt
# mutations per model = 2019-02-13-mutations-per-model.txt


cd ~
mkdir -p copy-number-and-SVs/data/
cd ~/copy-number-and-SVs/data/

# 1. Download Expression matrix2
wget --output-document='2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda' https://ndownloader.figshare.com/files/14452985

# 2. Clinical file
wget --output-document='pptc-pdx-clinical-web.txt' https://ndownloader.figshare.com/files/14508536

# 3. Approved Hugo Symbols file
wget --output-document='2019-02-14-Hugo-Symbols-approved.txt' https://ndownloader.figshare.com/files/14460317

# 4. Seg file 
wget --output-document='2019-02-10-252pdx-final-pptc-SNPRANK-noXY.seg' https://ndownloader.figshare.com/files/14414165

# 5. Hist colors file
wget --output-document='2019-02-09-all-hist-colors.txt' https://ndownloader.figshare.com/files/14508539


