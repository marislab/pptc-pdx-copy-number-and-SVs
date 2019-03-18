########################################################
# Author: Komal S Rathi
# Function: Get rawdata to run copynumber and SV scripts
# Date: 03/18/2019
########################################################

system("mkdir -p data/figures")
system("wget -- output-document='data/2019-02-09-pdx-clinical-final-for-paper.txt' \
       https://ndownloader.figshare.com/files/14372783")
system("wget --output-document='data/figures/2019-02-09-all-hist-colors.txt' \
       https://ndownloader.figshare.com/files/14372849")
system("wget --output-document='data/figures/2019-02-10-252pdx-final-pptc-SNPRANK-noXY.seg' \
       https://ndownloader.figshare.com/files/14414165")
system("wget --output-document='ata/figures/2019-02-13-mutations-per-model.txt' \
       https://ndownloader.figshare.com/files/14612012")