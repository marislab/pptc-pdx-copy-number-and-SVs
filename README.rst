
.. |date| date::

***************************************
PPTC PDX Focal Copy Number and SV plots
***************************************

:authors: Jo Lynne Rokita, Komal S. Rathi
:contact: Jo Lynne Rokita (rokita@email.chop.edu) or Komal S. Rathi (rathik@email.chop.edu)
:organization: CHOP
:status: Completed
:date: |date|

.. meta::
   :keywords: pdx, SNP array, copy number, chromothripsis 2019
   :description: code to create focal copy number matrix, SV plots, and breakpoint density plots



Pipeline
========

.. code-block:: bash

         # How to run:
         # Download github repository in your home directory (~/)
         # Make sure to not clone the repository inside any other repository in your home directory
         git clone https://github.com/marislab/pptc-pdx-copy-number-and-SVs.git
         
         # Run script to fetch data files
         ./data-fetch-pptc-pdx-copy-number-SV.sh
         
         # Run scripts
         Rscript focal-CN.R
         Rscript sv_figure.R
