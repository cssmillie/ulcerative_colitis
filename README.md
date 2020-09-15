Intra- and inter-cellular rewiring of the human colon during ulcerative colitis
Smillie, C.S., Biton, M.B., Ordovas-Montanes J., et al., Cell, 2019.
-------------------------------------------------------------------------------

This repository contains code for:
- Clustering single cells into cell subsets
- Detecting ambient RNA contamination
- Calculating significant changes in cell composition with disease
- Calculating significant changes in gene expression with disease

For additional code or questions please contact Chris Smillie (cssmillie@gmail.com)

## Getting started


1) Clone github repository and "cd" into the "ulcerative_colitis" directory:
```
git clone https://github.com/cssmillie/ulcerative_colitis.git
cd ulcerative_colitis
```


2) Install necessary R packages
```
Rscript install.r
```
If you have trouble installing any of these packages through CRAN, you will need to do so manually.



3) Download the expression data from the Single Cell Portal (download into the "ulcerative colitis" GitHub directory)
https://portals.broadinstitute.org/single_cell/study/SCP259

The files you need are:

**Metadata**
- all.meta2.txt

**Epithelial data**

- Epi.genes.tsv
- Epi.barcodes2.tsv
- gene_sorted-Epi.matrix.mtx

**Stromal data**

- Fib.genes.tsv
- Fib.barcodes2.tsv
- gene_sorted-Fib.matrix.mtx

**Immune data**

- Imm.genes.tsv
- Imm.barcodes2.tsv
- gene_sorted-Imm.matrix.mtx


4) Download Seurat objects for discovery cohort (download into the "ulcerative_colitis" GitHub directory):
https://www.dropbox.com/sh/dn4gwdww8pmfebf/AACXYu8rda5LoLwuCZ8aZXfma?dl=0

The files you need are:

- train.Epi.seur.rds
- train.Fib.seur.rds
- train.Imm.seur.rds


5) After everything has been downloaded into the same directory ("ulcerative_colitis"), you can follow the code example in the "run.r" script. This takes you through the steps of the analysis pipeline, including cell clustering, ambient RNA detection, estimating significant changes in cell frequencies with disease, and performing differential expression tests.
