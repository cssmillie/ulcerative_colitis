# -----------------------------------------------------------------------------------
# Intra- and inter-cellular rewiring of the human colon during ulcerative colitis
# Smillie, C.S., Biton, M.B., Ordovas-Montanes J., et al., Cell, 2019.
# -----------------------------------------------------------------------------------
#
# This repository contains code for:
# - clustering single cells into cell subsets
# - detecting ambient rna contamination
# - calculating significant changes in cell composition with disease
# - calculating significant changes in gene expression with disease
#
# For additional code or questions please contact Chris Smillie (cssmillie@gmail.com)


# --------------
# 0) Basic setup
# --------------


# Source "analysis.r", which loads the necessary R packages and functions
source('analysis.r')

# Read a list of cell subsets, including the group that each one belongs to
# Groups include: Epithelial, Endothelial, Fibroblasts, Glia, Myeloid, B, and T cells
cell_subsets = read.table('cell_subsets.txt', sep='\t', header=F, stringsAsFactors=F)



# ------------------------------------
# 1) Cell clustering and visualization
# ------------------------------------
# To cluster the cells, we follow the basic analysis steps:
# 1. calculate variable genes
# 2. batch correction (using ComBat)
# 3. pca
# 4. tsne
# 5. clustering (using phenograph)
# These steps are performed by the "run_analysis" function in analysis.r

# Alternatively, you can load the results from the discovery cohort using "skip_clustering = TRUE"
skip_clustering = FALSE


if(skip_clustering == FALSE){
    
    
    # -----------------------------------------------------------
    # Cluster the cells from the discovery and validation cohorts
    # -----------------------------------------------------------
    

    # read sparse counts matrices (rows = genes, columns = cell barcodes) and metadata tables
    # ---------------------------------------------------------------------------------------
    
    # epithelial counts
    epi.counts = readMM('gene_sorted-Epi.matrix.mtx')
    rownames(epi.counts) = readLines('Epi.genes.tsv')
    colnames(epi.counts) = readLines('Epi.barcodes2.tsv')
    
    # stromal counts
    fib.counts = readMM('gene_sorted-Fib.matrix.mtx')
    rownames(fib.counts) = readLines('Fib.genes.tsv')
    colnames(fib.counts) = readLines('Fib.barcodes2.tsv')
    
    # immune counts
    imm.counts = readMM('gene_sorted-Imm.matrix.mtx')
    rownames(imm.counts) = readLines('Imm.genes.tsv')
    colnames(imm.counts) = readLines('Imm.barcodes2.tsv')
    
    # Load metadata for discovery and validation cohorts
    meta = read.table('all.meta2.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
    
    
    # run sc-rna-seq analysis pipeline
    # --------------------------------
    # the arguments and parameters are as follows (see also STAR Methods):
    # - name = project name
    # - counts = subsampled expression matrix
    # - ming = remove cells with < x genes (ming = 250)
    # - var_regex = remove mitochondrial and ribosomal genes (see below)
    # - do.batch = perform batch correction with ComBat
    # - num_pcs = number of PCs to use (num_pcs = 20 for initial clustering)
    # - max_iter = number of tsne iterations
    # - k = parameter for phenograph clustering (epithelial k = 750, stromal k = 250, immune k = 250)
    # the output is a seurat object
    # because this pipeline takes a long time to run, we subsample 25,000 cells in this example
    
    # get sample names for batch correction
    epi.batch = setNames(meta[colnames(epi.counts), 'Sample'], colnames(epi.counts))
    fib.batch = setNames(meta[colnames(fib.counts), 'Sample'], colnames(fib.counts))
    imm.batch = setNames(meta[colnames(imm.counts), 'Sample'], colnames(imm.counts))
    
    # run analysis pipeline on data subset
    var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP' # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
    epi.seur = run_analysis(name='epi', counts=epi.counts[,sample(1:ncol(epi.counts), 25000)], ming=250, var_regex=var_regex, do.batch=TRUE, batch.use=epi.batch, num_pcs=20, max_iter=1000, k=750)
    fib.seur = run_analysis(name='fib', counts=fib.counts[,sample(1:ncol(fib.counts), 25000)], ming=250, var_regex=var_regex, do.batch=TRUE, batch.use=fib.batch, num_pcs=20, max_iter=1000, k=250)
    imm.seur = run_analysis(name='imm', counts=imm.counts[,sample(1:ncol(imm.counts), 25000)], ming=250, var_regex=var_regex, do.batch=TRUE, batch.use=imm.batch, num_pcs=20, max_iter=1000, k=250)
    
    # set cell clusters within seurat object
    epi.seur = SetIdent(epi.seur, ident.use=epi.seur@data.info$phenograph.k750)
    fib.seur = SetIdent(fib.seur, ident.use=fib.seur@data.info$phenograph.k250)
    imm.seur = SetIdent(imm.seur, ident.use=imm.seur@data.info$phenograph.k250)
    
} else {
    
    
    # Alternatively, load the discovery cohort data from our paper
    # ------------------------------------------------------------
    # The analysis pipeline has steps that are stochastic and can also take a long time to run
    # To reproduce the results in our paper, you can load our original Seurat object
    # Note that this object contains data from the "discovery" cohort and not the "validation" cohort
    
    # load seurat objects
    epi.seur = readRDS('train.Epi.seur.rds')
    fib.seur = readRDS('train.Fib.seur.rds')
    imm.seur = readRDS('train.Imm.seur.rds')
    
    # set counts matrices
    epi.counts = epi.seur@assays[['RNA']]@counts
    fib.counts = fib.seur@assays[['RNA']]@counts
    imm.counts = imm.seur@assays[['RNA']]@counts

}



# ------------------------
# 2) Ambient RNA detection
# ------------------------
# Ambient RNA contamination can be a problem for droplet-based single cell data
# We developed a method to flag putative contaminants within each cell subset (see STAR Methods)
# This method first learns a contamination model across major cell "groups" (i.e. Epithelial, Endothelial, Fibroblast, Glia, Myeloid, B, and T cells)
# It then applies this model to each cell subset ("idents")
# This method is implemented in the "detect_contamination" function, which has the following arguments:
# - tpm = TPM (or TP10K) expression matrix (rows = genes, columns = cell barcodes)
# - groups = coarse-level cell groups/ cell lineages (e.g. T cells) for each cell barcode
# - idents = fine-level cell subsets (e.g. CD8+ IELs, Tregs, etc.) for each cell barcode
# - samples = sample that each cell barcode was collected from
# The output of this function is an object containing information about the model fit, including the residuals (low residuals = putative contaminants)
# We later use these residuals to remove putative contaminants from our list of differentially expressed genes

# calculate tpm for all cells
data = sparse_cbind(list(epi.counts, fib.counts, imm.counts))
numi = colSums(data)
data = scaleMargins(data, cols=1e4/numi)

# get cell types and cell groups
idents = c(as.character(epi.seur@active.ident), as.character(fib.seur@active.ident), as.character(imm.seur@active.ident))
groups = data.frame(cell_subsets, row.names=1)[idents,1]
samples = c(epi.seur@meta.data$Sample, fib.seur@meta.data$Sample, imm.seur@meta.data$Sample)

# detect ambient rna contamination (output = object describing fit)
if(FALSE){
    contam = detect_contamination(tpm=data, groups=groups, idents=idents, samples=samples, do.plot=FALSE)
} else {
    contam = readRDS('train.contam.rds')
}

# get contamination residuals (low residuals = putative contaminants)
contam.res = as.data.table(as.data.frame(sapply(contam$res.idents, function(a) a$residuals)) %>% rownames_to_column('gene') %>% gather(ident, res, -gene))



# ---------------------------------------------
# 3) Changes in cell proportions during disease
# ---------------------------------------------


# Use dirichlet-multinomial regression to find significant changes in cell frequencies during disease
# ---------------------------------------------------------------------------------------------------

# Count each cell subset in every sample
epi.freq = as.matrix(as.data.frame.matrix(table(epi.seur@meta.data$Sample, epi.seur@meta.data$Cluster)))
fib.freq = as.matrix(as.data.frame.matrix(table(fib.seur@meta.data$Sample, fib.seur@meta.data$Cluster)))
imm.freq = as.matrix(as.data.frame.matrix(table(imm.seur@meta.data$Sample, imm.seur@meta.data$Cluster)))

# Combine counts into a single matrix
all.freq = sparse_cbind(list(epi.freq, fib.freq, imm.freq))

# For the validation cohort, we need to combine the replicate samples (because they are not independent)
# To construct a list of replicates, we remove "1", "2", "a", and "b" from the sample IDs
reps = gsub('[12]*[ab]*$', '', rownames(all.freq))
temp = as.matrix(data.frame(aggregate(as.matrix(all.freq), list(reps), sum), row.names=1))
colnames(temp) = colnames(all.freq)
all.freq = temp[,colSums(temp) > 0]

# Split matrix into "epithelial" and "lamina propria" cell subsets and samples
ep.ident = levels(epi.seur@active.ident)
lp.ident = c(levels(fib.seur@active.ident), levels(imm.seur@active.ident))
ep.freq = all.freq[grep('Epi', rownames(all.freq)), ep.ident]
lp.freq = all.freq[grep('LP', rownames(all.freq)), lp.ident]

# For the dirichlet-multinomial regression, we need to know the disease state for each sample
# We can get this from the metadata table as follows:
sample2health = data.frame(unique(data.frame(sample=gsub('[12]*[ab]*$', '', meta[,'Sample']), health=meta[,'Health'])), row.names=1)
ep.cov = data.frame(condition=factor(sample2health[rownames(ep.freq),1], levels=c('Healthy', 'Non-inflamed', 'Inflamed')), row.names=rownames(ep.freq))
lp.cov = data.frame(condition=factor(sample2health[rownames(lp.freq),1], levels=c('Healthy', 'Non-inflamed', 'Inflamed')), row.names=rownames(lp.freq))

# Calculate significant changes using dirichlet multinomial regression
# This returns a matrix of p-values for each cell type / disease state
ep.pvals = dirichlet_regression(counts=ep.freq, covariates=ep.cov, formula=counts ~ condition)$pvals
colnames(ep.pvals) = colnames(ep.freq)
lp.pvals = dirichlet_regression(counts=lp.freq, covariates=lp.cov, formula=counts ~ condition)$pvals
colnames(lp.pvals) = colnames(lp.freq)

# Plot epithelial cell proportions
ep.pct = 100*ep.freq/rowSums(ep.freq)
p1 = matrix_barplot(ep.pct, group_by=ep.cov$condition, pvals=ep.pvals, colors=set.colors)
save_plot(p1, file='train.Fig2A.epi_freqs.pdf', nrow=1, ncol=2.5)

# Plot lamina propria cell proportions
lp.pct = 100*lp.freq/rowSums(lp.freq)
p2 = matrix_barplot(lp.pct, group_by=lp.cov$condition, pvals=lp.pvals, colors=set.colors)
save_plot(p2, file='train.Fig2A.lp_freqs.pdf', nrow=1, ncol=2.5)



# ------------------------------------------------
# 4) Differentially expressed genes during disease
# ------------------------------------------------
# Use MAST to find differentially expressed genes during disease while controlling for technical covariates
# This code uses the function "p.find_markers" in "markers.r" with the following arguments:
# seur = seurat object to use
# ident.1 = name of cell subset to test
# ident.use = vector of cell subsets to use
# min_alpha = only test genes that are expressed by at least [min_alpha] fraction of cells
# min_fc = only test genes with a mean fold change of [min_fc] between health and disease
# max_cells = subsample [max_cells] from each cell subset
# covariates = matrix of covariates
# formula = regression formula for DE test
# After finding the differentially expressed genes, we remove putative contaminants, then make violin plots

# Build covariates matrices
epi.cov = data.frame(nGene = scale(epi.seur@meta.data$nGene), Health = factor(epi.seur@meta.data$Health, levels=c('Healthy', 'Uninflamed', 'Inflamed')), Sample = epi.seur@meta.data$Sample)
fib.cov = data.frame(nGene = scale(fib.seur@meta.data$nGene), Health = factor(fib.seur@meta.data$Health, levels=c('Healthy', 'Uninflamed', 'Inflamed')), Sample = fib.seur@meta.data$Sample)
imm.cov = data.frame(nGene = scale(imm.seur@meta.data$nGene), Health = factor(imm.seur@meta.data$Health, levels=c('Healthy', 'Uninflamed', 'Inflamed')), Sample = imm.seur@meta.data$Sample)

# Construct formula for DE test. This contains the following terms:
# - nGene = cell complexity (number of genes per cell)
# - ident = cell subset
# - ident:Health = interaction term for the effect of the disease in each cell subset
formula = '~ nGene + ident + ident:Health'

# Calculate marker genes with MAST
# --------------------------------
# This uses the function "p.find_markers" in "markers.r", which helps to run MAST on the Seurat object (see arguments above)
# We run this separately for each cell compartment (epithelial/stromal/immune), then combine the results with "rbind"
# The resulting table contains information about the differential expression of each gene, including the following columns:
# - gene = gene name
# - contrast = model coefficient that was tested
# - coefD, pvalD, padjD = coefficient, p-value, and adjusted p-value for discrete part of the model
# - coefC, pvalC, padjC = coefficient, p-value, and adjusted p-value for continuous part of the model
# - mastfc, pvalH, padjH = fold change, p-value, and adjusted p-value for combined model terms (discrete + continuous parts)
# - n = number of cells expressing the gene in target group
# - alpha = fraction of cells expressing the gene in target group
# - mu = mean non-zero expression level of the gene in target group
# - mean = mean expression level of the gene in target group
# - total = total expression level of the gene in target group
# - log2fc = log2 fold change of the gene in target group
# - ident = cell subset
# - ref, ref* = reference group that the "target" group is compared against

# epithelial cells
mast.epi = sapply(levels(epi.seur@active.ident), function(ident){
    ident.use = factor(ifelse(epi.seur@active.ident == ident, ident, 'Other'), levels=c('Other', ident))
    p.find_markers(seur=epi.seur, ident.1=ident, ident.use=ident.use, min_alpha=.001, min_fc=1.1, max_cells=2500, covariates=epi.cov, formula=formula, lrt_regex=ident)
}, simplify=F)
mast.epi = do.call(rbind, mast.epi)

# stromal cells
mast.fib = sapply(levels(fib.seur@active.ident), function(ident){
    ident.use = factor(ifelse(fib.seur@active.ident == ident, ident, 'Other'), levels=c('Other', ident))
    p.find_markers(seur=fib.seur, ident.1=ident, ident.use=ident.use, min_alpha=.001, min_fc=1.1, max_cells=2500, covariates=fib.cov, formula=formula, lrt_regex=ident)    
}, simplify=F)
mast.fib = do.call(rbind, mast.fib)

# immune cells
mast.imm = sapply(levels(imm.seur@active.ident), function(ident){
    ident.use = factor(ifelse(imm.seur@active.ident == ident, ident, 'Other'), levels=c('Other', ident))
    p.find_markers(seur=imm.seur, ident.1=ident, ident.use=ident.use, min_alpha=.001, min_fc=1.1, max_cells=2500, covariates=imm.cov, formula=formula, lrt_regex=ident)    
}, simplify=F)
mast.imm = do.call(rbind, mast.imm)

# concatenate marker lists into a single data.table
mast.all = rbind(mast.epi, mast.fib, mast.imm)


# Remove ambient RNA contamination from marker lists
# --------------------------------------------------
# Use the contamination model (see "Ambient RNA detection" above) to filter putative contaminants from the DE gene lists
# To do this, we first add the contamination residuals to the table of DE genes, then we filter based on these residuals

# Merge DE genes and contamination residuals along "ident" and "gene" columns
mast.all = merge(mast.all, contam.res, by=c('ident', 'gene'))

# Split markers into putative contaminants ("mast.con") and clean ("mast.all") lists
# We conservatively flag genes with residuals < 5 (representing a 32-fold increase in expression over the estimated contamination rate)
mast.con = mast.all[res <= 5]
mast.all = mast.all[res >  5]


# Make volcano plots of DE genes (for "contaminating" and "clean" marker tables)
# ------------------------------------------------------------------------------
# In this example, we use DE genes that were associated with enterocytes during inflammation

# Select DE genes that are putative contaminants ("m1") and non-contaminants ("m2")
m1 = mast.con[contrast == 'identEnterocytes:HealthInflamed']
m2 = mast.all[contrast == 'identEnterocytes:HealthInflamed']

# Make volcano plots of putative contaminants ("p1") and non-contaminants ("p2")
# This uses the "simple_scatter" function in "plot.r", which has the following arguments:
# - x = x coordinate for each point
# - y = y coordinate for each point
# - lab = label for each point
# - lab.n = number of points to automatically label
# - lab.size = label size
# - xlab = x axis label
# - ylab = y axis label
# - ggtitle = add title to plot using ggplot
p1 = simple_scatter(x=m1$coefD, y=-log10(m1$padjD), lab=m1$gene, lab.n=25, lab.size=3, xlab='Fold change (inflamed vs. healthy)', ylab='-log10(adjusted P-value)') + ggtitle('Contamination')
p2 = simple_scatter(x=m2$coefD, y=-log10(m2$padjD), lab=m2$gene, lab.n=25, lab.size=3, xlab='Fold change (inflamed vs. healthy)', ylab='-log10(adjusted P-value)') + ggtitle('Non-contamination')

# Save volcano plots to file
ps = plot_grid(p1, p2, nrow=1)
save_plot(ps, file='train.enterocytes_volcano.pdf', nrow=1.25, ncol=2.75)

# Save marker table to file
# -------------------------

# Split "cell type" and "disease state" markers
mast.ct = mast.all[grep('Health', contrast, invert=T)] # cell type markers
mast.du = mast.all[grep('HealthUninflamed', contrast, invert=T)] # disease uninflamed
mast.di = mast.all[grep('HealthInflamed', contrast, invert=T)] # disease inflamed
