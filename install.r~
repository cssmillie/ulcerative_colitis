# Load required packages and functions
# ------------------------------------
library(Seurat)
library(methods)

library(BiocGenerics, pos=length(search()))
library(S4Vectors, pos=length(search()))
library(DelayedArray, pos=length(search()))
library(MAST)

library(cccd)
library(cowplot)
library(data.table)
library(DirichletReg)
library(doParallel)
library(dplyr)
library(expm)
library(Hmisc)
library(ggbeeswarm)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(gtools)
library(MASS)
library(Matrix)
library(Matrix.utils)
library(rsvd)
library(Rtsne)
library(RColorBrewer)
library(sva)
library(tibble)
library(tidyverse)
library(tidyr)
library(wordspace)


source('contamination.r')
source('colors.r')
source('downsample.r')
source('markers.r')
source('mm_utils.r')
source('mtx.r')
source('parallel.r')
source('plot.r')
source('scores.r')
source('tpm.r')


mean_cv_loess = function(data, num_genes=1500, num_bins=20){
    
    # Get set of highly variable genes for a single sample
    # ----------------------------------------------------
    # data = log2(TP10K+1) expression matrix (rows = genes, columns = cell barcodes)
    # num_genes = number of variable genes to select
    # num_bins = number of expression bins to sample from
    
    # calculate mean and cv
    u = apply(data, 1, mean)
    v = apply(data, 1, var)
    i = u > 0 & v > 0
    u = u[i]
    v = v[i]
    cv = sqrt(v)/u
    
    # fit loess curve
    l = loess(log(cv) ~ log(u), family='symmetric')
    d = log(cv) - l$fitted
    
    # select variable genes from equal frequency bins
    k = as.integer(num_genes/num_bins)
    var_genes = as.character(unlist(tapply(d, cut2(u, g=num_bins), function(a) names(sort(a, decreasing=T)[1:k]))))
    
    var_genes

}


get_var_genes = function(data, samples, num_genes=1500, num_bins=20, min_cells=50){
    
    # Get consensus set of highly variable genes across multiple samples
    # ------------------------------------------------------------------
    # data = log2(TP10K+1) expression matrix (rows = genes, columns = cell barcodes)
    # num_genes = number of variable genes to select
    # num_bins = number of expression bins to sample from
    # remove genes that are expressed in < min_cells cells
    
    # Merge samples containing fewer than min_cells cells
    samples = as.factor(samples)
    levels(samples)[table(samples) <= min_cells] = 'Merge'

    # Calculate variable genes within each sample
    var_genes = sapply(levels(samples), function(i){
        
        # Subsample data
	data = data[,samples == i]
	data = data[rowSums(data > 0) >= min_cells,]
	mean_cv_loess(data, num_genes=num_genes, num_bins=num_bins)
    })
    
    # Select consensus gene set
    a = sort(table(as.character(unlist(var_genes))), decreasing=T)
    num_genes = min(length(a), num_genes)
    k = a[num_genes]
    u = names(a)[a > k]
    v = sample(names(a)[a == k], num_genes - length(u))
    var_genes = c(u,v)
    
    var_genes
}
    

run_phenograph = function(data, k=250){

    # Run phenograph (python code)
    # ----------------------------
    # data = input data (PCs)
    # k = number of nearest neighbors for knn graph
    
    # Write data to tempfile
    out = tempfile(pattern='phenograph.', tmpdir='.', fileext='.txt')
    write.table(data, file=out, sep='\t', quote=F)
    
    # Run phenograph
    system(paste0('python ~/code/single_cell/run_phenograph.py --data ', out, ' -k ', k, ' --metric cosine --out ', out))
    
    # Cleanup tempfile
    clusters = readLines(out)
    system(paste0('rm ', out))
    
    # Read and return clusters
    clusters = data.frame(x=as.character(clusters), row.names=rownames(data), stringsAsFactors=F)
    colnames(clusters) = c(paste0('phenograph.k', k))
    clusters
}


make_seurat = function(name, counts, minc=10, ming=250, maxg=1e6){
    
    # Filter expression matrix and make Seurat object
    # -----------------------------------------------
    # name = project name
    # counts = sparse counts matrix (rows = genes, columns = cell barcodes)
    # remove genes that appear in < minc cells
    # remove cells that contain < ming or > maxg genes
    
    # Filter genes
    j1 = colSums(counts > 0) >= ming
    j2 = colSums(counts > 0) <= maxg
    counts = counts[,(j1 & j2)]
    
    # Filter cells
    i = rowSums(counts > 0) >= minc
    counts = counts[i,]
    
    # Make Seurat object
    seur = CreateSeuratObject(counts=counts, project=name, min.cells=0, min.features=0, names.field=1, names.delim='\\.')
    
    # Calculate log2(TP10K+1)
    seur@assays[['RNA']]@data = calc_tpm(counts=seur@assays[['RNA']]@counts)
    seur@assays[['RNA']]@data@x = log2(seur@assays[['RNA']]@data@x + 1)
    
    return(seur)
}


run_analysis = function(name, counts, minc=10, ming=250, maxg=1e6, var_regex=NULL, num_genes=1500, do.batch=FALSE, batch.use=NULL, num_pcs=30, max_iter=1000, k=250){
    
    # Run basic single cell analysis on sparse counts matrix
    # ------------------------------------------------------
    # name = project name
    # counts = sparse counts matrix (rows = genes, columns = cell barcodes)
    # remove genes that appear in < minc cells
    # remove cells that contain < ming or > maxg genes
    # remove variable genes that match the pattern specified by var_regex
    # num_genes = number of variable genes to use
    # do.batch = perform batch correction with ComBat? (True or False)
    # batch.use = vector of batch labels for each cell
    # num_pcs = number of pcs to use
    # max_iter = number of tsne iterations
    # k = number of neighbors for knn graph clustering
    
    # Make Seurat object
    print('Making Seurat object')
    seur = make_seurat(name=name, counts=counts, minc=minc, ming=ming, maxg=maxg)
    
    # Select variable genes
    print('Selecting variable genes')
    var_genes = get_var_genes(seur@assays[['RNA']]@counts, samples=seur@active.ident, num_genes=num_genes)
    
    if(!is.null(var_regex)){
        var_genes = grep(var_regex, var_genes, invert=T, value=T)
    }
    seur@assays[['RNA']]@var.features = intersect(var_genes, rownames(seur@assays[['RNA']]@data))
    
    # Batch correction with ComBat
    if(do.batch == TRUE){

        print('Batch correcting with ComBat')
        
        # Get batch vector
        if(is.null(batch.use)){
	    batch.use = seur@active.ident
	}
	batch.use = batch.use[names(seur@active.ident)]
	
	# Batch correct with ComBat
	bc.data = ComBat(dat=as.matrix(seur@assays[['RNA']]@data), batch.use, par.prior=TRUE, prior.plots=FALSE)
	fwrite(as.data.table(bc.data), file=paste0(name, '.bc.data.txt'), sep='\t')	
	pc.data = t(scale(t(bc.data), center=F))

    } else {
        pc.data = seur@assays[['RNA']]@data
    }
    
    # Fast PCA
    print('Calculating PCA')
    pc.data = pc.data[intersect(var_genes, rownames(pc.data)),]
    pca.obj = rpca(t(pc.data), center=TRUE, scale=TRUE, retx=TRUE, k=num_pcs)
    seur@reductions[['pca']] = CreateDimReducObject(
        embeddings = pca.obj$x %*% diag(pca.obj$sdev**2),
	loadings = pca.obj$rotation,
	assay = 'RNA',
	stdev = pca.obj$sdev,
	key = 'PC_',
	misc = list(total.variance = sum(pca.obj$sdev))
    )
    
    # Calculate tSNE
    print('Calculating tSNE')
    tsne.rot = Rtsne(seur@reductions[['pca']]@cell.embeddings[,1:num_pcs], do.fast=TRUE, max_iter=max_iter, perplexity=25, verbose=T)$Y
    rownames(tsne.rot) = rownames(seur@reductions[['pca']]@cell.embeddings)
    seur@reductions[['tsne']] = CreateDimReducObject(
        embeddings = tsne.rot,
	assay = 'RNA',
	key = 'tSNE_'
    )
    
    # Cluster cells with Phenograph
    print('Clustering with Phenograph')
    clusters = run_phenograph(seur@reductions[['pca']]@cell.embeddings[,1:num_pcs], k=k)
    seur@meta.data[,colnames(clusters)] = clusters
    
    return(seur)
}


dirichlet_regression = function(counts, covariates, formula){

    # Dirichlet multinomial regression to detect changes in cell frequencies
    # formula is not quoted, example: counts ~ condition
    # counts is a [samples x cell types] matrix
    # covariates holds additional data to use in the regression
    #
    # Example:
    # counts = do.call(cbind, tapply(seur@data.info$orig.ident, seur@ident, table))
    # covariates = data.frame(condition=gsub('[12].*', '', rownames(counts)))
    # res = dirichlet_regression(counts, covariates, counts ~ condition)
    
    # Calculate regression
    counts = as.data.frame(counts)
    counts$counts = DR_data(counts)
    data = cbind(counts, covariates)
    fit = DirichReg(counts ~ condition, data)
    
    # Get p-values
    u = summary(fit)
    pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
    v = names(pvals)
    pvals = matrix(pvals, ncol=length(u$varnames))
    rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
    colnames(pvals) = u$varnames
    fit$pvals = pvals
    
    fit
}


