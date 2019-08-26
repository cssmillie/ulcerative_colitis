
load_mast = function(){
    library(BiocGenerics, pos=length(search()))
    library(S4Vectors, pos=length(search()))
    library(DelayedArray, pos=length(search()))
    library(MAST)
}


test_log_base = function(seur, base=2, total=1e4){
    # Test log base of seur@assays[['RNA']]@data
    j = sample(1:ncol(seur@assays[['RNA']]@data), 1)
    u = sum(base**seur@assays[['RNA']]@data[,j] - 1)
    (1e4 - 1e-2) <= u & u <= (1e4 + 1e-2)
}


unorder_factors = function(x){
    j = sapply(x, is.ordered)
    x[,j] = lapply(x[,j,drop=F], function(a) factor(as.character(a), levels=levels(a)))
    x
}


relevel_factors = function(x){
    j = sapply(x, is.factor)
    x[,j] = lapply(x[,j,drop=F], function(a) droplevels(a))
    j = sapply(x, function(a) length(unique(a)) > 1)
    x = x[,j,drop=F]
    x
}


get_data = function(seur, data.use='tpm', tpm=NULL, cells.use=NULL){
    
    # Retrieve data from a Seurat object
    # data.use can be: counts, tpm, log2, data, any matrix
    # can pre-calculate tpm for speed
    
    if(!is.character(data.use)){
        return(data.use)
    }
    
    if(data.use == 'counts'){
        data = seur@assays[['RNA']]@counts
    }
    
    if(data.use == 'tpm'){
        if(is.null(tpm)){
	    data = calc_tpm(seur, cells.use=cells.use)
	} else {
	    data = tpm
	}
    }
    
    if(data.use == 'log2'){
        if(!test_log_base(seur, base=2, total=1e4)){
	    stop('Error: seur@assays[["RNA"]]@data log base != 2')
        }
        data = seur@assays[['RNA']]@data
    }
    
    return(data)
}


select_cells = function(seur, covariates, batch.use=NULL, cells.use=NULL, max_cells=NULL){
    
    # Select cells from groups defined by covariates matrix
    # -----------------------------------------------------
    # 1. select cells.use
    # 2. build groups from covariates matrix
    # 3. select max_cells from each group
    # 4. sample evenly across batch.use
        
    cat('\nSelecting cells\n')
    
    # Get batch
    if(is.null(batch.use)){batch.use = rep('All', ncol(seur@assays[['RNA']]@data))}
    
    # Subset data
    i = apply(covariates, 1, function(a) !any(is.na(a)))
    covariates = covariates[i,,drop=F]
    batch.use = batch.use[i]
    
    # Get cells to use
    cells = rownames(covariates)
    if(!is.null(cells.use)){
        cells = intersect(cells, cells.use)
    }
    
    # Construct cell groups
    j = sapply(covariates, function(a) !is.numeric(a))
    groups = as.factor(apply(covariates[cells,j,drop=F], 1, function(a) paste(a, collapse='.')))
    
    # Select max_cells across all groups
    if(!is.null(max_cells)){
        
        # Combine batch.use and groups to sample evenly across batches
        batch.use = as.factor(paste0(groups, batch.use))

	# Sample [max_cells] evenly across all groups
	cells = simple_downsample(cells=cells, groups=batch.use, total_cells=max_cells)
	
	# Print subsampled group sizes
	groups = apply(covariates[cells,j,drop=F], 1, function(a) paste(a, collapse='.'))
    }
    return(cells)
}


select_genes = function(seur, stats, data.use=NULL, genes.use=NULL, min_cells=3, min_alpha=.025, min_fc=1.2, dir='both'){

    # Select genes from Seurat object by cells, genes.use, min_cells, min_alpha, and min_fc
    
    # Select genes by genes.use
    if(is.null(data.use)){data.use = seur@assays[['RNA']]@data}
    if(is.null(genes.use)){genes.use = rownames(data.use)}
    
    # Select genes by min_cells
    g1 = as.character(stats[, (max(n) >= min_cells) | (max(ref_n) >= min_cells), .(gene)][V1 == TRUE, gene])

    # Select genes by min_alpha
    g2 = as.character(stats[, (max(alpha) >= min_alpha) | (max(ref_alpha) >= min_alpha), .(gene)][V1 == TRUE, gene])
    
    # Select genes by min_fc
    if(dir == 'pos'){
        g3 = as.character(stats[, max(log2fc) >= log2(min_fc), .(gene)][V1 == TRUE, gene])
    } else {
        g3 = as.character(stats[, max(abs(log2fc)) >= log2(min_fc), .(gene)][V1 == TRUE, gene])
    }
    
    # Intersect and return genes.use
    genes.use = Reduce(intersect, list(genes.use, g1, g2, g3))
    return(genes.use)
}


expression_stats = function(tpm, covariates, formula, lrt_regex, genes.use=NULL, cells.use=NULL, invert_method='auto', invert_logic='last'){
    
    # Calculate expression statistics for groups specified by formula
    # each column of the model matrix is either a single term A or an interaction term A:B
    # if single term, then select cells A and ~A
    # if interaction, then select cells A:B and A:~B
    # for multi-level factors, ~A is the next highest level of A
    source('~/code/util/mm_utils.r')
    
    # Select cells
    if(!is.null(cells.use)){
        tpm = tpm[,cells.use]
	covariates = covariates[cells.use,,drop=F]
    }

    # Drop levels
    covariates = relevel_factors(covariates)

    # Select genes
    if(!is.null(genes.use)){tpm = tpm[genes.use,,drop=F]}
    total = rowSums(tpm)
    
    # Model matrix
    print('Model matrix')
    mm_formula = gsub('\\+ *\\S*\\|\\S+', '', as.character(formula), perl=T)
    mm = as.matrix(model.matrix(as.formula(mm_formula), data=unorder_factors(covariates)))
        
    # Invert matrix
    print('Invert matrix')
    u = mm_logical_not(mm, formula, covariates, method=invert_method, invert=invert_logic)
    MM = u$x
    refs = structure(u$names, names=colnames(MM))
        
    # For every column that matches lrt_regex
    print('Expression stats')
    stats = lapply(grep(lrt_regex, colnames(mm), value=T), function(a){print(a)
        
	# cell indices
	i = as.logical(mm[,a])
	j = as.logical(MM[,a])
	ref = refs[[a]]
		
	# number of expressing cells
	n1 = rowSums(tpm[,i,drop=F] > 0)
	n2 = rowSums(tpm[,j,drop=F] > 0)
	
	# total expression over all cells
	s1 = rowSums(tpm[,i,drop=F])
	s2 = rowSums(tpm[,j,drop=F])
	
	# fraction of expressing cells (alpha)
	a1 = n1/sum(i)
	a2 = n2/sum(j)
	
	# mean over expressing cells (mu)
	m1 = ifelse(n1 > 0, s1/n1, 0)
	m2 = ifelse(n2 > 0, s2/n2, 0)
	
	# mean over all cells
	u1 = s1/sum(i)
	u2 = s2/sum(j)
	
	# fraction of total expression
	t1 = s1/total
	t2 = s2/total
	
	# fix zeros for logs
	if(class(tpm) == 'dgCMatrix'){
	    zero = .5*min(tpm@x)/(sum(i) + sum(j)) # fast min (sparse matrix)
	} else {
	    zero = .5*min(tpm)/(sum(i) + sum(j)) # slow min (other matrix)
	}
	m1 = m1 + .5*zero
	m2 = m2 + .5*zero
	u1 = u1 + .5*zero
	u2 = u2 + .5*zero
	
	# log fold change
	log2fc = log2(u1) - log2(u2)
	
	# combine in data frame
	res = data.frame(gene=rownames(tpm), contrast=a, ref=ref, n=n1, ref_n=n2, alpha=a1, ref_alpha=a2, mu=log2(m1), ref_mu=log2(m2), mean=log2(u1), ref_mean=log2(u2), total=t1, ref_total=t2, log2fc=log2fc)
	return(res)
    })
    stats = as.data.table(do.call(rbind, stats))
    return(stats)
}


p.find_markers = function(seur, ident.1=NULL, ident.2=NULL, ident.use=NULL, tpm.use='tpm', data.use='log2', genes.use=NULL, cells.use=NULL, test.use='mast', min_cells=3, min_alpha=.05, min_fc=1.25,
                          max_cells=1000, batch.use=NULL, dir='pos', tpm=NULL, covariates=NULL, formula='~ ident', lrt_regex='ident', invert_method='auto', invert_logic='last', n.cores=1){
    
    # Get cell identities
    if(is.null(ident.use)){ident.use = seur@active.ident}
    
    # Build covariates
    print(c(ident.1, ident.2))
    if(!is.null(ident.1)){
	if(is.null(ident.2)){
	    ident.use = as.factor(ifelse(ident.use == ident.1, ident.1, 'Other'))
	    ident.use = relevel(ident.use, 'Other')
	} else {
	    ident.use = factor(ifelse(ident.use %in% c(ident.1, ident.2), as.character(ident.use), NA), levels=c(ident.2, ident.1))
	}
	if(is.null(covariates)){
	    covariates = data.frame(ident=ident.use)
	} else {
	    covariates$ident = ident.use
	}
    }
    rownames(covariates) = colnames(seur@assays[['RNA']]@data)
    
    # Check covariates
    q = sapply(covariates, typeof)
    if('character' %in% q){print(q); stop('error: invalid covariates type')}
    
    # Select cells
    cells.use = select_cells(seur, covariates, cells.use=cells.use, max_cells=max_cells, batch.use=batch.use)

    # TPM for log fold changes [genes x cells]
    tpm = get_data(seur, data.use=tpm.use, tpm=tpm, cells.use=cells.use)
    
    # Data for DE test [genes x cells]
    data = get_data(seur, data.use=data.use, tpm=tpm, cells.use=cells.use)
    
    # Calculate expression statistics
    lrt_regex = escapeRegex(lrt_regex)
    stats = expression_stats(tpm, covariates, formula, lrt_regex, genes.use=genes.use, cells.use=cells.use, invert_method=invert_method, invert_logic=invert_logic)

    # Select genes
    genes.use = select_genes(seur, stats, data.use=data, genes.use=genes.use, min_cells=min_cells, min_alpha=min_alpha, min_fc=min_fc, dir=dir)
    if(length(genes.use) == 0){return(c())}
    
    # Subset data
    print(paste('Testing', length(genes.use), 'genes in', length(cells.use), 'cells'))
    data = data[genes.use, cells.use, drop=F]
    data = rbind(data, rnorm(ncol(data)))
    covariates = unorder_factors(covariates[cells.use, , drop=F])
    covariates = relevel_factors(covariates)
    
    # Run marker tests
    labels = covariates[,1]
    markers = de.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, n.cores=n.cores)
    
    # Add cluster information
    if(! 'contrast' %in% colnames(markers)){
        markers$contrast = paste0('ident', levels(labels)[nlevels(labels)])
    }
    
    # Merge results
    markers = as.data.table(markers)
    setkey(markers, gene, contrast)
    if(!is.null(stats)){
        setkey(stats, gene, contrast)
        markers = markers[stats,]    
    }
    
    # Sort markers
    markers = markers[order(contrast, pvalH),]
    
    # Add ident
    markers$ident = paste(ident.1, ident.2, sep=';')
    markers$ident = gsub(';$', '', markers$ident)
    
    # Return marker genes
    return(markers)
}


p.find_all_markers = function(seur, ident.use=NULL, data.use='log2', tpm=NULL, do.precalc=T, n.cores=1, ...){
    
    # Get cell identities
    if(is.null(ident.use)){ident.use = seur@active.ident}
    ident.use = as.factor(ident.use)
    idents = as.character(levels(ident.use))
        
    # Pre-calculate TPM and data
    if(do.precalc == TRUE){
        tpm = get_data(seur, data.use='tpm', tpm=tpm)
        data.use = get_data(seur, data.use=data.use, tpm=tpm)
    }
    
    # Find marker genes
    run_parallel(
	foreach(i=idents, .combine=rbind) %dopar% {
	    print(i)
	    p.find_markers(seur, ident.1=i, ident.use=ident.use, tpm=tpm, data.use=data.use, n.cores=1, ...)
	},
	n.cores = n.cores
    )
}


de.mast = function(data, covariates, formula=NULL, lrt_regex=TRUE, n.cores=1){
    
    load_mast()
    options(mc.cores=n.cores)
        
    # Make single cell assay (SCA) object
    fdata = data.frame(matrix(rep(1, nrow(data))))
    covariates = as.data.frame(covariates)
    sca = MAST::FromMatrix(as.matrix(data), covariates, fdata)    
    
    # Fit MAST hurdle model
    if(is.null(formula)){
        formula = paste('~' , paste(colnames(covariates), collapse=' + '))
    }
    formula = as.formula(formula)
    zlm.obj = zlm(formula, sca, force=TRUE)
    
    # Likelihood ratio test
    if(is.logical(lrt_regex)){lrt_regex = colnames(covariates)}
    contrasts = grep(paste(lrt_regex, collapse='|'), colnames(zlm.obj@coefC), value=T, perl=T)
    res = summary(zlm.obj, doLRT=contrasts)$datatable
    
    # Get component information
    res.f = res[res$component == 'logFC', .(primerid, contrast, coef)]
    res.d = res[res$component == 'D', .(primerid, contrast, coef, `Pr(>Chisq)`)]
    res.c = res[res$component == 'C', .(primerid, contrast, coef, `Pr(>Chisq)`)]
    res.h = res[res$component == 'H', .(primerid, contrast, `Pr(>Chisq)`)]
    
    # Combine results
    res = merge(res.d, res.c, by=c('primerid', 'contrast'), all=T, suffixes=c('D', 'C'))
    res = Reduce(function(...) merge(..., by=c('primerid', 'contrast'), all=T), list(res, res.f, res.h))
    res = data.frame(subset(res, !is.na(`Pr(>Chisq)`)), stringsAsFactors=F)
    
    # Cleanup results
    colnames(res) = c('gene', 'contrast', 'coefD', 'pvalD', 'coefC', 'pvalC', 'mastfc', 'pvalH')
    res = res[res$gene != '',]
    res = res[order(res$contrast, res$pvalH),]
    
    # Replace NAs in mastfc with maximum
    res = as.data.table(res)
    res[, mastfc := ifelse(is.na(mastfc), max(mastfc, na.rm=T), mastfc), .(contrast)]
    
    # Adjust p-values
    res[, padjD := p.adjust(pvalD, 'fdr'), .(contrast)]
    res[, padjC := p.adjust(pvalC, 'fdr'), .(contrast)]
    res[, padjH := p.adjust(pvalH, 'fdr'), .(contrast)]
    
    options(mc.cores=1)
    return(res)
}

