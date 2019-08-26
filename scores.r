

predict_dge_type = function(x, bases.use=c(2, exp(1), 10), tol=1e-4){
    # Returns counts, tpm, or log base
    
    # Select data
    u = as.numeric(x[,1])
    v = as.numeric(x[,2])
    
    # Check for integers
    if(abs(sum(u - round(u))) <= tol & abs(sum(v - round(v))) <= tol){
        return('counts')
    }
    
    # Check for constant sum
    if(abs(sum(u) - sum(v)) <= tol){
        return('tpm')
    }
    
    # Check for big numbers
    if(max(u) > 32 | max(v) > 32){
        print('predict_log_base: found big numbers, so guessing data type = imputed TPM')
        return('tpm')
    }

    # Check for logTPM
    u = sapply(bases.use, function(base){
        abs(sum(base**u) - sum(base**v))
    })
    i = which.min(u)
    
    # Test logTPM deviance
    if(u[[i]] <= tol){
        return(bases.use[[i]])
    } else {
        print('predict_log_base: confused, so guessing data type = TPM???')
    }
}


get_data = function(seur, data.use='tpm', tpm=NULL, genes.use=NULL, cells.use=NULL){
    
    # Retrieve data from a Seurat object
    # data.use can be: counts, tpm, log2, data, any matrix
    # optionally, pre-calculate tpm for speed
    
    if(!is.character(data.use)){
        data = data.use
    } else if(data.use == 'counts'){
        data = seur@assays[['RNA']]@counts
    } else if(data.use == 'tpm'){
        if(is.null(tpm)){
	    data = calc_tpm(seur, genes.use=genes.use, cells.use=cells.use)
	} else {
            data = tpm
        }
    } else if(data.use == 'log2'){
        if(predict_dge_type(seur@assays[['RNA']]@data, bases.use=c(2)) != 2){
	    print('Warning: seur@assays[["RNA"]]@data is not log base 2')
	}
	data = seur@assays[['RNA']]@data
    } else if(data.use == 'scale'){
        if(!is.null(genes.use)){
	    genes.use = intersect(genes.use, rownames(seur@assays[['RNA']]@data))
	} else {
	    genes.use = rownames(seur@assays[['RNA']]@data)
	}
	data = t(scale(t(seur@assays[['RNA']]@data[genes.use,])))
    } else {
        stop('Error: get_data invalid argument')
    }
    
    # Subset genes and cells
    if(!is.null(genes.use)){data = data[intersect(genes.use, rownames(data)),]}
    if(!is.null(cells.use)){data = data[,intersect(cells.use, colnames(data))]}

    return(data)
}


map_names = function(seur=NULL, names=NULL){
    
    # Map "names" to genes or features in a Seurat object
    # ---------------------------------------------------
    # seur = seurat object
    # names = list of names (genes or features)
    # returns list(genes=c(GENES), feats=(FEATS))
    
    # Initialize variables
    names = as.list(names)
    genes = c()
    feats = c()
    
    # Get data and metadata
    data = seur@assays[['RNA']]@data
    meta = seur@meta.data
    
    # Map names
    if(!is.null(names)){
	genes = sapply(names, function(a){intersect(a, rownames(data))}, simplify=F)
	feats = sapply(names, function(a){intersect(a, colnames(meta))}, simplify=F)
    }
    
    # Filter genes and feats
    genes = genes[lengths(genes) > 0]
    feats = feats[lengths(feats) > 0]
        
    # Fix gene names
    if(length(genes) > 0){
        if(is.null(names(genes))){names(genes) = sapply(genes, paste, collapse='.')}
    }
    
    # Fix feat names
    if(length(feats) > 0){
        if(is.null(names(feats))){names(feats) = sapply(feats, paste, collapse='.')}
    }
    
    return(list(genes=genes, feats=feats))
}


score_cells = function(seur=NULL, names=NULL, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL){
    
    # Score genes and features across cells and optionally aggregate
    # The steps are:
    # - calculate mean expression across genes (combine_genes = 'sum', 'mean', 'scale', 'scale2')
    # - calculate mean expression within each cell type (groups, group_stat = 'mean', 'alpha', 'mu')
    
    # Fix input arguments and get data for name mapping
    data = seur@assays[['RNA']]@data
    meta = seur@meta.data
    if(!is.null(groups)){groups = setNames(groups, colnames(seur@assays[['RNA']]@data))}
    scores = NULL
    
    # Map genes and feats
    res = map_names(seur=seur, names=names)
    genes = res$genes
    feats = res$feats
    genes.use = unique(do.call(c, genes))
    
    # Subset cells
    if(!is.null(cells.use)){
        data = data[,cells.use,drop=F]
        meta = meta[cells.use,,drop=F]
        if(!is.null(groups)){groups = groups[cells.use]}
    }
        
    group_genes = function(x, method){
        
        # combine expression data across genes within a signature
	# x = [genes x cells] matrix
	# method = 'sum', 'mean', 'scale'
	# returns [genes x cells] or [1 x cells] matrix
	
	if(nrow(x) == 1){return(x[1,,drop=F])}
	if(method == 'sum'){
	    t(colSums(x))
	} else if(method == 'mean'){
	    t(colMeans(x, na.rm=T))
	} else if(method == 'scale'){
	    x = t(scale(t(x)))
	    t(colMeans(x, na.rm=T))
	} else if(method == 'scale2'){
	    x = t(scale(t(x), center=F))
	    t(colMeans(x, na.rm=T))
	} else if(method == 'none'){
	    x
	} else {
	    stop('Error: invalid combine_genes method')
	}
    }
    
    group_cells = function(x, groups, method){
        
        # combine expression data across cells
	# x = [genes x cells] matrix
	# group_stat = 'alpha', 'mu', or 'mean'
	# returns [genes x groups] matrix
	
	if(is.null(groups)){return(x)}
	if(method %in% c('n', 'sum')){
	    if(method == 'n'){x = x > 0}
	    x = t(data.frame(aggregate(t(x), list(groups), sum, na.rm=T), row.names=1))
	} else {
	    if(method == 'alpha'){x = x > 0}
	    if(method == 'mu'){x[x == 0] = NA}
	    x = t(data.frame(aggregate(t(x), list(groups), mean, na.rm=T), row.names=1))
	}
	x[is.na(x)] = 0
	x
    }
    
    # Calculate scores
    names.use = unique(c(names(genes), names(feats)))
    
    # Speed improvements (fast indexing for flat structures)
    name_map = sapply(names.use, function(a) c(genes[[a]], feats[[a]]), simplify=F)
    do.flat = all(lengths(name_map) == 1)
    if(do.flat == TRUE){
        genes[['flat']] = do.call(c, genes)
	feats[['flat']] = do.call(c, feats)
	names.iter = 'flat'
	combine_genes = 'none'
    } else {
        names.iter = names.use
    }
    
    backup = scores
    scores = lapply(names.iter, function(name){
        
        # Combine data and metadata
	if(name %in% names(genes)){
	    si = data[genes[[name]],,drop=F]
	} else {
	    si = c()
	}
	if(name %in% names(feats)){
	    if(is.null(si)){
	        si = t(meta[,feats[[name]],drop=F])
	    } else {
	        si = rBind(si, t(meta[,feats[[name]],drop=F]))
	    }
	}
	si = as.matrix(as.data.frame(si))
	si = group_genes(si, method=combine_genes)
	si = group_cells(si, groups=groups, method=group_stat)
	si = data.frame(t(si))
    })
    
    # Collapse scores
    if(do.flat == TRUE){
	scores = scores[[1]][,make.names(name_map[names.use]),drop=F]
    } else {
        do.collapse = all(lapply(scores, ncol) == 1)
	if(do.collapse == TRUE){
	    scores = as.data.frame(do.call(cbind, scores))
	}
    }
    
    # Fix names
    names.use = make.names(names.use)
    names(scores) = names.use
    
    # Combine data
    if(!is.null(backup)){
        if(is.data.frame(scores)){
            scores = cbind.data.frame(scores, backup)
        } else {
            scores = c(scores, backup)
        }
    }

    if(nrow(scores) == 0){scores = NULL}
    return(scores)
}
