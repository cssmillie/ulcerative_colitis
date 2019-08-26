predict_log_base = function(x, n=10, bases.use=c(2, exp(1), 10), tol=1e-4){
    
    # Predict the log base for a vector of log-transformed integers
    # returns 'counts', 'tpm', or log base
    
    # Get values
    u = as.numeric(x[,1])
    v = as.numeric(x[,2])
    
    # Test for counts
    if(abs(sum(u - round(u))) <= tol & abs(sum(v - round(v))) <= tol){
        return('counts')
    }
    
    # Test for TPM
    if(abs(sum(u) - sum(v)) <= tol){
        return('tpm')
    }
    
    # Test for imputed TPM
    if(max(u) > 32 | max(v) > 32){
        print('predict_log_base: guessing data type = imputed TPM')
	stop('no good')
        return('tpm')
    }
    
    # Test for log(TPM)
    u = sapply(bases.use, function(base){
        abs(sum(base**u) - sum(base**v))
    })
    i = which.min(u)
    
    # If deviance is small, return log base
    if(u[[i]] <= tol){
        return(bases.use[[i]])
    } else {
        stop('predict_log_base: unknown data type')
    }
}


calc_tpm = function(seur=NULL, counts=NULL, tpm=NULL, data=NULL, genes.use=NULL, cells.use=NULL, total=1e4){
    require(wordspace)
    
    # Calculate TPM from @raw.data or data argument
    
    # Get counts data for calculation
    if(!is.null(counts)){data = counts}
    if(!is.null(tpm)){data = counts}
    if(is.null(data)){data = seur@assays[['RNA']]@counts}
    if(is.null(genes.use)){genes.use = rownames(data)}
    if(is.null(cells.use)){cells.use = colnames(data)}
        
    # Intersect genes.use and cells.use
    genes.use = intersect(genes.use, rownames(data))
    cells.use = intersect(cells.use, colnames(data))

    # Predict log base
    if(!is.null(counts)){base = 'counts'}
    if(!is.null(tpm)){base = 'tpm'}
    base = predict_log_base(data[, sample(1:ncol(data), 2)])
    
    if(base == 'counts'){
        numi = colSums(data[,cells.use])
	data = scaleMargins(data[genes.use, cells.use, drop=F], cols=total/numi)
    }
    if(base == 'tpm'){
        # Do nothing
    }
    if(is.numeric(base)){
        data = base**data[genes.use, cells.use, drop=F]
	data = data - min(data)
    }
    
    return(data)
}


