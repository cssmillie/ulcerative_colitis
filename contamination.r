
# -----------------------------------------------------------------
# Detect ambient RNA contamination in droplet-based sc-RNA-seq data
# -----------------------------------------------------------------


detect_contamination = function(tpm, groups, idents, samples, anno=NULL, fit.n=50, do.plot=TRUE, lab.use=NULL, prefix='test'){
    
    # Fit a contamination model to single cell expression data
    # --------------------------------------------------------
    # See Smillie, C.S., Biton, M.B., Ordovas, J.O., et al. (Cell, 2019) for more details
    #
    # Fits a contamination model using the specified "groups" (coarse-level clusters, e.g. Epithelial/Myeloid/Lymphoid/etc)
    # Then applies this model to each cell subset in "idents" (fine-level clusters, e.g. Enterocytes/Macrophages/Tregs/etc)
    #
    # Arguments:
    # - tpm = TPM expression matrix (rows = genes, columns = cell barcodes)
    # - groups = vector of "groups" each cell barcode is assigned to (coarse-level clusters, e.g. Epithelial/Myeloid/Lymphoid/etc)
    # - idents = vector of "idents" each cell barcode is assigned to (fine-level clusters, e.g. Enterocytes/Macrophages/Tregs/etc)
    # - samples = vector of samples that each cell barcode is associated with
    # - anno = named list mapping each "group" to its constituent "idents" (e.g. Myeloid -> DCs/Macrophages/Mast/etc)
    # - fit.n = number of genes to use for fit
    # - do.plot = plot regression? (TRUE/FALSE)
    # - lab.use = genes to label in plot
    # - prefix = output prefix (for plots)
    #
    # Returns object containing information about the fit, including the residuals (lower residuals = more likely contamination)
    
    # Fit models to cell groups
    print('Detecting contamination')
    
    print('Fitting model')
    res.groups = fit_contamination(tpm, idents=groups, samples=samples, anno=NULL, fit.n=fit.n, do.plot=do.plot, lab.use=lab.use, prefix=prefix)
    
    # Average model across groups
    coefs = sapply(res.groups, function(a) a$coefs)
    global_coefs = apply(coefs, 1, median)
    
    # Run average model on group
    print('Applying model to groups')
    res.groups = fit_contamination(tpm, groups, samples, anno=NULL, fit.n=fit.n, do.plot=do.plot, lab.use=lab.use, prefix=prefix, coefs=global_coefs)
    
    # Run average model on idents
    print('Applying model to idents')
    res.idents = fit_contamination(tpm, idents, samples, anno=anno, fit.n=fit.n, do.plot=do.plot, lab.use=lab.use, prefix=prefix, coefs=global_coefs)
    
    # Return data
    return(list(res.groups=res.groups, res.idents=res.idents))
}



fit_contamination = function(tpm, idents, samples, anno=NULL, coefs.use=NULL, fit.n=50, do.plot=TRUE, lab.use=NULL, prefix='out'){

    # Fit contamination model to single cell expression data
    # ------------------------------------------------------
    # See "detect_contamination" for arguments and additional information
    
    # initialize variables
    if(is.null(anno)){anno = structure(unique(as.character(idents)), names=unique(as.character(idents)))}
    if(any(! idents %in% anno)){stop('! idents %in% anno')}
    groups = names(anno)
    
    # summarize data
    cat('\n\nDetecting ambient contamination\n\n')
        
    # iterate over groups
    res = sapply(groups, function(group){
        print(group)
        flush.console()
		
        # output file
        out = paste(prefix, group, 'fit.pdf', sep='.')
    	
        # subset data
	i = idents %in% anno[[group]]
	j = idents %in% setdiff(idents, i)
		
	# sample frequencies
	f = table(as.factor(samples)[i])
	f = as.matrix(f/sum(f))
	
	# group mean
	u1 = rowMeans(tpm[,i])
	
	# other mean
	u2 = sapply(unique(samples), function(a){
	    rowSums(tpm[,j & (samples == a)])/sum(samples == a)
	})
	
	stopifnot(colnames(u2) == colnames(f))
	u2 = (u2 %*% f)[,1]
	
	# log-transform
	nice_log2 = function(x){y = log2(x); y[is.infinite(y)] = NA; y}
	l1 = nice_log2(u1)
	l2 = nice_log2(u2)
	
        # fit boundaries
        lo = quantile(l2[u1 == 0], .9, na.rm=T)
        hi = sort(l2, decreasing=T)[100]
        cat(paste0('\n\tLo Cutoff = ', lo, '\n\tHi Cutoff = ', hi))
        exclude = list(c(-Inf, lo), c(hi, Inf))
	
	# select points for regression
	lab.fit = names(select_points(l2, l1, n=fit.n, dir='down', nbins=10, loess=TRUE, exclude=exclude))
	cat(paste0('\n\tGenes for regression: ', paste(lab.fit, collapse=', ')))
	    
	# robust linear model
	cat('\n\tFitting rlm')
	fit = rlm(l1[lab.fit] ~ l2[lab.fit])	
	coefs = as.matrix(coef(fit))
	print(coefs)
	
	if(!is.null(coefs.use)){coefs = cbind(coefs.use, coefs)}
	
	# calculate residuals
	residuals = l1 - (coefs[2,1]*l2 + coefs[1,1])
	lab.con = names(which(residuals < 2))
	cat(paste0('\n\tLikely contaminants: ', paste(lab.con, collapse=', ')))
	if(do.plot == TRUE){plot_contamination(u1, u2, coefs, residuals, lab.use=lab.use, lab.fit=lab.fit, fit.cutoff=fit.cutoff, out=out, lab.n=20)}
	
	# update results
	list(u1=u1, u2=u2, fit=fit, coefs=coefs, residuals=residuals, lab.use=lab.use, lab.fit=lab.fit, lab.con=lab.con)

    }, simplify=F)
    
    names(res) = groups
    cat(paste('\ndetect_contamination: finished\n', names(res), '\n'))
    return(res)
}


loess_regression = function(...){
    # fit loess curve and get fitted values
    fit = loess(...)
    p.x = fit$x[order(fit$x)]
    p.y = fit$fitted[order(fit$x)]
    return(list(fit=fit, x=p.x, y=p.y))
}


select_points = function(x, y, n, dir='both', loess=FALSE, nbins=25, bin_type='equal_width', exclude=c(0,0)){
    
    # fix inputs
    if(!is.list(exclude)){exclude=list(exclude)}
    i = ((is.na(x) | is.na(y)) | (is.infinite(x) | is.infinite(y)))
    i = which(!i)
    xi = x[i]
    yi = y[i]
    
    # de-trend
    if(loess == TRUE){
        l = loess_regression(yi ~ xi, family='symmetric')
	yi = l$fit$residuals
    }
    
    # exclude data
    j = apply(sapply(exclude, function(e) (e[[1]] < xi) & (xi < e[[2]])), 1, any)
    j = which(!j)
    xi = xi[j]
    yi = yi[j]
    i = i[j]
    
    # bin x-axis
    if(bin_type == 'equal_width'){
        groups = cut2(xi, cuts=seq(from=min(xi, na.rm=T), to=max(xi, na.rm=T), length.out=nbins), m=2*n/nbins)
    } else {
        groups = cut2(xi, g=nbins)
    }
    
    # points
    j = c()
    
    # up points
    if(dir %in% c('up', 'both')){
        j = c(j, as.numeric(simple_downsample(cells=1:length(xi), groups=groups, ngene=yi, total_cells=n)))
    }
    
    # down points
    if(dir %in% c('down', 'both')){
        j = c(j, as.numeric(simple_downsample(cells=1:length(xi), groups=groups, ngene=-1*yi, total_cells=n)))
    }
    
    return(i[j])
}


plot_contamination = function(u1, u2, coefs, residuals, lab.use=NULL, lab.fit=NULL, lab.n=NULL, fit.cutoff=2, lab.name='Group', out=NULL){

    # Plot the regression model, labeling select genes and putative contaminants
    # --------------------------------------------------------------------------
    # u1 = mean TPM expression (in-group)
    # u2 = mean TPM expression (out-group)
    # coefs = coefficients of linear model
    # residuals = residuals for each gene
    # lab.use = additional genes to label
    # lab.fit = genes used for fit
    # lab.n = number of genes to label
    # fit.cutoff = residuals cutoff to use
    # out = output filename
    
    # log-transform
    l1 = log2(u1 + .5*min(u1[u1 > 0]))
    l2 = log2(u2 + .5*min(u2[u2 > 0]))
    
    # contamination
    lab.con = names(which(residuals < fit.cutoff))
    
    # scatterplot data
    d = data.frame(x=l2, y=l1, lab=ifelse(names(l1) %in% lab.fit, names(l1), ''), Type=rep('Other', length(l1)), stringsAsFactors=F)
    d[lab.con, 'Type'] = 'Contamination'
    d[lab.fit, 'Type'] = 'Fit'
    lab.use = intersect(rownames(d), lab.use)
    d[lab.use, 'Type'] = 'Label'
    d[lab.use, 'lab'] = lab.use
    
    # select subset to label
    if(!is.null(lab.n)){
        lab.sub = sample(lab.fit, min(lab.n, length(lab.fit)))
	d[(!d$lab %in% lab.use & !d$lab %in% lab.sub), 'lab'] = ''
    }
        
    # make plot
    if(!is.null(out)){alpha=.25} else {alpha=1}
    p = ggplot(d, aes(x=x, y=y)) +
        geom_point(aes(colour=Type)) +
   	geom_text_repel(aes(label=lab), size=3, segment.color='grey') +
	xlab('log2(<TPM>) (non-group)') +
    	ylab('log2(<TPM>) (group)') +
    	scale_colour_manual(values=c('lightcoral', 'black', 'steelblue3', 'lightgray')) +	
    	theme_cowplot()

    # add regression lines
    coefs = as.matrix(coefs)
    for(j in 1:ncol(coefs)){
        x0 = (min(l1, na.rm=T) - coefs[1,j] - fit.cutoff)/coefs[2,j]
	x1 = max(l2, na.rm=T)
	d.line = data.frame(x=c(x0, x1))
        d.line$y = coefs[2,j]*d.line$x + coefs[1,j] + fit.cutoff
	if(j == 1){lty = 'longdash'} else {lty = 'dotted'}
	p = p + geom_line(data=d.line, aes(x=x, y=y), lty=lty)
    }
    
    # save or display plot
    if(!is.null(out)){
        save_plot(p, file=out, nrow=2.25, ncol=2.5)
    } else {
        p
    }
}
