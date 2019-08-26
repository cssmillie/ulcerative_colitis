

qtrim = function(x, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, rescale=NULL){
    
    # Trim by value
    x[x < vmin] = vmin
    x[x > vmax] = vmax
    
    # Trim by quantile
    u = quantile(x, qmin, na.rm=T)
    v = quantile(x, qmax, na.rm=T)
    x[x < u] = u
    x[x > v] = v

    return(x)
}


plot_tsne = function(seur=NULL, names=NULL, scores=NULL, coords=NULL, data=NULL, meta=NULL, regex=NULL, files=NULL, file.cols=NULL, file.regex=NULL, top=NULL, ident=TRUE, data.use='log2',
                     combine_genes='mean', cells.use=NULL, ymin=0, ymax=1, num_col='auto', pt.size=.75, font.size=11, do.label=T, label.size=5, do.title=TRUE, title.use=NULL, pal=NULL,
	             do.legend=TRUE, legend.title='log2(TP10K+1)', share_legend=FALSE, legend_width=.05, legend.cols=NULL, vmin=NA, vmax=NA, na.value='transparent', out=NULL, nrow=1.5, ncol=1.5, ...){
    
    
    # TSNE coordinates
    if(is.null(coords)){
        d = structure(as.data.frame(seur@reductions[['tsne']]@cell.embeddings[,1:2]), names=c('x', 'y'))
    } else {
        d = structure(as.data.frame(coords[,1:2]), names=c('x', 'y'))
    }
    
    # Cell identities
    if(!is.logical(ident)){
        d$Identity = ident
    } else if(ident & !is.null(seur)){
        d$Identity = seur@active.ident
    }

    # Cell scores
    if(!is.null(cells.use)){cells.use = intersect(cells.use, rownames(d))}
    scores = score_cells(seur=seur, names=names, combine_genes=combine_genes, cells.use=cells.use)
    if(!is.null(scores)){
        i = intersect(rownames(d), rownames(scores))
	ni = c(names(d), names(scores))
	d = cbind.data.frame(d[i,], scores[i,])
	colnames(d) = ni
    }
        
    # Subset cells
    if(is.null(cells.use)){cells.use = rownames(d)}
    d = data.frame(d[cells.use,])
        
    # Initialize plotlist
    cat('\nPlotting:', paste(colnames(subset(d, select=-c(x,y))), collapse=', '), '\n')
    ps = list()

    # Shuffle point order
    d = d[sample(1:nrow(d)),]
    
    # Get limits for shared legend
    if(share_legend == TRUE){
        j = names(which(sapply(subset(d, select=-c(x,y)), is.numeric)))
	cat('\nShared limits:', paste(j, collapse=', '), '\n')
        if(is.na(vmin)){vmin = na.omit(min(d[,j]))}
	if(is.na(vmax)){vmax = na.omit(max(d[,j]))}
	cat('> vmin =', vmin, '\n> vmax =', vmax, '\n')
    }
    
    for(col in setdiff(colnames(d), c('x', 'y'))){

        # plot NAs first
        d = d[c(which(is.na(d[,col])), which(!is.na(d[,col]))),]
	    	
	if(is.numeric(d[,col])){
		    
	    # Continuous plot
	    d[,col] = qtrim(d[,col], qmin=ymin, qmax=ymax)
	    p = ggplot(d) +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		scale_colour_gradientn(colours=material.heat(50), guide=guide_colourbar(barwidth=.5, title=legend.title), na.value=na.value, limits=c(vmin, vmax)) + 
		theme_cowplot(font_size=font.size) +
		xlab('TSNE 1') + ylab('TSNE 2')
		
	} else {
	    
	    # Discrete plot

	    # Get colors
	    if(do.label == TRUE){tsne.colors = set.colors} else {tsne.colors = set.colors}
	    if(!is.null(pal)){tsne.colors = pal}
	    
	    p = ggplot(d) +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		theme_cowplot(font_size=font.size) +
		xlab('TSNE 1') + ylab('TSNE 2') +
		scale_colour_manual(values=tsne.colors, na.value=na.value, drop=F) + 
		theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
		guides(colour=guide_legend(ncol=legend.cols, title=legend.title))
	    
	    if(do.label == T){
	        t = aggregate(d[,c('x', 'y')], list(d[,col]), median)
		colnames(t) = c('l', 'x', 'y')
		p = p + geom_text_repel(data=t, aes(x=x, y=y, label=l, lineheight=.8), point.padding=NA, size=label.size, family='Helvetica') + theme(legend.position='none')
	    }
	}
	
	if(do.title == TRUE){
	    if(is.null(title.use)){title = col} else {title = title.use}
	    p = p + ggtitle(title)
	}
	if(do.legend == FALSE){p = p + theme(legend.position='none')}
	if(legend.title == ''){p = p + theme(legend.title=element_blank())}

	ps[[col]] = p
    }
    
    if(num_col == 'auto'){
        num_col = ceiling(sqrt(length(ps)))
    }
    
    ps = make_compact(plotlist=ps, num_col=num_col)
    if(length(ps) > 1){
        if(share_legend == TRUE){
            p = share_legend(ps, num_col=num_col, width=legend_width)
        } else {
            p = plot_grid(plotlist=ps, ncol=num_col, align='h')
        }
    }
    
    if(is.null(out)){
        p
    } else {
        save_plot(file=out, p, nrow=nrow, ncol=ncol)
    }
}


make_compact = function(plotlist, num_col, labels=TRUE, ticks=TRUE){
    
    # Make a "plot_grid" plotlist compact by removing axes from interior plots

    # x-axis
    if(length(plotlist) > num_col){
        i = setdiff(1:length(plotlist), rev(1:length(plotlist))[1:min(num_col, length(plotlist))])
	if(labels == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.title.x=element_blank()))}
	if(ticks == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.x=element_blank()))}
    }
    
    # y-axis
    if(num_col > 1){
        i = setdiff(1:length(plotlist), which(1:length(plotlist) %% num_col == 1))
	if(labels == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.title.y=element_blank()))}
	if(ticks == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.y=element_blank()))}
    }
    
    return(plotlist)
}


share_legend = function(plotlist, num_col, rel_widths=NULL, width=.1){
    
    # Get first legend in plotlist
    i = min(which(sapply(plotlist, function(p) 'guide-box' %in% ggplotGrob(p)$layout$name)))
    cat(paste('\nUsing shared legend:', names(plotlist)[i], '\n'))
    legend = get_legend(plotlist[[i]])
    
    # Remove all legends
    plotlist = lapply(plotlist, function(p) p + theme(legend.position='none'))
    
    # Make combined plot
    if(is.null(rel_widths)){rel_widths = rep(1, length(plotlist))}
    p = plot_grid(plotlist=plotlist, ncol=num_col, align='h', rel_widths=rel_widths)
    p = plot_grid(p, legend, ncol=2, rel_widths=c(1-width, width))
    
    return(p)
}


ggheatmap = function(data, Rowv='hclust', Colv='hclust', xlab='', ylab='', xsec=FALSE, ysec=FALSE, xstag=FALSE, xstag_space=.15, ystag=FALSE, ystag_space=.15, title='', legend.title='',
                     pal='nmf', do.legend=TRUE, font_size=7, 
                     out=NULL, nrow=1.25, ncol=1.25, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, symm=FALSE, xstrip=NULL, ystrip=NULL, hclust_met='complete', border='#cccccc', replace_na=NA,
		     labRow=NULL, labCol=NULL, ret.order=FALSE, ret.legend=FALSE, pvals=NULL, max_pval=1, pval_border='black', pval_width=.25){
    
    # Scale values
    data[is.na(data)] = replace_na
    data = qtrim(data, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax)
    
    # Convert to long format
    x = as.data.frame(data) %>% rownames_to_column('row') %>% gather(col, value, -row)
    x$value = as.numeric(x$value)
    
    # Merge data and p-values
    if(is.null(pvals)){
        x$pval = border
    } else {
        if(ncol(pvals) == 3){
	    colnames(pvals) = c('row', 'col', 'pval')
	    pvals$row = as.character(pvals$row)
	    pvals$col = as.character(pvals$col)
	} else {
	    pvals = as.data.frame(pvals) %>% rownames_to_column('row') %>% gather(col, pval, -row)
	    pvals$pval = as.numeric(pvals$pval)
	}
	if(length(intersect(x$row, pvals$row)) == 0){
	    colnames(pvals) = c('col', 'row', 'pval')
	}
	x = as.data.frame(merge(as.data.table(x), as.data.table(pvals), by=c('row', 'col'), all.x=TRUE))
	x$pval = ifelse(x$pval <= max_pval, pval_border, NA)
	x$pval[is.na(x$pval)] = border
    }
    
    # Order rows
    if(length(Rowv) > 1){rowv = Rowv; Rowv = 'Rowv'} else {rowv = rev(rownames(data))}
    if(length(Colv) > 1){colv = Colv; Colv = 'Colv'} else {colv = colnames(data)}
    if(nrow(data) <= 2){Rowv = 'none'}
    if(ncol(data) <= 2){Colv = 'none'}
    if(Rowv == 'hclust'){
        rowv = rev(rownames(data)[hclust(dist(data), method=hclust_met)$order])
    }
    if(Colv == 'hclust'){
        colv = colnames(data)[hclust(dist(t(data)), method=hclust_met)$order]
    }
    if(Rowv == 'none'){
        rowv = rev(rownames(data))
    }
    if(Colv == 'none'){
        colv = colnames(data)
    }
    if(Rowv == 'min'){
        rowv = rev(rownames(data)[order(apply(data, 1, which.min))])
    }
    if(Rowv == 'max'){
    	i = match(colnames(data)[apply(data, 1, which.max)], colv)
	rowv = rev(rownames(data)[order(i)])
    }
    if(Colv == 'min'){
        i = match(rownames(data)[apply(data, 2, which.min)], rowv)
	colv = rev(colnames(data)[order(i)])        
    }
    if(Colv == 'max'){
	i = match(rownames(data)[apply(data, 2, which.max)], rowv)
	colv = rev(colnames(data)[order(i)])
    }
    Rowv = rowv
    Colv = colv

    # Set order of row and column labels
    x$row = factor(x$row, levels=Rowv)
    x$col = factor(x$col, levels=Colv)
    
    # Get odd/even indices
    r1 = seq(1, length(Rowv), by=2)
    r2 = seq(2, length(Rowv), by=2)
    c1 = seq(1, length(Colv), by=2)
    c2 = seq(2, length(Colv), by=2)
    
    # Get plot data
    if(length(pal)==1){if(pal == 'nmf'){pal = rev(colorRampPalette(nmf.colors)(101))[10:101]} else {pal = colorRampPalette(brewer.pal(9, pal))(101)}}
    
    # Plot significant boxes last
    x = x[rev(order(x$pval != 'black')),]
        
    # Plot with geom_tile
    p = ggplot(x) +
        geom_tile(aes(x=as.numeric(col), y=as.numeric(row), fill=value), color=x$pval, size=pval_width) +
	labs(x=xlab, y=ylab, title=title, fill=legend.title) +
	theme_cowplot(font_size=font_size) +
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.line=element_blank())
    
    # Set scale
    if(is.infinite(vmin)){vmin = min(x$value)}
    if(is.infinite(vmax)){vmax = max(x$value)}
    limits = c(vmin, vmax)
    if(symm == TRUE){
        values = c(min(x$value), 0, max(x$value))
	p = p + scale_fill_gradientn(colours=pal, values=scales::rescale(values))
    } else {
        p = p + scale_fill_gradientn(colours=pal, limits=limits)
    }
    
    # Secondary x-axis
    if(xsec == FALSE){
        p = p + scale_x_continuous(breaks=1:length(Colv), labels=Colv, expand=c(0,0))
    } else {
	p = p + scale_x_continuous(breaks=c1, labels=Colv[c1], sec.axis=dup_axis(breaks=c2, labels=Colv[c2]), expand=c(0,0)) + theme(axis.text.x.top= element_text(angle=90, hjust=0, vjust=.5))
    }
    
    # Secondary y-axis
    if(ysec == FALSE){
        p = p + scale_y_continuous(breaks=1:length(Rowv), labels=Rowv, expand=c(0,0))
    } else {
	p = p + scale_y_continuous(breaks=r1, labels=Rowv[r1], sec.axis=dup_axis(breaks=r2, labels=Rowv[r2]), expand=c(0,0))
    }
    
    # Add axes with ggrepel (use both sides to fit as many labels as possible)
    if(!is.null(labRow)){
        p = p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	if(!is.logical(labRow)){
	    lab.use = intersect(Rowv, labRow)
	    lab.l = ifelse(Rowv %in% lab.use[seq(from=1, to=length(lab.use), by=2)], Rowv, '')
	    lab.r = ifelse(Rowv %in% lab.use[seq(from=2, to=length(lab.use), by=2)], Rowv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none', plot.margin=margin(l=-10, unit='pt'))
	    axis.l = add_ggrepel_axis(plot=p, lab=lab.l, type='left', font_size=font_size)
	    axis.r = add_ggrepel_axis(plot=p, lab=lab.r, type='right', font_size=font_size)
	    p = plot_grid(axis.l, p, axis.r, nrow=1, rel_widths=c(.15,1,.15), align='h', axis='tb')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))	    
	}
    }
    
    if(!is.null(labCol)){
        p = p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	if(xsec == TRUE & !is.logical(labCol)){
	    lab.use = intersect(Colv, labCol)
	    lab.u = ifelse(Colv %in% lab.use[seq(from=1, to=length(lab.use), by=2)], Colv, '')
	    lab.d = ifelse(Colv %in% lab.use[seq(from=2, to=length(lab.use), by=2)], Colv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none', plot.margin=margin(t=-10, b=-12.5, unit='pt'))
	    axis.u = add_ggrepel_axis(plot=p, lab=lab.u, type='up', font_size=font_size, do.combine=FALSE)
	    axis.d = add_ggrepel_axis(plot=p, lab=lab.d, type='down', font_size=font_size, do.combine=FALSE)
	    p = plot_grid(axis.u, p, axis.d, nrow=3, rel_heights=c(.15,1,.15), align='v', axis='lr')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))
	}
	if(xsec == FALSE & !is.logical(labCol)){
	    lab.use = ifelse(Colv %in% intersect(Colv, labCol), Colv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none')
	    axis.d = add_ggrepel_axis(plot=p, lab=lab.use, type='down', font_size=font_size, do.combine=FALSE)
	    p = plot_grid(p, axis.d, nrow=2, rel_heights=c(1, .15), align='v', axis='lr')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))
	}
    }
    
    if(xstag == TRUE){
        nudge = rep(c(0, 1), length.out=length(Colv))
	p = add_ggrepel_axis(plot=p, lab=Colv, dir='x', type='down', font_size=font_size, force=0, axis.width=xstag_space, nudge=nudge, ret.legend=ret.legend)
    }
    
    if(ystag == TRUE){
        nudge = rep(c(0, 1), length.out=length(Rowv))
	p = add_ggrepel_axis(plot=p, lab=Rowv, dir='y', type='left', font_size=font_size, force=0, axis.width=ystag_space, nudge=nudge, ret.legend=ret.legend)
    }
        
    if(do.legend == FALSE){p = p + theme(legend.position='none')}
    
    # Save plot
    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }
    
    if(ret.order == FALSE){p} else {list(p=p, Rowv=Rowv, Colv=Colv)}
}


add_ggrepel_axis = function(plot, lab, dir='both', type='left', lab.pos=NULL, lab.lim=NULL, nudge=0, force=1, axis.width=.10, legend.width=.075, font_size=8, ret.legend=FALSE, do.combine=TRUE){
    
    # Fix input arguments
    if(is.null(lab.pos)){lab.pos = 1:length(lab)}
    if(is.null(lab.lim)){lab.lim = c(min(lab.pos)-1, max(lab.pos)+1)}
    
    # Set label directions
    if(type == 'left'){x=0; y=lab.pos; xlim=c(-1,0); ylim=lab.lim; angle=0; nudge_x=-1*nudge; nudge_y=0}
    if(type == 'up'){x=lab.pos; y=0; xlim=lab.lim; ylim=c(0,1); angle=90; nudge_x=0; nudge_y=nudge}
    if(type == 'right'){x=0; y=lab.pos; xlim=c(0,1); ylim=lab.lim; angle=0; nudge_x=nudge; nudge_y=0}
    if(type == 'down'){x=lab.pos; y=0; xlim=lab.lim; ylim=c(-1,0); angle=90; nudge_x=0; nudge_y=-1*nudge}
    
    # Get data for ggplot
    d = data.frame(x=x, y=y, lab=lab, dir=dir, nudge_y=nudge_y)
    
    # Make ggrepel axis
    axis = ggplot(d, aes(x=x, y=y, label=lab)) +
           geom_text_repel(min.segment.length=grid::unit(0,'pt'), color='grey30', size=(font_size-1)/.pt, angle=angle, segment.color='#cccccc', segment.size=.15,
	                   direction=dir, nudge_x=nudge_x, nudge_y=nudge_y, force=force) +
	   scale_x_continuous(limits=xlim, expand=c(0,0), breaks=NULL, labels=NULL, name=NULL) +
	   scale_y_continuous(limits=ylim, expand=c(0,0), breaks=NULL, labels=NULL, name=NULL) +
    	   theme(panel.background = element_blank(), plot.margin = margin(0, 0, 0, 0, 'pt'))
    
    if(do.combine == FALSE){return(axis)}
    
    # Get plot legend
    legend = get_legend(plot)
    plot = plot + theme(legend.position='none')
    
    # Combine plots
    if(type == 'left'){
        plot = plot + scale_y_continuous(limits=lab.lim, expand=c(0,0))    
        plot = plot + theme(plot.margin=margin(l=-12.5, unit='pt')) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	p = plot_grid(axis, plot, nrow=1, align='h', axis='tb', rel_widths=c(axis.width, 1-axis.width))
	if(ret.legend == FALSE){p = plot_grid(legend, p, nrow=1, rel_widths=c(legend.width, 1-legend.width))}
    }
    if(type == 'up'){
        plot = plot + scale_x_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(t=-10, unit='pt')) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	p = plot_grid(axis, plot, nrow=2, align='v', axis='lr', rel_heights=c(axis.width, 1-axis.width))
	if(ret.legend == FALSE){p = plot_grid(legend, p, nrow=1, rel_widths=c(legend.width, 1-legend.width))}
    }
    if(type == 'right'){
        plot = plot + scale_y_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(r=-10, unit='pt')) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	p = plot_grid(plot, axis, nrow=1, align='h', axis='tb', rel_widths=c(1-axis.width, axis.width))
	if(ret.legend == FALSE){p = plot_grid(p, legend, nrow=1, rel_widths=c(1-legend.width, legend.width))}
    }
    if(type == 'down'){
        plot = plot + scale_x_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(b=-11.5, t=5, unit='pt')) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	p = plot_grid(plot, axis, nrow=2, align='v', axis='lr', rel_heights=c(1-axis.width, axis.width))
	if(ret.legend == FALSE){p = plot_grid(p, legend, nrow=1, rel_widths=c(1-legend.width, legend.width))}
    }
    
    if(ret.legend == FALSE){p} else {list(p=p, legend=legend)}
}


plot_dots = function(x1=NULL, y1=NULL, x2=NULL, y2=NULL, de=NULL, big_first=FALSE, reorder=FALSE, replace_na=0, na.value='transparent', symm=FALSE, vmin=NULL, vmax=NULL, qmin=0, qmax=1, max_size=5,
                     fill_title='', size_title='', Rowv=NULL, Colv=NULL, pal=rev(brewer.pal(7, 'RdBu')), out=NULL, nrow=1, ncol=1, coord_flip=FALSE, sig_only=TRUE, xlab='', ylab='', font_size=10,
		     pvals=NULL){
    
    # input arguments
    # x1 = size matrix or list(row, col, size)
    # x2 = size matrix or list(row, col, size)
    # y1 = color matrix or list(row, col, color)
    # y2 = color matrix or list(row, col, color)
    # xanno = annotation track (named list)
    # yanno = annotation track (named list)
    # pvals = data.frame(row, col)
    
    # convert to long format
    if(!is.null(de)){
        de = as.data.table(de)
        if(sig_only == TRUE){de = de[padjD < .05]}
	x1 = de[,.(col=ident, row=gene, size=alpha)]
	x2 = de[,.(col=ident, row=gene, size=ref_alpha)]
	y1 = de[,.(col=ident, row=gene, color=coefD)]
	y2 = de[,.(col=ident, row=gene, color=coefD)]
    }
    if(!is.null(x1)){
        if(! all(c('row', 'col', 'size') %in% colnames(x1))){
	    cnames = rev(colnames(x1))
	    x1 = x1 %>% rownames_to_column('row')
	    x1$row = factor(x1$row, levels=x1$row)
	    x1 = x1 %>% gather(col, size, -row)
	    x1$col = factor(x1$col, levels=cnames)	    
	}
        x1 = as.data.table(x1)
    }
    if(!is.null(y1)){
        if(! all(c('row', 'col', 'color') %in% colnames(y1))){
	    cnames = rev(colnames(y1))
	    y1 = y1 %>% rownames_to_column('row')
	    y1$row = factor(y1$row, levels=y1$row)
	    y1 = y1 %>% gather(col, color, -row)
	    y1$col = factor(y1$col, levels=cnames)	    
	}
	y1 = as.data.table(y1)	
    }
    if(!is.null(x2)){
        if(! all(c('row', 'col', 'size') %in% colnames(x2))){
	    cnames = rev(colnames(x2))
	    x2 = x2 %>% rownames_to_column('row')
	    x2$row = factor(x2$row, levels=x2$row)
	    x2 = x2 %>% gather(col, size, -row)
	    x2$col = factor(x2$col, levels=cnames)	    
	}
	colnames(x2)[colnames(x2) == 'size'] = 'size2'
	if(is.null(x1)){x1 = x2} else {x1 = merge(x1, x2, by=c('row', 'col'), all=T)}
        x2 = as.data.table(x2)	
    }
    if(!is.null(y2)){
        if(! all(c('row', 'col', 'color') %in% colnames(y2))){
	    cnames = rev(colnames(y2))
	    y2 = y2 %>% rownames_to_column('row')
	    y2$row = factor(y2$row, levels=y2$row)
	    y2 = y2 %>% gather(col, color, -row)
	    y2$col = factor(y2$col, levels=cnames)	    
	}
	colnames(y2)[colnames(y2) == 'color'] = 'color2'
	if(is.null(y1)){y1 = y2} else {y1 = merge(y1, y2, by=c('row', 'col'), all=T)}
        y2 = as.data.table(y2)	
    }
        
    # merge data
    d = merge(x1, y1, by=c('row', 'col'), all=T)
    for(j in c('size', 'color', 'size2', 'color2')){if(! j %in% colnames(d)){d[,j] = NA}}
    d = as.data.table(d)
        
    # replace missing values
    d[,'size'][is.na(d[,'size'])] = replace_na
    d[,'size2'][is.na(d[,'size2'])] = replace_na    
    
    # fix plot order
    if(big_first == TRUE & all(c('size', 'size2') %in% colnames(d))){
        d[, order := size <= size2]
	d[order == TRUE, c('size', 'color', 'size2', 'color2') := .(size2, color2, size, color)]
    }
    else{
        
    }
    d = as.data.table(d)
        
    if(reorder != FALSE){
        if(! reorder %in% c('size', 'color', 'size2', 'color2')){stop('reorder != FALSE,size,color,size2,color2 - must specify attribute to use')}
	u = d[,.SD[which.max(abs(color)),row],col]
	Colv = u[order(factor(u[,V1], levels=Rowv)),col]
    }
    
    if(!is.null(Rowv)){
	d = d[d$row %in% Rowv]
        d$row = factor(d$row, levels=Rowv)
    }
    if(!is.null(Colv)){
	d = d[d$col %in% Colv,]
        d$col = factor(d$col, levels=rev(Colv))
    }
    
    # get palette
    if(length(pal) == 1){pal = brewer.pal(7, pal)}
    ci = c(d$color, d$color2)    
    if(is.null(vmin)){vmin = min(ci, na.rm=T)}
    if(is.null(vmax)){vmax = max(ci, na.rm=T)}
    d$color = qtrim(d$color, vmin=vmin, vmax=vmax, qmin=qmin, qmax=qmax)
    d$color2 = qtrim(d$color2, vmin=vmin, vmax=vmax, qmin=qmin, qmax=qmax)
    border1 = ifelse(d$order, '#000000', '#999999')
    border2 = ifelse(d$order, '#999999', '#000000')
    d$sig = 'P > .05'
    
    # add p-values
    if(!is.null(pvals)){
        pvals$row = make.names(pvals$row)
	pvals$col = make.names(pvals$col)
        setkeyv(d, c('row', 'col'))
        d[pvals[,.(row, col)], 'sig'] = 'P <= .05'
	d$sig = factor(d$sig, levels=c('P <= .05', 'P > .05'))
	p = ggplot(d) + geom_point(aes(x=col, y=row, size=size, fill=as.numeric(color), color=sig), pch=21, stroke=.25) + scale_color_manual('', values=c('black', 'grey'))
    } else {
        p = ggplot(d) + geom_point(aes(x=col, y=row, size=size, fill=as.numeric(color)), pch=21, color=border1, stroke=.25)
    }
    
    p = p + theme_cowplot() + xlab(xlab) + ylab(ylab) +
	scale_size_area(size_title, max_size=max_size) +
	theme(panel.grid.major=element_line(colour='black'), axis.text.x=element_text(angle=-45, hjust=0, vjust=.5)) + panel_border() + background_grid(major='xy') 
    
    if(any(d$size2 != 0)){p = p + geom_point(aes(x=col, y=row, size=size2, fill=as.numeric(color2)), pch=21, color=border2, stroke=.25)}
    
    if(symm == TRUE){
        values = c(vmin, 0, vmax)
	p = p + scale_fill_gradientn(fill_title, colours=pal, values=scales::rescale(values), na.value=na.value, limits=c(vmin, vmax))
    } else {
        p = p + scale_fill_gradientn(fill_title, colours=pal, na.value=na.value, limits=c(vmin, vmax))
    }
    if(coord_flip == TRUE){print('flip')
        p = p + coord_flip() + scale_x_discrete(limits=rev(levels(d$row))) + scale_y_discrete(limits=rev(levels(d$col)))
    } else {
        p = p + scale_x_discrete(limits=levels(d$col)) + scale_y_discrete(limits=levels(d$row))
    } 
    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }
        
    p
}


plot_violin = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL,
                       regex=NULL, files=NULL, file.regex=NULL, file.cols=NULL, top=NULL, type='mean', group_by=NULL, color_by=NULL, pt.size=.25,
	               do.facet=FALSE, facet_by=NULL, facet_genes=NULL, facet_formula=NULL, facet_scales='free_y', qmin=0, qmax=1, vmin=NULL, vmax=NULL, do.scale=FALSE, num_col='auto', resort=NULL,
                       ident=TRUE, cells.use=NULL, do.title=TRUE, do.legend=FALSE, xlab='', ylab='log2(TP10K+1)', out=NULL, nrow=1.5, ncol=1.5, legend.title='Group',
		       coord_flip=FALSE, alpha=1, order=NULL, font_size=10, do.pvals=FALSE, combine_genes='scale2', ret.data=FALSE,
		       mar=unit(c(6,6,6,6),'pt')){
    
    # Initialize plots
    ps = list()
    
    # Get facet formula
    if(is.null(facet_formula)){
    if(do.facet == FALSE){
        facet_formula = ifelse(is.null(facet_by), '~ .', 'Facet ~ .')
    } else {
        do.title = FALSE
        facet_formula = ifelse(is.null(facet_by), 'Feature ~ .', 'Feature ~ Facet + .')
    }
    } else {facet_formula = as.formula(facet_formula)}
    
    # Fix input arguments
    if(is.null(group_by)){group_by = seur@active.ident}
    if(is.null(color_by)){color_by = group_by}
    if(is.null(facet_by)){facet_by = rep('', ncol(seur@assays[['RNA']]@data))}
    if(is.null(cells.use)){cells.use = colnames(seur@assays[['RNA']]@data)}
    
    # Plot data
    d = data.frame(Group=as.factor(group_by), Color=as.factor(color_by), Facet=as.factor(facet_by), row.names=colnames(seur@assays[['RNA']]@data))
    d = d[cells.use,,drop=F]
    if(coord_flip == TRUE){d$Group = factor(d$Group, levels=rev(levels(d$Group)))}
    
    # Cell scores
    scores = score_cells(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, files=files, top=top, file.regex=file.regex, file.cols=file.cols, cells.use=cells.use,
                         combine_genes=combine_genes)
    scores = scores[cells.use,,drop=F]
    names = colnames(scores)
    d = cbind.data.frame(d, scores)
    
    # Fix NAs
    d = d[!is.na(d$Facet),,drop=F]
    
    # Scale data
    j = sapply(d, function(a) !is.factor(a))
    d[,j] = apply(d[,j,drop=F], 2, function(a) qtrim(a, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax))
    
    # Facet data
    if(do.facet == TRUE){
	d = gather(d, Feature, Value, -Group, -Color, -Facet)
	d$Feature = factor(d$Feature, levels=names)
    }
    
    if(!is.null(facet_genes)){
        d$Facet = facet_genes[d$Feature]
    }
    
    # Calculate p-values
    d$Group = droplevels(d$Group)
    if(do.pvals == TRUE){
    pvals = do.call(rbind, lapply(levels(d$Feature), function(fi) {do.call(rbind, 
        lapply(levels(d$Group), function(gi) {do.call(rbind, 
	    lapply(setdiff(levels(d$Facet), levels(d$Facet)[[1]]), function(hi){
                u = d[(d$Feature == fi & d$Group == gi) & d$Facet == levels(d$Facet)[[1]], 'Value']
	        v = d[(d$Feature == fi & d$Group == gi) & d$Facet == hi, 'Value']
                c(fi, gi, hi, max(v), tryCatch({wilcox.test(u,v,use='pairwise.complete.obs')$p.value}, error=function(e){1}))
	    })
	)})
    )}))
    
    pvals = data.frame(pvals, stringsAsFactors=F)
    colnames(pvals) = c('Feature', 'Group', 'Facet', 'Value', 'Pval')
    pvals$Value = as.numeric(pvals$Value)
    pvals$Pval = p.adjust(as.numeric(pvals$Pval), 'fdr')
    pvals$Pval[is.na(pvals$Pval)] = 1
    pvals$Label = ifelse(pvals$Pval <= 1e-3, '***', ifelse(pvals$Pval <= 1e-2, '**', ifelse(pvals$Pval <= .05, '*', '')))
    pvals$Color = pvals$Group
    pvals$Facet = factor(pvals$Facet, levels=levels(d$Facet))
    }
    
    # Violin plots
    for(col in setdiff(colnames(d), c('Group', 'Color', 'Facet', 'Value'))){ 
        
        # Convert to numeric?
	if(!is.factor(d[,col]) & !is.character(d[,col])){d[,col] = as.numeric(d[,col])}
	
        # Facet if necessary
        if(do.facet == TRUE){
	    if(!is.null(resort)){
	        i = d$Feature == levels(d$Feature)[[1]]
	        group_order = names(sort(tapply(d$Value[i], d$Group[i], resort)))
		d$Group = factor(d$Group, levels=group_order)
	    }
	    p = ggplot(data=d, aes_string(x='Group', y='Value', fill='Color'))
	} else {
	    if(!is.null(resort)){
	        group_order = names(sort(tapply(d$Group, d[,col], resort)))
		d$Group = factor(d$Group, levels=group_order)
	    }
	    p = ggplot(data=d, aes_string(x='Group', y=col, fill='Color'))
	}
	
	# Make violin plot
	p = p +
	    #geom_point(position=position_jitterdodge(dodge.width=0.6, jitter.width=1), size=pt.size, show.legend=F) +
	    geom_quasirandom(position=position_dodge(), size=pt.size, show.legend=F, method='pseudorandom') +
	    geom_violin(scale='width', alpha=alpha, size=.25) +
	    scale_fill_manual(values=set.colors) + theme_cowplot(font_size=font_size) +
	    xlab(xlab) + ylab(ylab) + labs(fill=legend.title) +
	    stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom='crossbar', width=.5, show.legend=F, size=.25) +
	    theme(axis.text.x = element_text(angle = -45, hjust = 0))
	
	if(do.title == TRUE){
	    p = p + ggtitle(col)
	}

	if(do.legend == FALSE){
	    p = p + guides(fill=FALSE)
	}

	if(coord_flip == TRUE){
	    p = p + coord_flip()
	}

	if(facet_formula != '~ .'){
	    p = p + facet_grid(as.formula(facet_formula), scales=facet_scales)
	}

	if(do.pvals == TRUE){
	    p = p + geom_text(data=pvals, aes(x=Group, y=Value + .01*max(Value), label=Label))
	}
	    
	ps[[col]] = p
    }

    if(num_col == 'auto'){
        num_col = ceiling(sqrt(length(ps)))
    }
    ps = make_compact(ps, num_col=num_col)
    p = plot_grid(plotlist=ps, ncol=num_col)
    p = p + theme(plot.margin=mar)
    if(!is.null(out)){
        save_plot(file=out, p, nrow=nrow, ncol=ncol)
    }
    p
}


matrix_barplot = function(data, group_by=NULL, pvals=NULL, xlab='', ylab='Frequency', value='mean', error='se', legend.title='Groups', colors='Paired', pos='dodge', border=NA,
                          out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE, sig_only=F, do.facet=F){
    
    # Plot barplot of [M x N] matrix
    # x-axis = matrix columns (e.g. cell types)
    # y-axis = matrix values (e.g. frequencies)
    # fill = matrix rows (e.g. samples) or groups (e.g. conditions)
    
    # Arguments:
    # group.by = the group of each row
    # pvals = [G x N] matrix of p-values for each group and column
    # error = sd, se, or none
    
    # Groups (default = rows)
    if(is.null(group_by)){group_by = rownames(data)}
    if(nlevels(group_by) == 0){group_by = as.factor(group_by)}
    
    # Select significant comparisons
    if(sig_only == TRUE){
        j = apply(pvals, 2, min) <= .05
	if(sum(j) == 0){return(NULL)}
	data = data[,j,drop=F]
	pvals = pvals[,j,drop=F]
    }
        
    # Construct input data
    names = colnames(data)
    data = data.frame(group=group_by, data)
    group_levels = levels(group_by)
    colnames(data)[2:ncol(data)] = names
    data = as.data.table(gather_(data, 'x', 'y', setdiff(colnames(data), 'group')))
    
    # Value function
    if(value == 'mean'){vf = mean} else if(value == 'median'){vf = median} else {stop()}
    
    # Error function
    se = function(x, na.rm=T){sd(x, na.rm=na.rm)/sqrt(length(x))}    
    if(error == 'sd'){ef = sd} else if(error == 'se'){ef = se} else {ef = function(x, ...){0}}
    
    # Estimate error bars
    data = data[,.(u=vf(y, na.rm=T), s=ef(y, na.rm=T)),.(group, x)]
    
    # Add p-values 1
    if(!is.null(pvals)){
        pvals = as.data.frame(pvals) %>% rownames_to_column('group') %>% gather(x, pval, -group) %>% as.data.table()
	setkeyv(data, c('x', 'group'))
	setkeyv(pvals, c('x', 'group'))
	data = merge(data, pvals, all=T)
	data$lab1 = ifelse(data$pval <= .001, '**', ifelse(data$pval <= .05, '*', ''))
    }

    if(coord_flip == TRUE){names = rev(names); group_levels=rev(group_levels)}
    data$x = factor(data$x, levels=names)    
    data$group = factor(data$group, levels=group_levels)
    
    # Get colors
    if(length(colors) == 1){colors = set.colors[1:length(group_levels)]}
    
    # Plot data
    if(pos == 'stack'){
        p = ggplot(data) + geom_bar(aes(x=x, y=u, fill=group), colour=border, size=.25, stat='identity')
	if(error %in% c('sd', 'se')){p = p + geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=group), stat='identity', width=.25)}
    } else {
        pos = position_dodge(.9)
        p = ggplot(data) + geom_bar(aes(x=x, y=u, fill=group), colour=border, size=.25, stat='identity', position=pos)
	if(error %in% c('sd', 'se')){p = p + geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=group), stat='identity', position=pos, width=.25)}
    }

    p = p + 
        scale_fill_manual(values=colors, name=legend.title) + xlab(xlab) + ylab(ylab) +
        scale_color_manual('', values=c('#000000', '#999999', '#cccccc'), guide='none')
    
    # Facet wrap
    if(do.facet == TRUE){
        p = p + facet_grid(group ~ ., scales='free')
    }

    dy = max(data$u + data$s, na.rm=T)*.01
    if(coord_flip == FALSE){
        p = p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
	if(!is.null(pvals)){p = p + geom_text(aes(x=x, y=u+s+dy, label=lab1, group=group), hjust='center', vjust=0, size=5, angle=0, position=pos)}
    } else {
        p = p + coord_flip()
	if(!is.null(pvals)){p = p + geom_text(aes(x=x, y=u+s+dy, label=lab1, group=group), hjust='center', vjust=1, size=5, angle=90, position=pos)}
    }
    
    # Save plot
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    p
}


simple_scatter = function(x, y, lab=NA, sig=FALSE, col=NULL, col.title='', size=NULL, size.title='', lab.use=NULL, lab.sig=FALSE, lab.near=0,  lab.n=0, lab.g=0, groups=NULL, lab.size=4, lab.type='up',
                          palette=NULL, legend.fontsize=10, border=F, edges=NULL, na.value='#cccccc',
                          xlab=NULL, ylab=NULL, out=NULL, nrow=1, ncol=1, min_size=1, max_size=3, xskip=c(0,0), yskip=c(0,0), xlim=NULL, ylim=NULL, alpha=1, unlab.grey=FALSE, auto_color='multi',
			  xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, xmin.n=1, xmax.n=1, ymin.n=1, ymax.n=1){
    
    if(is.null(col)){col = rep('', length(x))}
    if(is.null(sig)){sig = rep(NA, length(x))}
    if(is.null(size)){size = rep(1, length(x))}

    if(!is.null(edges)){colnames(edges) = c('x', 'y', 'xend', 'yend')}

    # adjust range
    xmin = max(na.omit(sort(x))[[xmin.n]], xmin)
    xmax = min(na.omit(sort(x, decreasing=T))[[xmax.n]], xmax)
    ymin = max(na.omit(sort(y))[[ymin.n]], ymin)
    ymax = min(na.omit(sort(y, decreasing=T))[[ymax.n]], ymax)
    x[x < xmin] = xmin
    x[x > xmax] = xmax
    y[y < ymin] = ymin
    y[y > ymax] = ymax
    
    d = data.frame(x=x, y=y, lab=lab, col=col, size=size, flag='', sig=sig, stringsAsFactors=FALSE)
    i.lab = !is.na(d$lab)
    di = d[i.lab,]
    
    if(is.null(xlab)){xlab = deparse(substitute(x))}
    if(is.null(ylab)){ylab = deparse(substitute(y))}
    
    if(lab.n > 0 | lab.g > 0){

        i = c()
	
	if('up' %in% lab.type | 'down' %in% lab.type){

            # get breaks
	    if(!is.null(groups)){groups.use = groups} else {
	        breaks = seq(from=min(di$x, na.rm=T), to=max(di$x, na.rm=T), length.out=min(20, lab.n))
	        groups.use = cut(di$x, breaks=breaks, include.lowest=TRUE)
		groups.use[xskip[[1]] < di$x & di$x < xskip[[2]]] = NA
	    }
	    
	    # get cells
	    if('up' %in% lab.type){
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=di$y, total_cells=lab.n)))
		i = c(i, order(-1*di$y)[1:lab.g])
	    }
	    if('down' %in% lab.type){
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=-1*di$y, total_cells=lab.n)))
		i = c(i, order(di$y)[1:lab.g])
	    }
	}
	
	if('right' %in% lab.type | 'left' %in% lab.type){

	    # get breaks
	    if(!is.null(groups)){groups.use = groups} else {
	        breaks = seq(from=min(di$y, na.rm=T), to=max(di$y, na.rm=T), length.out=min(20, lab.n))
	        groups.use = cut(di$y, breaks=breaks, include.lowest=TRUE)
		groups.use[yskip[[1]] < di$y & di$y < yskip[[2]]] = NA
	    }
	    
	    # get cells
	    if('right' %in% lab.type){
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=di$x, total_cells=lab.n)))
		i = c(i, order(-1*di$x)[1:lab.g])
	    }
	    
	    if('left' %in% lab.type){
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=-1*di$x, total_cells=lab.n)))
		i = c(i, order(di$x)[1:lab.g])
	    }
	}
	di[unique(i), 'flag'] = 'lab.n'
    }
    d[i.lab,] = di
    
    if(!is.null(lab.use)){
        
        # label neighbors
	if(lab.near > 0){
	    u = as.matrix(d[, c('x','y')])
	    v = as.matrix(d[lab %in% lab.use, c('x','y')])
	    i = unique(sort(unlist(apply(dist.matrix(u,v,skip.missing=TRUE,method='euclidean'), 2, order)[1:(lab.near + 1),])))
	    d$flag[i] = 'lab.near'
	}
	
        # label points
        d$flag[d$lab %in% lab.use] = 'lab.use'
    }
    
    if(lab.sig == TRUE){d$flag[d$sig == TRUE] = 'lab.use'; d$sig = FALSE}
        
    d$lab[d$flag == ''] = ''
    d = d[order(d$flag),]
    d$flag = factor(d$flag, levels=c('', 'lab.n', 'lab.use', 'lab.near'), ordered=T)

    if(auto_color == 'none'){d$col = ''}
    if(auto_color == 'bw' & all(col == '')){d$col = ifelse(d$flag == '', '', 'lab.n')}
    if(auto_color == 'multi' & all(col == '')){d$col = d$flag}
        
    # plot labeled points last
    d = d[order(d$flag),]
    i = is.na(d$col) | d$col == ''
    d = rbind(d[i,], d[!i,])
    
    if(unlab.grey == TRUE){d[d$lab == '',]$col = ''}

    p = ggplot(d)

    if(!is.null(edges)){p = p + geom_segment(data=edges, aes(x=x,y=y,xend=xend,yend=yend), color='black')}
    
    if(border == FALSE){
        p = p + geom_point(aes(x=x, y=y, col=col, size=size), alpha=alpha)
    } else {
        p = p + geom_point(aes(x=x, y=y, fill=col, size=size, stroke=ifelse(sig == TRUE, stroke, 0)), alpha=alpha, pch=21, colour=border)
    }
    p = p + geom_text_repel(data=d[(d$lab != '') | (1:nrow(d) %in% sample(1:nrow(d), min(nrow(d), 1000))),], aes(x=x, y=y, label=lab), size=lab.size, segment.color='grey') +
	theme_cowplot() + xlab(xlab) + ylab(ylab) + theme(legend.title=element_text(size=legend.fontsize), legend.text=element_text(size=legend.fontsize))
        
    if(all(size == 1)){p = p + scale_size(guide = 'none', range=c(min_size, max_size))} else {p = p + scale_size(name=size.title, range=c(min_size, max_size))}
    
    if(!is.null(xlim)){p = p + xlim(xlim)}
    if(!is.null(ylim)){p = p + ylim(ylim)}
    
    if(border == TRUE){
        if(!is.null(palette)){
            p = p + scale_fill_manual(name=col.title, values=palette, na.value=na.value, breaks=levels(as.factor(d$col)))
        } else {
            if(!is.numeric(d$col)){
	        p = p + scale_fill_manual(name=col.title, values=c('lightgrey', 'black', 'red', 'pink')) + theme(legend.position='none')
	    } else {
	        p = p + scale_fill_gradientn(name=col.title, colors=material.heat(100), na.value=na.value)
	    }
        }    
    } else {
        if(!is.null(palette)){
            p = p + scale_color_manual(name=col.title, values=palette, na.value=na.value, breaks=levels(as.factor(d$col)))
        } else {
            if(!is.numeric(d$col)){
	        p = p + scale_color_manual(name=col.title, values=c('lightgrey', 'black', 'red', 'pink')) + theme(legend.position='none')
	    } else {
	        p = p + scale_color_gradientn(name=col.title, colors=material.heat(100), na.value=na.value)
	    }
        }
    }
    
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}

    return(p)
}
