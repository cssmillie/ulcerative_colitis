

num_cells_per_group = function(groups, total_cells=NULL, cells_per_group=NULL){

    # Calculate number of cells to downsample from each group
    # -------------------------------------------------------
    # Arguments:
    # - groups = vector of groups for each cell barcode
    # - total_cells = total number of cells to downsample
    # - cells_per_group = total number of cells to select per group
    
    num_cells = sort(table(groups))
    
    if(!is.null(cells_per_group)){
        num_cells[num_cells > cells_per_group] = cells_per_group
    } else {
        n = sort(table(groups))
	if(length(n) == 1){
	    num_cells = total_cells
	    names(num_cells) = names(n)
	} else {
	    u = c(0, cumsum(n)[1:(length(n)-1)])
	    i = (total_cells - u)/seq(length(n), 1, -1) < n
	    if(sum(i) > 0){
	        num_cells[i] = as.integer(ceiling((total_cells - sum(n[!i]))/sum(i)))
	    }
	}
    }
    
    num_cells
}

resample = function(x,...){if(length(x)==1) x else sample(x,...)} 

simple_downsample = function(cells, groups, ngene=NULL, total_cells=NULL, cells_per_group=NULL){

    # Downsample cells evenly across groups
    # -------------------------------------
    # Arguments:
    # - cells = vector of cell names
    # - groups = vector of cell groups
    # - ngene = cell quality (if provided, select highest quality cells from each group - otherwise, select randomly)
    # - total_cells = total number of cells to downsample
    # - cells_per_group = total number of cells to downsample from each group
        
    # Set ngene (default = 1)
    if(is.null(ngene)){
        ngene = structure(rep(1, length(cells)), names=cells)
    }
    if(is.null(names(ngene))){names(ngene) = cells}
    
    # Calculate group sizes
    groups = as.factor(groups)
    num_cells_per_group = num_cells_per_group(groups=groups, total_cells=total_cells, cells_per_group=cells_per_group)
        
    # Downsample cells within each group
    ds.cells = sapply(levels(groups), function(a){

        # Shuffle cells within group
        cells = resample(cells[groups == a])
	
	# Select by highest ngene
	cells[order(ngene[cells], decreasing=T)[1:num_cells_per_group[[a]]]]
    })
    ds.cells = as.character(na.omit(unname(unlist(ds.cells))))
    return(ds.cells)
}
