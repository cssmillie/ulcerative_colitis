run_parallel = function(f, n.cores=1){    

    # Execute parallel function on cluster
    # Example:
    # g = run_parallel(foreach(i=1:10, .combine=cbind), %dopar% rnorm(100, i), n.cores=10)
    
    if(n.cores > 1){    
        
        # Register cluster
	cluster = makeCluster(n.cores, type='FORK', outfile='')
	registerDoParallel(cluster)
	
        # Run f in parallel
    	y = f
	
	# Register sequential
	stopCluster(cluster)
    	registerDoSEQ()

    } else {
        
	# Run f sequentially
        y = f
    }

    return(y)
}
