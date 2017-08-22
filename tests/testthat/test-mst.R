context('mst')

test_that('Building a minimum spanning tree', {
	adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE)
	expect_true(!is.null(adjl))
	adjl2 <- contract_mst(adjl, n_fold = 2)
	expect_true(!is.null(adjl2))
	ret <- gen_progindex(adjl, snap_start = 21)
	expect_true(!is.null(ret))
	ret2 <- gen_annotation(ret, snap_start = 21)
	mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = TRUE, 
	             normalize_d = FALSE, clu_radius = 10, clu_hardcut = 100, min_span_tree = TRUE)
	expect_true(file.exists("MST_DUMPLING.nc"))
	ret <- gen_progindex(nsnaps = 100, read_from_netcdf = T, snap_start = 21)
	# expect_failure(gen_progindex(nsnaps = 10000, read_from_netcdf = T, snap_start = 21))
	expect_true(!is.null(ret))
	ret2 <- gen_annotation(ret, snap_start = 21, local_cut_width = 10)
	expect_true(!is.null(ret2))

	
})

