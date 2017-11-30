context('mst')

test_that('Building a minimum spanning tree', {
	expect_error(adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), distance_method = 1, dump_to_netcdf = FALSE), NA)
	expect_error(adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), distance_method = 1, dump_to_netcdf = TRUE), NA)
	expect_error(adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE), NA)
	expect_true(!is.null(adjl))
	expect_error(adjl2 <- contract_mst(adjl, n_fold = 2), NA)
	expect_true(!is.null(adjl2))
	expect_error(ret <- gen_progindex(adjl, snap_start = 21), NA)
	expect_true(!is.null(ret))
	expect_error(ret2 <- gen_annotation(ret, snap_start = 21), NA)
	expect_error(mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = TRUE, 
	                          normalize_d = FALSE, clu_radius = 10, clu_hardcut = 100, min_span_tree = TRUE), NA)
	expect_true(file.exists("MST_DUMPLING.nc"))
	expect_error(ret <- gen_progindex(nsnaps = 100, read_from_netcdf = T, snap_start = 21), NA)
	# expect_failure(gen_progindex(nsnaps = 10000, read_from_netcdf = T, snap_start = 21))
	expect_true(!is.null(ret))
	expect_error(ret2 <- gen_annotation(ret, snap_start = 21, local_cut_width = 10), NA)
	expect_true(!is.null(ret2))

	if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
	if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})

