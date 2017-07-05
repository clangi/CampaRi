context('mst')

test_that('Building a minimum spanning tree', {
	adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE)
	expect_true(!is.null(adjl))
})
