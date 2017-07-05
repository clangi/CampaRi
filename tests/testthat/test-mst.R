context('mst')

test_that('Building a minimum spanning tree', {
	adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE)
	expect_true(!is.null(adjl))
	ret <- gen_progindex(nsnaps = 1000, snap_start = 21)
	expect_true(!is.null(ret))
	ret2 <- gen_annotation(ret, snap_start = 21)
	expect_true(!is.null(ret2))
	plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', title = "CAMPARI WRAPPER - MST", return_plot = TRUE)
	expect_true(!is.null(plottt))
})
