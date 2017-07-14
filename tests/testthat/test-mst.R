context('mst')

test_that('Building a minimum spanning tree', {
	adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE)
	expect_true(!is.null(adjl))
	adjl2 <- contract_mst(adjl, n_fold = 2)
	expect_true(!is.null(adjl2))
	ret <- gen_progindex(adjl, snap_start = 21)
	expect_true(!is.null(ret))
	ret2 <- gen_annotation(ret, snap_start = 21)
	expect_true(!is.null(ret2))
	plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', title = "CAMPARI WRAPPER - MST", return_plot = TRUE)
	expect_true(!is.null(plottt))
	plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', title = "CAMPARI WRAPPER - MST", return_plot = TRUE, only_timeline = TRUE)
	expect_true(!is.null(plottt))
	mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = TRUE, 
	             normalize_d = FALSE, clu_radius = 10, clu_hardcut = 100, min_span_tree = TRUE)
	expect_true(file.exists("MST_DUMPLING.nc"))
	ret <- gen_progindex(nsnaps = 100, read_from_netcdf = T, snap_start = 21)
	# expect_failure(gen_progindex(nsnaps = 10000, read_from_netcdf = T, snap_start = 21))
	expect_true(!is.null(ret))
	ret2 <- gen_annotation(ret, snap_start = 21, local_cut_width = 10)
	expect_true(!is.null(ret2))
	plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', timeline = T, title = "CAMPARI WRAPPER - MST", 
	                        return_plot = T, ann_trace = c(rep(2,50), rep(1,50)), timeline_proportion = 1.1, 
	                        rescaling_ann_col=TRUE, annotate_snap_dist= TRUE, horiz_lines_on_timeline= c(10,20,30), 
	                        reorder_annotation=TRUE,reorder_horizline_on_timeline=TRUE, points_on_timeline=c(2,10))
	expect_true(!is.null(plottt))
	
})

test_that('Keyfile handling', {
  a <- keywords_from_keyfile(key_file_input='KEYWORD 123', return_table=TRUE, return_string_of_arguments=FALSE,
                        keyword_list_first=TRUE, key_file_is_keywords=TRUE)
  expect_true(!is.null(a))
  a <- keywords_from_keyfile(key_file_input='KEYWORD 123', return_table=FALSE, return_string_of_arguments=TRUE,
                        keyword_list_first=FALSE, key_file_is_keywords=TRUE)
  expect_true(!is.null(a))
})