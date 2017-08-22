context('SAPPHIRE plot')

test_that('Various plotting options', {
adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE)
ret <- gen_progindex(adjl, snap_start = 21)
ret2 <- gen_annotation(ret, snap_start = 21)


plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', timeline = T, title = "CAMPARI WRAPPER - MST", 
                        return_plot = T, ann_trace = c(rep(2,50), rep(1,50)), timeline_proportion = 1.1, 
                        rescaling_ann_col=TRUE, annotate_snap_dist= TRUE, horiz_lines_on_timeline= c(10,20,30), 
                        reorder_annotation=TRUE,reorder_horizline_on_timeline=TRUE, points_on_timeline=c(2,10))
expect_true(!is.null(plottt))

# Classic plots 
# ------------------------------
expect_that(sapphire_plot('REPIX_000000000021.dat'), not(throws_error())) # baseline and local cut
expect_that(sapphire_plot('REPIX_000000000021.dat', timeline = TRUE), not(throws_error())) # with the timeline at the bottom # the dot size can be reduced internally
expect_that(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE), not(throws_error())) # without the local cut

# Annotation preprocessing
# ------------------------------
ann <- matrix(sample(c(rep(1,100), rep(2,100))), nrow = 2, ncol = 100)
# lets plot it!!!
expect_that(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = ann, reorder_annotation = T), not(throws_error())) # with ANNOTATION!!
expect_that(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = ann, annotate_snap_dist = T), not(throws_error())) # with distance values

# the advanced option -> adding layers to the gg object
# ------------------------------
gg <- sapphire_plot('REPIX_000000000021.dat', return_plot = TRUE, local_cut = FALSE, ann_trace = ann) # without the local cut
expect_true(!is.null(gg))

gg <- sapphire_plot('REPIX_000000000021.dat', return_plot = TRUE, local_cut = FALSE, ann_trace = ann, 
                    rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"))  # ADD a specific palette!
expect_true(!is.null(gg))

# add theme classic and text
gg + 
  ggplot2::theme_classic() + # changing the theme
  ggplot2::annotate("text", x = rep(100, 3), y = c(6.55, 6.9, 7.28), label = c("CCCH", "HCCC", "CCCC")) # annotate some text

# Showing only the temporal annotation
expect_that(sapphire_plot('REPIX_000000000021.dat', only_timeline = T), not(throws_error()))
expect_that(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, horiz_lines_on_timeline = c(20,40,60,80)), not(throws_error()))
expect_that(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, horiz_lines_on_timeline = c(20,40,60,80), horiz_colored_areas = c(1,1,2,1,2)), not(throws_error()))
expect_that(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, sub_sampling_factor = 3), not(throws_error())) # subsampling the timeline of a factor 3
expect_that(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, points_on_timeline = c(20,40,60,80)), not(throws_error())) # green points on the timeline
})












# setwd('../../CampaRi/vignettes/to_d/')




sapphire_plot('PROGIDX_000000000022.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T,
rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"),
annotation_type = 'continuous', legend_title = "Nalegenda", legend_labels = c('marimba', 'surrioga', 'katanga', 'supergyoza'))  # ADD a specific palette!

sapphire_plot('PROGIDX_000000000022.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T,
rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"),
annotation_type = 'discrete', legend_title = "Nalegenda", legend_labels = c('marimba', 'katanga', 'supergyoza'))  # ADD a specific palette!
traceback()


xa <- seq(1,20000)
ya <- rep(0, length(xa))
anna <- as.factor(sample(c(1,2), size = length(xa), replace = T))
ggplot() + geom_segment(aes(x = xa, y = ya, xend = xa, yend = ya + 1,
                            col = ann), 
                        size = 0.1) +
  scale_color_manual(name = "tit", values = c("#b47b00", "#000000"),
                     labels = c('bels1', 'bels2'))






