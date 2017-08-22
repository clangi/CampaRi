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
# ann <- matrix(sample(c(rep(1,100), rep(2,100))), nrow = 2, ncol = 100)
ann <- matrix(sample(c(1, 2, 3), size = 200, replace = T), nrow = 2, ncol = 100)
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

# legend stuff
expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T,
                          rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"),
                          annotation_type = 'continuous', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))
expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T,
                          rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"),
                          annotation_type = 'continuous', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))
expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
              rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000", "#C70039"),
              annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))

expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
              rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"),
              annotation_type = 'continuous', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))
expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
              rescaling_ann_col=F, annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))  
expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
              rescaling_ann_col=T, annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))  

expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
              rescaling_ann_col=F, annotation_type = 'continuous', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))

# can I do it with only charachters?
char_annotation <- factor(ann[1,])
levels(char_annotation) <- c("#b47b00", "#000000", "#C70039")
char_annotation <- as.character(char_annotation)
expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = char_annotation, timeline = T,
              rescaling_ann_col=T, annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))

# what is winning? rescaling_ann_col, charachter annotation or palette? # palette over charachter over rescaling
char_annotation <- factor(ann[1,])
levels(char_annotation) <- c("#b47b00", "#02AB1E", "#C70039")
char_annotation <- as.character(char_annotation)
expect_that(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = char_annotation, timeline = T,
              rescaling_ann_col=T,specific_palette_annotation= c("#b47b00", "#000000", "#C70039"),
              annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), not(throws_error()))
})