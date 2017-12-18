context('SAPPHIRE plot')

test_that('Various plotting options', {
  # general setting and initialization
  adjl <- mst_from_trj(trj = matrix(rnorm(50000), nrow = 5000, ncol = 10), dump_to_netcdf = FALSE)
  ret <- gen_progindex(adjl, snap_start = 21)
  ret2 <- gen_annotation(ret, snap_start = 21)
  
  
  plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', timeline = T, title = "CAMPARI WRAPPER - MST", 
                          return_plot = T, ann_trace = c(rep(2,50), rep(1,50)), timeline_proportion = 1.1, 
                          rescaling_ann_col=TRUE, annotate_snap_dist= TRUE, horiz_lines_on_timeline= c(10,20,30), 
                          reorder_annotation=TRUE,reorder_horizline_on_timeline=TRUE, points_on_timeline=c(2,10))
  expect_true(!is.null(plottt))
  expect_error(sapphire_plot('REPIsadsasdadsasda12312382733569'))  
  
  
  # Classic plots 
  # ------------------------------
  expect_error(sapphire_plot('REPIX_000000000021.dat'), NA) # baseline and local cut
  expect_error(sapphire_plot('REPIX_000000000021.dat', ann_trace = T), NA) # baseline and local cut
  expect_error(sapphire_plot('REPIX_000000000021.dat', ann_trace = 7), NA) # baseline and local cut
  expect_error(sapphire_plot('REPIX_000000000021.dat', timeline = TRUE), NA) # with the timeline at the bottom # the dot size can be reduced internally
  expect_error(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE), NA) # without the local cut
  expect_error(sapphire_plot('REPIX_000000000021.dat', sub_sampling_factor = 10), NA) # w/subsampling
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_ann_trace = T)) # To check better
  expect_error(sapphire_plot('REPIX_000000000021.dat', title = 'ragnarok'), NA) # To check betterc("#b47b00", "#000000","#000000")
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = T), NA) # To check better
  
  # Annotation preprocessing
  # ------------------------------
  # with splitting in 2 parts -> ann_trace = T
  expect_warning(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = T, reorder_annotation = T)) 
  expect_error(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = T, reorder_annotation = F), NA) 
  expect_error(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = 5, reorder_annotation = F), NA) 
  
  # defining an annotation
  ann <- matrix(sample(c(1, 2, 3), size = 200, replace = T), nrow = 2, ncol = 5000)
  # lets plot it!!!
  expect_error(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = ann, reorder_annotation = T), NA) # with ANNOTATION!!
  expect_error(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = ann, reorder_annotation = F), NA) # with ANNOTATION!!
  
  # the advanced option -> adding layers to the gg object
  # ------------------------------
  gg <- sapphire_plot('REPIX_000000000021.dat', return_plot = TRUE, local_cut = FALSE, ann_trace = ann) # without the local cut
  expect_true(!is.null(gg))
  expect_error(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = ann, 
                      rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000")))  # ADD a specific palette!
  expect_error(sapphire_plot('REPIX_000000000021.dat', local_cut = FALSE, ann_trace = ann, 
                      rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000","#000000")), NA)  # ADD a specific palette!
  # add theme classic and text
  expect_error(gg + 
    ggplot2::theme_classic() + # changing the theme
    ggplot2::annotate("text", x = rep(100, 3), y = c(6.55, 6.9, 7.28), label = c("CCCH", "HCCC", "CCCC")), NA) # annotate some text
  
  # Showing only the temporal annotation
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, horiz_lines_on_timeline = c(20,40,60,80)), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, horiz_lines_on_timeline = c(20,40,60,80), horiz_colored_areas = c(1,1,2,1,2)), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, sub_sampling_factor = 3), NA) # subsampling the timeline of a factor 3
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, points_on_timeline = c(20,40,60,80)), NA) # green points on the timeline
  
  # testing temporal and geometric annotation combinations:
  # -------------------------------------------------------
  ann_timeline <- as.vector(matrix(sample(c(1, 2, 3), size = 200, replace = T), nrow = 1, ncol = 5000))
  # with geometric:
  expect_error(sapphire_plot('REPIX_000000000021.dat', timeline = T, ann_trace = ann, timeline_trace = ann_timeline), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', timeline = T, ann_trace = ann, uniform_color_timeline = T,
                             which_uniform_color_timeline = "blue"), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', timeline = T, ann_trace = ann, uniform_color_timeline = F,
                             which_uniform_color_timeline = "annotation", timeline_trace = ann_timeline), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', timeline = T, ann_trace = ann, uniform_color_timeline = F,
                             which_uniform_color_timeline = "annotation", timeline_trace = ann_timeline, timeline_annotation_type = "discrete"), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', timeline = T, ann_trace = ann, uniform_color_timeline = F,
                             which_uniform_color_timeline = "annotation", timeline_trace = ann_timeline, 
                             specific_palette_annotation= c("#b47b00", "#000000","#000000"), timeline_annotation_type = "continuous"), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', timeline = T, ann_trace = ann, uniform_color_timeline = F,
                             which_uniform_color_timeline = "annotation", timeline_trace = ann_timeline, 
                             specific_palette_annotation= c("#b47b00", "#000000","#000000"), timeline_annotation_type = "discrete",
                             specific_palette_timeline=c("#d7191c", "#ffffbf", "#2c7bb6")), NA)
  # only_timeline:
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, timeline_trace = ann_timeline), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, timeline_trace = ann_timeline, which_uniform_color_timeline = "annotation"), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, timeline_trace = ann_timeline, 
                             specific_palette_timeline = c("#d7191c", "#3cafb8", "#2c7bb6"), which_uniform_color_timeline = "annotation"), NA)
  
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, timeline_trace = ann_timeline, ann_trace = ann,
                             which_uniform_color_timeline = "annotation"), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, timeline_trace = ann_timeline, ann_trace = ann,
                             which_uniform_color_timeline = "annotation", timeline_annotation_type = "continuous", rescaling_ann_col = F,
                             specific_palette_timeline = c("#b47b00", "#000000","#000000")), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, timeline_trace = ann_timeline, ann_trace = ann,
                             which_uniform_color_timeline = "annotation", timeline_annotation_type = "continuous", rescaling_ann_col = F), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, timeline_trace = ann_timeline, ann_trace = ann,
                             which_uniform_color_timeline = "annotation", timeline_annotation_type = "continuous", rescaling_ann_col = F,
                             specific_palette_timeline = rev(c("#d7191c", "#3cafb8", "#2c7bb6"))), NA)
  
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, timeline_trace = ann_timeline, timeline_annotation_type = "discrete",
                             size_points_on_timeline = 1.5), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, which_uniform_color_timeline = "annotation", ann_trace = ann,
                             timeline_trace = ann_timeline, timeline_annotation_type = "continuous", rescaling_ann_col = F), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', only_timeline = T, which_uniform_color_timeline = "annotation", ann_trace = ann,
                            timeline_trace = ann_timeline, timeline_annotation_type = "continuous", rescaling_ann_col = F, 
                            specific_palette_timeline=c("#000000", "#3cafb8", "#2c7bb6"), size_points_on_timeline = 1.5), NA)
  
  # legend stuff
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T,
                            rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"),
                            annotation_type = 'continuous', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T,
                            rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"),
                            annotation_type = 'continuous', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
                rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000", "#C70039"),
                annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)
  
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
                rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000"),
                annotation_type = 'continuous', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
                rescaling_ann_col=F, annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)  
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
                rescaling_ann_col=T, annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)  
  
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann[1,], timeline = T,
                rescaling_ann_col=F, annotation_type = 'continuous', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)
  
  # can I do it with only charachters?
  char_annotation <- factor(ann[1,])
  levels(char_annotation) <- c("#b47b00", "#000000", "#C70039")
  char_annotation <- as.character(char_annotation)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = char_annotation, timeline = T,
                rescaling_ann_col=T, annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)
  
  # what is winning? rescaling_ann_col, charachter annotation or palette? # palette over charachter over rescaling
  char_annotation <- factor(ann[1,])
  levels(char_annotation) <- c("#b47b00", "#02AB1E", "#C70039")
  char_annotation <- as.character(char_annotation)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = char_annotation, timeline = T,
                rescaling_ann_col=T,specific_palette_annotation= c("#b47b00", "#000000", "#C70039"),
                annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = char_annotation, timeline = T, uniform_color_timeline=T,
                            rescaling_ann_col=T,specific_palette_annotation= c("#b47b00", "#000000", "#C70039"),
                            annotation_type = 'discrete', legend_title = "Pot en min", legend_labels = c('gauche+', 'gauche-', 'anti')), NA)
  
  # new stuff ------- for errors in subsampling factor
  ann <- matrix(sample(seq(-10, 10, length.out = 10000)), nrow = 2, ncol = 5000)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T, plot_legend = T, sub_sampling_factor = 1,
                            rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000", "#C70039")), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T, plot_legend = T, sub_sampling_factor = 7,
                            rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000", "#C70039")), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', return_plot = F, local_cut = T, ann_trace = ann, timeline = T, plot_legend = T, sub_sampling_factor = 5000,
                            rescaling_ann_col=F, specific_palette_annotation= c("#b47b00", "#000000", "#C70039")))
  
  # checking the snapdist
  ann <- matrix(sample(c(1, 2, 3), size = 200, replace = T), nrow = 2, ncol = 5000)
  expect_error(sapphire_plot('REPIX_000000000021.dat', annotate_snap_dist = T), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', ann_tr = T, annotate_snap_dist = T), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', ann_tr = T, snap_dist_annotation_height = 0.5), NA)
  
  # checking plotly
  expect_error(sapphire_plot('REPIX_000000000021.dat', annotate_snap_dist = T, use_plotly = T), NA)
  expect_error(sapphire_plot('REPIX_000000000021.dat', annotate_snap_dist = T, timeline = T, use_plotly = T, size_points_on_timeline = 0.1), NA)
  
  # checking the parabolic_subtraction
  expect_error(sapphire_plot('REPIX_000000000021.dat', ann_trace = ann, parabolic_subtraction = T, timeline = T, plot_legend = T), NA)  
  
  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})
  