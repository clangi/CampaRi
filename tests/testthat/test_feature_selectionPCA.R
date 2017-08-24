context('feature selection')


test_that('feature selection and its plotting', {
  data_f <- matrix(rnorm(n = 1000), nrow = 100, ncol = 10)
  expect_true(!is.null(data_f))
  slected_elements <- select_features(data_f, feature_selection = 'pca') # automatically selected the first 2 components
  expect_true(!is.null(slected_elements))
  slected_elements_robust <- select_features(data_f, feature_selection = 'pca', pca_method = 'robust') # projection pursuit 2 components
  expect_true(!is.null(slected_elements_robust))
  
  
  # plotting
  expect_that(select_features(data_f, feature_selection = 'pca', plotit = T), not(throws_error()))
  expect_that(select_features(data_f, feature_selection = 'pca', plotit = T, 
                  cluster_vector = sample(c(1,2,3), size = nrow(data_f), replace = T)), not(throws_error()))
  
  # let's use a better clustering definition:
  clu_vector <- sample(c(1,2), size = nrow(data_f), replace = T)
  expect_true(!is.null(clu_vector))
  # the selected_elements* objects are data_frames with 2 columns (selected variables)
  # If I want to select more variables:
  slected_elements <- select_features(data_f, feature_selection = 'pca', n_princ_comp = 4)
  expect_true(!is.null(slected_elements))
  expect_that(select_features(data_f, feature_selection = 'pca', plotit = T, cluster_vector = clu_vector), not(throws_error()))
  
  expect_that(select_features(data_f, feature_selection = 'pca', plotit = T, frameit = T, cluster_vector = clu_vector), not(throws_error()))
  
  plot1 <- select_features(data_f, feature_selection = 'pca', plotit = T, frameit = T, return_plot = T, cluster_vector = clu_vector)
  expect_true(!is.null(plot1))
  expect_that(select_features(data_f, feature_selection = 'pca', plotit = T, frameit = F, plotly_it = T, points_size = 1.3, cluster_vector = clu_vector), 
              not(throws_error()))
  
  # adding the legend?
  expect_that(select_features(data_f, feature_selection = 'pca', 
                  plotit = T, frameit = T, plotly_it = F, points_size = 1.3, cluster_vector = clu_vector,
                  plot_legend = T), not(throws_error()))
  expect_that(select_features(data_f, feature_selection = 'pca', 
                  plotit = T, frameit = F, plotly_it = F, points_size = 1.3, cluster_vector = clu_vector,
                  plot_legend = T, specific_palette = c("#b47b00", "#000000")), not(throws_error()))
  expect_that(select_features(data_f, feature_selection = 'pca', 
                  plotit = T, frameit = T, points_size = 1.3, cluster_vector = clu_vector,
                  plot_legend = T), not(throws_error()))
  if(file.exists('selected_pca.out')) file.remove('selected_pca.out')
})
