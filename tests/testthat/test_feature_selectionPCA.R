context('feature selection')


test_that('feature selection and its plotting', {
  
  silent <- T
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)} 
  
  
  data_f <- matrix(rnorm(n = 1000), nrow = 100, ncol = 10)
  expect_error(slected_elements <- select_features(data_f, feature_selection = 'pca', silent = T), NA) # automatically selected the first 2 components
  expect_error(slected_elements_robust <- select_features(data_f, feature_selection = 'pca', pca_method = 'robust', silent = T), NA) # projection pursuit 2 components
  
  
  # plotting
  # expect_error(select_features(data_f, feature_selection = 'pca', plotit = T, silent = T), NA)
  a <- capture_output_lines(expect_error(select_features(data_f, feature_selection = 'pca', plotit = T, 
                  cluster_vector = sample(c(1,2,3), size = nrow(data_f), replace = T), silent = F), NA))
  
  # let's use a better clustering definition:
  clu_vector <- sample(c(1,2), size = nrow(data_f), replace = T)
  # the selected_elements* objects are data_frames with 2 columns (selected variables)
  # If I want to select more variables:
  expect_error(slected_elements <- select_features(data_f, feature_selection = 'pca', n_princ_comp = 4, silent = T), NA)
  expect_error(select_features(data_f, feature_selection = 'pca', plotit = T, cluster_vector = clu_vector, silent = T), NA)
  
  expect_error(select_features(data_f, feature_selection = 'pca', plotit = T, frameit = T, cluster_vector = clu_vector, silent = T), NA)
  
  expect_error(plot1 <- select_features(data_f, feature_selection = 'pca', plotit = T, frameit = T, return_plot = T, cluster_vector = clu_vector, silent = T), NA)
  
  # adding the legend?
  expect_error(select_features(data_f, feature_selection = 'pca', 
                  plotit = T, frameit = T, plotly_it = F, points_size = 1.3, cluster_vector = clu_vector,
                  plot_legend = T, silent = T), NA)
  expect_error(select_features(data_f, feature_selection = 'pca', 
                  plotit = T, frameit = F, plotly_it = F, points_size = 1.3, cluster_vector = clu_vector,
                  plot_legend = T, specific_palette = c("#b47b00", "#000000"), silent = T), NA)
  expect_error(select_features(data_f, feature_selection = 'pca', 
                  plotit = T, frameit = T, points_size = 1.3, cluster_vector = clu_vector,
                  plot_legend = T, silent = T), NA)
  if(file.exists('selected_pca.out')) file.remove('selected_pca.out')
})
