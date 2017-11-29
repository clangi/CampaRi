context('plot_traces')

test_that('Test simple traces plotting', {
  expect_error(plot_traces(matrix(rnorm(1000), nrow = 10, ncol = 100), return_plot = T,
                           sub_sampling_factor = NULL, title = "a", xlab = "b", ylab = "c", highlight_rows = c(1,2),
                           size_line = 1.3, nrows_to_plot = NULL), NA)
  expect_error(plot_traces(matrix(rnorm(1000), nrow = 10, ncol = 100), return_plot = F,
                           sub_sampling_factor = NULL, title = "a", xlab = "b", ylab = "c", highlight_rows = c(1,2),
                           size_line = 1.3, nrows_to_plot = NULL), NA)
  expect_error(plot_traces(matrix(rnorm(1000), nrow = 10, ncol = 100), return_plot = F,
                           sub_sampling_factor = NULL, title = "a", xlab = "b", ylab = "c", highlight_rows = c(1,2),
                           size_line = 1.3, nrows_to_plot = c(1,5)), NA) #nrows_to_plot does do nothing
  expect_error(plot_traces(matrix(rnorm(1000), nrow = 10, ncol = 100), return_plot = F,
                           sub_sampling_factor = 4, title = "a", xlab = "b", ylab = "c", highlight_rows = c(1,2),
                           size_line = 1.3, nrows_to_plot = c(1,5)), NA) 
  expect_warning(plot_traces(matrix(rnorm(1000), nrow = 10, ncol = 100), return_plot = F,
                           sub_sampling_factor = 4.2, title = "a", xlab = "b", ylab = "c", highlight_rows = c(1,2),
                           size_line = 1.3, nrows_to_plot = c(1,5))) 
})