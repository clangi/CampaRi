context('show_clustering_summary')

test_that('Test for showing the clustering summary', {
  file_fyc <- system.file("extdata", "NBU.fyc", package = "CampaRi")
  file_fyc2 <- system.file("extdata", "FYC_example.dat", package = "CampaRi")
  file_log <- system.file("extdata", "NBU.log", package = "CampaRi")
  
  # specific cluster option
  expect_error(show_clustering_summary(log_file = file_log, fyc_file = file_fyc, which_cluster = 10, return_centers = T), NA)
  expect_error(show_clustering_summary(log_file = file_log, fyc_file = file_fyc2, which_cluster = 10, return_centers = T))
  
  # first 11 clusters
  expect_error(show_clustering_summary(log_file = file_log, fyc_file = file_fyc, which_first_clusters = 11, return_centers = T), NA)
  expect_error(show_clustering_summary(log_file = file_log, fyc_file = file_fyc2, which_first_clusters = 11, return_centers = T))


  ahoh <- system(paste0("printf \"%s\n\" `head -n 1 ", file_fyc, " | sed -e 's/||/ /g' | sed -e 's/|//g' | sed -e 's/I 0/I0/g' | sed -e 's/C 0/C0/g' | \
                        sed -e 's/\\([A-Z]\\) OME/\\1_OME/g' | sed -e 's/\\([A-Z]\\) PHI/\\1_PHI/g' | sed -e 's/\\([A-Z]\\) PSI/\\1_PSI/g' | \
                        sed -e 's/\\([A-Z]\\) CHI/\\1_CHI/g' | \
                        sed -e 's/\\([A-Z]\\) NUC/\\1_NUC/g'` | awk -v rs=0 '{if ((NR > 1) && (length($1) > 5)) {rs = rs + 1; split($1, resname, \"_\"); \
                        printf(\"%5d : %10s %5d\\n\", NR, $1, rs)}; if ((NR > 1) && (length($1) <= 5)) {printf(\"%5d : %10s %5d\\n\", NR, resname[1]\"_\"$1, rs);}}'"),
                 intern = TRUE)
  expect_true(!is.null(ahoh))
})