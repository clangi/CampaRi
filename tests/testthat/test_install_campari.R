context('install_campari')

test_that('Test campari original installation from inst directory', {
  expect_error(install_campari(), NA)
  # expect_error(install_campari(install_ncminer = T), NA)
  # expect_error(install_campari(install_threads = T, silent_built = F), NA)
  # expect_error(install_campari(install_threads = T, install_ncminer = T, silent_built = F), NA)
  # expect_error(install_campari(install_threads = T, install_mpi = T, silent_built = F), NA)
  # expect_error(install_campari(install_threads = T, install_mpi = T, silent_built = F), NA)
  expect_error(install_campari(install_ncminer = T, install_threads = T, install_mpi = T, silent_built = F))
  expect_error(install_campari(installation_location = 'sadsdasdsadsadsadsadsadsafdfsdf', install_ncminer = T, install_threads = T, install_mpi = T, silent_built = F))
})